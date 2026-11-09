/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particles.h>

#include <exaStamp/analysis_particle/grains/grain_segmentation_algo.h>

#include <ptm_constants.h>

// Grain segmentation, pass (ii): host-side, sequential (same as OVITO's own engine -- neither the
// Node-Pair-Sampling dendrogram construction nor Dijkstra orphan adoption is GPU-parallelizable,
// each merge/settle step depends on prior state). Consumes compute_grain_bond_misorientation's own
// per-atom bond arrays (pass i, GPU-compatible) plus PTM's own struct_type/orientation fields, and
// writes a per-atom grain_id field (-1 default = "not part of any grain") plus a per-grain RGB color field looked up from the final
// grain id.
//
// Algorithm (see grain_segmentation_algo.cpp for the actual implementation, ported from OVITO's
// GrainSegmentationEngine1/2, crystalanalysis/modifier/grains):
//  1. Build a graph over every candidate crystalline bond (compute_grain_bond_misorientation's own
//     4-degree-gated disorientation), edge weight exp(-disorientation^2/3).
//  2. Node-Pair-Sampling: reciprocal-nearest-neighbor chain agglomerative clustering (mutual-nearest
//     -neighbor pairs contract first, union-by-degree), producing a dendrogram with each atom's own
//     running weighted quaternion sum accumulated incrementally at every merge.
//  3. Auto-threshold: robust (IRLS) log-log power-law fit of merge distance vs. merge size, cutting
//     at the largest "inlier" distance -- picks up the anomalous jump where within-grain noise gives
//     way to a genuine inter-grain merge. This is what gives robustness against thermal noise (the
//     motivating use case for choosing Node-Pair-Sampling + auto-threshold over a fixed cutoff).
//  4. Cut the dendrogram at that threshold, dissolve any surviving cluster under min_grain_atom_count
//     back to "no grain", assign contiguous ids largest-to-smallest, random per-grain HSV color
//     (fixed seed 1, matching OVITO exactly).
//  5. Optional orphan-atom adoption: multi-source Dijkstra over every recorded bond (not just
//     crystalline candidates), real bond length as edge cost -- assigns each still-unassigned atom
//     (PTM "Other"/disordered, or a dissolved too-small cluster) to whichever grain is geometrically
//     closest by cumulative real-space path length.
//
// NOT implemented (see doc/grain_segmentation_workflow.tex for the full comparison against
// OVITO): coherent-interface (twin/stacking-fault) reclassification; cross-rank MPI grain stitching
// (each rank clusters its own local owned+ghost view independently -- a grain spanning a rank
// boundary is NOT currently reconciled into one id).
namespace exaStamp
{
  using namespace exanb;

  struct GrainIdFieldFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const int32_t * const __restrict__ m_atom_grain_id = nullptr;
    const double * const __restrict__ m_grain_color = nullptr; // 3/grain, indexed grain_id-1
    inline void operator () ( size_t cell, unsigned int part, double& id_out, Vec3d& color_out ) const
    {
      const size_t i = m_cell_particle_offset[cell] + part;
      const int32_t g = m_atom_grain_id[i];
      id_out = static_cast<double>(g);
      color_out = ( g > 0 ) ? Vec3d{ m_grain_color[3*(g-1)+0], m_grain_color[3*(g-1)+1], m_grain_color[3*(g-1)+2] } : Vec3d{0.8,0.8,0.8};
    }
  };

  template<class GridT>
  class ComputeGrainClusters : public OperatorNode
  {
    ADD_SLOT( GridT   , grid , INPUT_OUTPUT );
    ADD_SLOT( std::string , struct_field , INPUT , std::string("ptm_type") , DocString{"Per-particle PTM structure-type field name (see ptm_fields)"} );
    ADD_SLOT( std::string , orient_field , INPUT , std::string("ptm_orientation") , DocString{"Per-particle PTM lattice-orientation rotation-tensor field name (see ptm_fields)"} );
    ADD_SLOT( onika::memory::CudaMMVector<uint64_t> , grain_bond_id            , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double>   , grain_bond_distance      , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double>   , grain_bond_disorientation , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<int>      , grain_bond_count         , INPUT , REQUIRED );
    ADD_SLOT( bool   , auto_threshold  , INPUT , true , DocString{"true (default, matches OVITO's own GraphClusteringAutomatic): auto-select the merge threshold via a robust log-log power-law fit -- the recommended, noise-robust choice for real/thermal-noise data. false: use manual_threshold_log directly (OVITO's own internal log-distance unit, NOT degrees -- see this file's own header comment)."} );
    ADD_SLOT( double , manual_threshold_log , INPUT , 0.0 , DocString{"Only used when auto_threshold=false -- OVITO's own internal log-distance unit for Node-Pair-Sampling mode, not a degree value"} );
    ADD_SLOT( long   , min_grain_atom_count , INPUT , 100 , DocString{"Clusters below this atom count are dissolved back to \"no grain\" -- same default as OVITO's own GrainSegmentationModifier"} );
    ADD_SLOT( bool   , orphan_adoption , INPUT , true , DocString{"Adopt PTM \"Other\"/disordered atoms (and dissolved too-small clusters) into the geometrically closest grain via multi-source Dijkstra over real bond lengths -- matches OVITO's own default (on)"} );
    ADD_SLOT( unsigned int , color_seed , INPUT , 1u , DocString{"RNG seed for per-grain random color assignment -- default 1 matches OVITO's own hardcoded seed exactly"} );
    ADD_SLOT( std::string , grain_id_field , INPUT , std::string("grain_id") , DocString{"Name of the resulting per-particle grain-id field (-1 = not part of any grain)"} );
    ADD_SLOT( std::string , grain_color_field , INPUT , std::string("grain_color") , DocString{"Name of the resulting per-particle RGB grain-color field"} );
    ADD_SLOT( long , n_grains , OUTPUT );
    ADD_SLOT( long , n_orphans_adopted , OUTPUT );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) { *n_grains=0; *n_orphans_adopted=0; return; }
      if( ! grid->has_allocated_field( field::mk_generic_real( *struct_field ) ) )
      {
        fatal_error() << "compute_grain_clusters: input field '" << *struct_field << "' does not exist (run compute_ptm + ptm_fields first)" << std::endl;
      }
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *orient_field ) ) )
      {
        fatal_error() << "compute_grain_clusters: input field '" << *orient_field << "' does not exist (run compute_ptm + ptm_fields first)" << std::endl;
      }

      const size_t n_particles = grid->number_of_particles();
      const size_t n_cells = grid->number_of_cells();
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      const int max_neighbors = static_cast<int>( grain_bond_count->size()>0 ? (grain_bond_id->size()/grain_bond_count->size()) : 0 );

      // Cluster OWNED atoms only. Ghost atoms are periodic/MPI images sharing the global id of a real
      // atom: as separate graph nodes they split boundary atoms in two (bond targets resolved to
      // whichever copy came last) and get counted in cluster sizes. Bond buffers are only filled for
      // owned atoms anyway. A bond to a ghost resolves, by global id, to the owned atom it images.
      std::vector<size_t> owned; // compact index -> flat particle index
      owned.reserve( n_particles );
      std::vector<double> struct_type;
      std::vector<double> orientation_mat3;
      std::vector<uint64_t> global_id;
      std::vector<uint64_t> bond_id;
      std::vector<double> bond_distance, bond_disorientation;
      std::vector<int> bond_count;
      {
        auto cells = grid->cells_accessor();
        auto type_acc = grid->field_accessor( field::mk_generic_real( *struct_field ) );
        auto orient_acc = grid->field_accessor( field::mk_generic_mat3( *orient_field ) );
        for(size_t c=0;c<n_cells;c++)
        {
          if( grid->is_ghost_cell(c) ) { continue; }
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            const size_t i = cell_particle_offset[c] + p;
            owned.push_back( i );
            struct_type.push_back( cells[c][type_acc][p] );
            const Mat3d R = cells[c][orient_acc][p];
            const double m[9] = { R.m11,R.m12,R.m13, R.m21,R.m22,R.m23, R.m31,R.m32,R.m33 };
            orientation_mat3.insert( orientation_mat3.end(), m, m+9 );
            global_id.push_back( cells[c][field::id][p] );
            const size_t base = i*static_cast<size_t>(max_neighbors);
            bond_id.insert( bond_id.end(), grain_bond_id->data()+base, grain_bond_id->data()+base+max_neighbors );
            bond_distance.insert( bond_distance.end(), grain_bond_distance->data()+base, grain_bond_distance->data()+base+max_neighbors );
            bond_disorientation.insert( bond_disorientation.end(), grain_bond_disorientation->data()+base, grain_bond_disorientation->data()+base+max_neighbors );
            bond_count.push_back( (*grain_bond_count)[i] );
          }
        }
      }

      GrainSegmentationResult result;
      grain_segmentation_nps(
        owned.size(), struct_type.data(), orientation_mat3.data(), global_id.data(),
        bond_id.data(), bond_distance.data(), bond_disorientation.data(), bond_count.data(),
        max_neighbors, *auto_threshold, *manual_threshold_log, *min_grain_atom_count, *orphan_adoption, *color_seed,
        result );

      *n_grains = result.n_grains;
      *n_orphans_adopted = result.n_orphans_adopted;
      lout << "compute_grain_clusters: " << result.n_grains << " grains found ("
           << *min_grain_atom_count << "+ atoms each), " << result.n_orphans_adopted
           << " orphan atoms adopted, merge threshold (log) = " << result.merge_threshold_log << std::endl;

      auto id_acc = grid->field_accessor( field::mk_generic_real( *grain_id_field ) );
      auto color_acc = grid->field_accessor( field::mk_generic_vec3( *grain_color_field ) );
      // atom_grain_id is 0-based-"none"/1-based-grain internally (matches union-find root indexing
      // convention); shift to -1="none" for the OUTPUT field.
      // owned atoms get their own result, ghost copies the result of the owned atom with the same id
      // (-1 if that atom is owned by another MPI rank)
      std::unordered_map<uint64_t,int32_t> id_to_grain;
      id_to_grain.reserve( owned.size()*2 );
      std::vector<int32_t> shifted( n_particles, -1 );
      for(size_t k=0;k<owned.size();k++) { shifted[owned[k]] = result.atom_grain_id[k] - 1; id_to_grain[global_id[k]] = shifted[owned[k]]; }
      {
        auto cells = grid->cells_accessor();
        for(size_t c=0;c<n_cells;c++)
        {
          if( ! grid->is_ghost_cell(c) ) { continue; }
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            auto it = id_to_grain.find( cells[c][field::id][p] );
            if( it != id_to_grain.end() ) { shifted[ cell_particle_offset[c] + p ] = it->second; }
          }
        }
      }
      GrainIdFieldFunctor func = { cell_particle_offset, shifted.data(), result.grain_color.empty() ? nullptr : &result.grain_color[0].x };
      compute_cell_particles( *grid, false, func, onika::make_flat_tuple( id_acc, color_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Grain segmentation pass (ii): Node-Pair-Sampling graph clustering of PTM-classified atoms by
misorientation, with auto-threshold selection and optional Dijkstra-by-distance orphan-atom
adoption -- see this file's own header comment for the full algorithm and OVITO fidelity notes.
Writes a per-particle grain_id field (-1 = none) and an RGB grain_color field.

Usage example:

compute_ptm: { rcut: 3.6 ang }
ptm_fields: { struct_field: ptm_type, orient_field: ptm_orientation }
ghost_update_opt: { opt_fields: [ ptm_type, ptm_orientation ] }
compute_grain_bond_misorientation: { rcut: 3.6 ang }
compute_grain_clusters: { min_grain_atom_count: 100, orphan_adoption: true }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_grain_clusters)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_grain_clusters", make_grid_variant_operator< ComputeGrainClusters > );
  }

}
