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

#include <exaStamp/grains/grain_segmentation_algo.h>

#include <ptm_constants.h>

// Grain segmentation, pass (ii): host-side, sequential (same as OVITO's own engine -- neither the
// Node-Pair-Sampling dendrogram construction nor Dijkstra orphan adoption is GPU-parallelizable,
// each merge/settle step depends on prior state). Consumes compute_grain_bond_misorientation's own
// per-atom bond arrays (pass i, GPU-compatible) plus PTM's own struct_type/orientation fields, and
// writes a per-atom grain_id field (-1 default, matching this project's own dxa_mark_core_atoms
// convention for "not part of anything") plus a per-grain RGB color field looked up from the final
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
// NOT implemented (see this feature's own investigation notes for the full comparison against
// OVITO): coherent-interface (twin/stacking-fault) reclassification; cross-rank MPI grain stitching
// (each rank clusters its own local owned+ghost view independently -- a grain spanning a rank
// boundary is NOT currently reconciled into one id, same known scope gap as compute_dxa_lattice_
// clusters' own per-rank BFS).
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

      std::vector<double> struct_type( n_particles );
      std::vector<double> orientation_mat3( n_particles*9 );
      std::vector<uint64_t> global_id( n_particles );
      {
        auto cells = grid->cells_accessor();
        auto type_acc = grid->field_accessor( field::mk_generic_real( *struct_field ) );
        auto orient_acc = grid->field_accessor( field::mk_generic_mat3( *orient_field ) );
        for(size_t c=0;c<n_cells;c++)
        {
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            const size_t i = cell_particle_offset[c] + p;
            struct_type[i] = cells[c][type_acc][p];
            const Mat3d R = cells[c][orient_acc][p];
            orientation_mat3[9*i+0]=R.m11; orientation_mat3[9*i+1]=R.m12; orientation_mat3[9*i+2]=R.m13;
            orientation_mat3[9*i+3]=R.m21; orientation_mat3[9*i+4]=R.m22; orientation_mat3[9*i+5]=R.m23;
            orientation_mat3[9*i+6]=R.m31; orientation_mat3[9*i+7]=R.m32; orientation_mat3[9*i+8]=R.m33;
            global_id[i] = cells[c][field::id][p];
          }
        }
      }

      const int max_neighbors = static_cast<int>( grain_bond_count->size()>0 ? (grain_bond_id->size()/grain_bond_count->size()) : 0 );

      GrainSegmentationResult result;
      grain_segmentation_nps(
        n_particles, struct_type.data(), orientation_mat3.data(), global_id.data(),
        grain_bond_id->data(), grain_bond_distance->data(), grain_bond_disorientation->data(), grain_bond_count->data(),
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
      // convention); shift to -1="none" for the OUTPUT field, matching dxa_mark_core_atoms' own
      // "-1 = not part of anything" convention used elsewhere in this codebase.
      std::vector<int32_t> shifted( n_particles );
      for(size_t i=0;i<n_particles;i++) { shifted[i] = result.atom_grain_id[i] - 1; }
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
