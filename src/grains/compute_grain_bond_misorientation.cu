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
#include <onika/memory/allocator.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <ptm_constants.h> // PTM_MATCH_* only -- no link dependency, see this plugin's CMakeLists.txt
#include <exaStamp/grains/grain_quat_math.h>

#include <cmath>
#include <cstdint>
#include <limits>

// GPU-compatible pass (i) of grain segmentation, following OVITO's GrainSegmentationEngine1
// (crystalanalysis/modifier/grains): for every chunk-neighbor bond of a PTM-classified atom, records
// the neighbor's global id, real bond length, and -- if both atoms share the same (supported)
// structure type -- the crystal-symmetry-aware disorientation angle. This is exactly OVITO's
// `_neighborBonds` list, reused downstream for BOTH the misorientation-threshold graph-clustering
// step (compute_grain_clusters) and Dijkstra-by-distance orphan-atom adoption -- same reuse OVITO's
// own engine makes (`_neighborBonds`/`noncrystallineBonds` are the same underlying array).
//
// Disorientation formula (`quat_rot`/`rotate_quaternion_into_fundamental_zone`/`quat_disorientation_
// cubic`/`quat_disorientation_hcp_conventional`, and the 24-element cubic / 12-element "conventional"
// hcp point-group quaternion tables) is ported VERBATIM from the vendored PTM library's own
// src/ptm/lib/ptm_quat.cpp (PM Larsen, MIT license, the exact routine OVITO's own grains modifier
// calls) -- duplicated here rather than linked, specifically so this pass can be genuinely
// CudaCompatible=true: the vendored ptm:: functions are plain host C++, not annotated for device
// compilation, and ptm_type/ptm_orientation (compute_ptm.cu) are themselves already a hard host-only
// bottleneck upstream, so linking against them would have silently made "GPU-compatible" a fiction.
//
// Scope matches OVITO's own GrainSegmentation: only FCC/BCC/SC (cubic, 24-element symmetry group)
// and HCP (12-element "conventional" hexagonal group) support disorientation/clustering -- ICO has
// no supported symmetry group there either (same gap OVITO itself has, see this feature's own
// investigation notes), diamond/graphene aren't checked by compute_ptm.cu in the first place.
// Coherent-interface handling (OVITO's optional HCP-in-FCC twin/stacking-fault reclassification,
// `handleCoherentInterfaces`) is NOT implemented -- not requested, and a materially separate feature
// (its own Dijkstra flood-fill + hcp<->fcc quaternion mapping); every HCP-inside-FCC boundary is
// treated as a genuine grain boundary for now.
namespace exaStamp
{
  using namespace exanb;

  static constexpr int MAX_GRAIN_NEIGHBORS = 24;
  static constexpr uint64_t GRAIN_BOND_EMPTY = std::numeric_limits<uint64_t>::max();

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) GrainBondExtStorage
  {
    int central_type = PTM_MATCH_NONE;
    double central_quat[4] = {1.,0.,0.,0.};
    uint64_t nbh_id[MAX_GRAIN_NEIGHBORS];
    double nbh_distance[MAX_GRAIN_NEIGHBORS] = {};
    double nbh_disorientation[MAX_GRAIN_NEIGHBORS] = {};
    int count = 0;

    ONIKA_HOST_DEVICE_FUNC inline void reset() { count = 0; }
  };

  template<class TypeFieldT, class OrientFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) GrainBondMisorientationFunctor
  {
    const double m_threshold_deg = 4.0;
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    TypeFieldT m_type_field = {};
    OrientFieldT m_orient_field = {};
    uint64_t * const __restrict__ m_out_id = nullptr;
    double   * const __restrict__ m_out_distance = nullptr;
    double   * const __restrict__ m_out_disorientation = nullptr;
    int      * const __restrict__ m_out_count = nullptr;

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart ) const
    {
      ctx.ext.reset();
      ctx.ext.central_type = static_cast<int>( cells[cell_a][m_type_field][p_a] );
      if( ctx.ext.central_type != PTM_MATCH_NONE )
      {
        const Mat3d Ra = cells[cell_a][m_orient_field][p_a];
        gb_matrix_to_quat( Ra, ctx.ext.central_quat );
      }
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
       ComputeBufferT& ctx
      , const Vec3d& dr, double d2
      , LocalCellsT cells, size_t cell_b, size_t p_b
      , double /*scale*/) const
    {
      const double dist = sqrt(d2);
      int idx;
      if( ctx.ext.count < MAX_GRAIN_NEIGHBORS )
      {
        idx = ctx.ext.count++;
      }
      else
      {
        // Full: keep the MAX_GRAIN_NEIGHBORS(24) NEAREST bonds seen so far, not just the first ones
        // encountered in traversal order (same "sorted by true distance" principle compute_ptm.cu's
        // own PTMFunctor already applies) -- an unsorted cap can silently drop the closest bonds,
        // which matters both for clustering fidelity and, more importantly, for orphan-adoption
        // reachability (a dropped near bond can be the only short path across a boundary).
        int farthest = 0; double farthest_d = ctx.ext.nbh_distance[0];
        for(int k=1;k<MAX_GRAIN_NEIGHBORS;k++) { if( ctx.ext.nbh_distance[k] > farthest_d ) { farthest_d = ctx.ext.nbh_distance[k]; farthest = k; } }
        if( dist >= farthest_d ) { return; }
        idx = farthest;
      }
      ctx.ext.nbh_id[idx] = cells[cell_b][field::id][p_b];
      ctx.ext.nbh_distance[idx] = dist;

      double disorientation = -1.0; // sentinel: not a valid crystalline-clustering candidate bond
      const int type_b = static_cast<int>( cells[cell_b][m_type_field][p_b] );
      if( ctx.ext.central_type != PTM_MATCH_NONE && type_b == ctx.ext.central_type )
      {
        const Mat3d Rb = cells[cell_b][m_orient_field][p_b];
        double qb[4]; gb_matrix_to_quat( Rb, qb );
        const double deg = gb_disorientation_deg( ctx.ext.central_type, ctx.ext.central_quat, qb );
        if( deg >= 0.0 && deg < m_threshold_deg ) { disorientation = deg; }
      }
      ctx.ext.nbh_disorientation[idx] = disorientation;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ComputeBufferT& ctx, LocalCellsT /*cells*/, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
    {
      const size_t i = m_cell_particle_offset[cell_a] + p_a;
      const size_t base = i * static_cast<size_t>(MAX_GRAIN_NEIGHBORS);
      for(int k=0;k<ctx.ext.count;k++)
      {
        m_out_id[base+k] = ctx.ext.nbh_id[k];
        m_out_distance[base+k] = ctx.ext.nbh_distance[k];
        m_out_disorientation[base+k] = ctx.ext.nbh_disorientation[k];
      }
      for(int k=ctx.ext.count;k<MAX_GRAIN_NEIGHBORS;k++) { m_out_id[base+k] = GRAIN_BOND_EMPTY; }
      m_out_count[i] = ctx.ext.count;
    }
  };
}

namespace exanb
{
  template<class TypeFieldT, class OrientFieldT>
  struct ComputePairTraits< exaStamp::GrainBondMisorientationFunctor<TypeFieldT,OrientFieldT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };
}

namespace exaStamp
{
  template<class GridT>
  class ComputeGrainBondMisorientation : public OperatorNode
  {
    ADD_SLOT( GridT                     , grid             , INPUT , REQUIRED );
    ADD_SLOT( Domain                    , domain           , INPUT , REQUIRED );
    ADD_SLOT( double                    , rcut             , INPUT , REQUIRED , DocString{"Neighbor search cutoff -- same convention as compute_ptm's own rcut, should match/exceed it"} );
    ADD_SLOT( double                    , misorientation_threshold , INPUT , 4.0 , DocString{"Hard ceiling (degrees) on a candidate crystalline-clustering bond -- matches OVITO's own hardcoded 4 degree GrainSegmentation ceiling, applied here (not downstream) so the clustering pass only ever sees already-filtered candidates"} );
    ADD_SLOT( std::string               , struct_field     , INPUT , std::string("ptm_type") , DocString{"Name of the per-particle PTM structure-type field (see ptm_fields). Must have fresh ghost data -- run ghost_update_opt on it first."} );
    ADD_SLOT( std::string               , orient_field     , INPUT , std::string("ptm_orientation") , DocString{"Name of the per-particle PTM lattice-orientation rotation-tensor field (see ptm_fields). Must have fresh ghost data -- run ghost_update_opt on it first."} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors  , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double                    , rcut_max         , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( onika::memory::CudaMMVector<uint64_t> , grain_bond_id            , OUTPUT , DocString{"Flat per-particle neighbor-bond global-id buffer, MAX_GRAIN_NEIGHBORS(24) slots/particle, indexed [cell_particle_offset[cell]+particle]*24+k; empty slots hold UINT64_MAX"} );
    ADD_SLOT( onika::memory::CudaMMVector<double>   , grain_bond_distance      , OUTPUT , DocString{"Same shape as grain_bond_id: real Euclidean bond length, for every recorded bond (regardless of crystalline-candidate status)"} );
    ADD_SLOT( onika::memory::CudaMMVector<double>   , grain_bond_disorientation , OUTPUT , DocString{"Same shape as grain_bond_id: disorientation angle in degrees for a same-structure-type bond under the threshold, -1.0 otherwise (not a valid clustering candidate -- still usable for orphan adoption via grain_bond_distance)"} );
    ADD_SLOT( onika::memory::CudaMMVector<int>      , grain_bond_count         , OUTPUT , DocString{"Per-particle number of valid entries actually written (<=24)"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( grid->number_of_cells() == 0 ) return;

      if( ! grid->has_allocated_field( field::mk_generic_real( *struct_field ) ) )
      {
        fatal_error() << "compute_grain_bond_misorientation: input field '" << *struct_field << "' does not exist (run compute_ptm + ptm_fields first)" << std::endl;
      }
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *orient_field ) ) )
      {
        fatal_error() << "compute_grain_bond_misorientation: input field '" << *orient_field << "' does not exist (run compute_ptm + ptm_fields first)" << std::endl;
      }

      const size_t total_particles = grid->number_of_particles();
      grain_bond_id->resize( total_particles * MAX_GRAIN_NEIGHBORS );
      grain_bond_distance->resize( total_particles * MAX_GRAIN_NEIGHBORS );
      grain_bond_disorientation->resize( total_particles * MAX_GRAIN_NEIGHBORS );
      grain_bond_count->resize( total_particles );

      auto type_acc = grid->field_accessor( field::mk_generic_real( *struct_field ) );
      auto orient_acc = grid->field_accessor( field::mk_generic_mat3( *orient_field ) );

      using ComputeBuffer = ComputePairBuffer2<false,false,GrainBondExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      GrainBondMisorientationFunctor<decltype(type_acc),decltype(orient_acc)> compute_op =
        { *misorientation_threshold, grid->cell_particle_offset_data(), type_acc, orient_acc
        , grain_bond_id->data(), grain_bond_distance->data(), grain_bond_disorientation->data(), grain_bond_count->data() };

      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{}, cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr DefaultPositionFields posfields = {};
      static constexpr std::integral_constant<bool,true> force_use_cells_accessor = {};
      compute_cell_particle_pairs2( *grid, *rcut, false, optional, compute_buf, compute_op, compute_field_set
                                   , posfields, parallel_execution_context(), force_use_cells_accessor );

      long n_candidate_bonds = 0, n_total_bonds = 0;
      for(size_t i=0;i<total_particles;i++)
      {
        n_total_bonds += (*grain_bond_count)[i];
        for(int k=0;k<(*grain_bond_count)[i];k++) { if( (*grain_bond_disorientation)[i*MAX_GRAIN_NEIGHBORS+k] >= 0.0 ) { ++n_candidate_bonds; } }
      }
      lout << "compute_grain_bond_misorientation: " << n_total_bonds << " neighbor bonds recorded, "
           << n_candidate_bonds << " under the " << (*misorientation_threshold) << " deg threshold (candidate crystalline bonds)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Pass (i) of grain segmentation: for every chunk-neighbor bond, records the neighbor's global id, real
distance, and (if both atoms share the same supported PTM structure type -- FCC/BCC/SC or HCP) the
crystal-symmetry-aware disorientation angle in degrees, gated at misorientation_threshold (OVITO's own
hardcoded 4 degree ceiling by default). GPU-compatible. Feeds compute_grain_clusters (Node-Pair-
Sampling clustering) and its own orphan-atom-adoption pass.

Usage example:

compute_ptm: { rcut: 3.6 ang }
ptm_fields: { struct_field: ptm_type, orient_field: ptm_orientation }
ghost_update_opt: { opt_fields: [ ptm_type, ptm_orientation ] }
compute_grain_bond_misorientation: { rcut: 3.6 ang }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_grain_bond_misorientation)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_grain_bond_misorientation", make_grid_variant_operator< ComputeGrainBondMisorientation > );
  }

}
