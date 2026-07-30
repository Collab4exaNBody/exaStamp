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
#include <onika/math/basic_types.h>
#include <onika/memory/allocator.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_cell_particle_pairs_chunk.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <ptm_functions.h>
#include <ptm_constants.h>
#include <ptm_initialize_data.h>

#include <omp.h>
#include <vector>
#include <cmath>
#include <algorithm>

// Per-particle structure identification + local lattice orientation via PTM (Polyhedral Template
// Matching, Larsen, Schmidt & Schiotz, Model. Simul. Mater. Sci. Eng. 24 (2016) 055007), vendored
// from OVITO's bundled copy (src/ptm/lib, src/ptm/include -- see PTM_LICENSE). Supersedes the
// earlier hand-rolled correspondence-refinement approach: PTM's own quaternion-based rotation fit
// (Theobald's QCP, ptm_polar.cpp) plus its combinatorial template/graph matching is the
// established, validated method for this, not something to re-derive from scratch. This gives
// structure type, orientation quaternion and RMSD in one call per particle -- the orientation is
// PTM's own lattice-registration result, unrelated to (and not derived from) any
// deformation-gradient/polar-decomposition computation elsewhere in this codebase.
//
// PTM itself is host-only (combinatorial graph/template matching, not GPU-portable), so this
// operator's functor is intentionally NOT CudaCompatible -- runs on CPU/OpenMP only, one
// ptm_local_handle_t scratch object per OpenMP thread (PTM's handles are not thread-safe to share).
// Uses the buffer-based compute-pair pattern (all jnum neighbors available at once via buf.drx/
// dry/drz/d2), same shape already used for host/CudaCompatible=false pair ops in this codebase
// (see exaStamp/potential/eam_potential_template/eam_force_op_singlemat.h).
//
// Output goes to flat onika::memory::CudaMMVector buffers indexed via grid->cell_particle_offset_
// data() (exactly compute_bispectrum.cu's pattern), NOT named grid fields -- keeps this operator
// off compute_cell_particle_pairs2's deep cells_accessor() path entirely.
//
// compute_cell_particle_pairs2 itself is NOT used here, even so: its own CS=1/CS=VARIMPL dispatch
// macro (compute_cell_particle_pairs.h) unconditionally compiles BOTH branches for any chunk-
// neighbor caller, and the CS=VARIMPL (runtime chunk size) branch cannot compile for a
// CudaCompatible=false functor on this exaNBody version -- none of compute_cell_particle_pairs_
// cell's three overloads (chunk.h/chunk_cs1.h/chunk_scb.h) accept a genuinely runtime chunk size,
// only the compile-time onika::UIntConst<CS>. This is apparently untested for CudaCompatible=false
// on this GPU build: eam_force_op_singlemat.h's own CudaCompatible=false branch is dead code here
// (USTAMP_POTENTIAL_CUDA_COMPATIBLE resolves true), and compute_bispectrum.cu's CudaCompatible=true
// functor never exercises this code path either. Confirmed CudaCompatible=false is what triggers
// it (not the cells-accessor question above): both attempts hit the identical compile error.
//
// Fix: call compute_cell_particle_pairs_cell directly with a hardcoded onika::UIntConst<1>{},
// mirroring exactly what ComputeParticlePairFunctor::operator() does internally on the host path,
// but skipping compute_cell_particle_pairs2's broken outer dispatch. This requires chunk_neighbors
// to actually be built with chunk_size==1 (checked at runtime below, fatal_error otherwise) --
// set `random_access: true` in the .msp's chunk_neighbors config to guarantee this (see
// chunk_neighbors_config.h: random_access forces chunk_size to 1).
//
// Checks FCC, HCP, BCC, ICO and SC (PTM's single-shot-callback structures, see PTMFunctor::
// m_flags below). DCUB/DHEX/GRAPHENE are deliberately excluded -- they need a two-shell neighbour
// query our callback doesn't implement (see the comment on m_flags).
namespace exaStamp
{
  using namespace exanb;

  // PTM's neighbour_ordering for BCC/FCC/HCP/ICO/SC (num_outer==0) calls this callback exactly
  // once, with the atom's own request -- since we've already gathered and sorted this particle's
  // neighbors before calling ptm_index(), just hand back the pre-built environment as-is.
  static int ptm_get_neighbours_from_prebuilt_env(void* vdata, size_t /*unused*/, size_t /*atom_index*/, int /*num*/, ptm_atomicenv_t* env)
  {
    *env = *reinterpret_cast<const ptm_atomicenv_t*>(vdata);
    return 0;
  }

  struct PTMFunctor
  {
    // FCC/HCP/BCC/ICO/SC all go through PTM's single-shot neighbour_ordering (num_outer==0,
    // ptm_index.cpp) -- the one our ptm_get_neighbours_from_prebuilt_env callback above actually
    // supports. DCUB/DHEX/GRAPHENE use a genuinely different *two-shell* ordering (num_inner/
    // num_outer both nonzero) that calls the neighbour callback multiple times for distinct
    // inner/outer shell requests; our callback ignores that and always hands back the same single
    // prebuilt env, so enabling those flags would silently produce wrong results, not just
    // unsupported ones. Left out until the callback is taught the two-shell protocol.
    static constexpr int32_t ENABLED_FLAGS = PTM_CHECK_SC | PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC;
    const int32_t m_flags = ENABLED_FLAGS;
    const double m_rmsd_cutoff = 0.1; // OVITO's own PTM modifier default; <=0 disables the cutoff
    ptm_local_handle_t * const m_thread_handles = nullptr; // one per OpenMP thread, sized by caller
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    double * const __restrict__ m_struct_type_out = nullptr; // 1 per particle
    double * const __restrict__ m_orientation_out = nullptr; // 4 per particle (quaternion w,x,y,z)
    double * const __restrict__ m_rmsd_out = nullptr;        // 1 per particle

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () ( int jnum, ComputeBufferT& buf, CellParticlesT /*cells*/) const
    {
      // PTM wants its own fixed candidate-neighbor pool -- not "however many happen to lie
      // within rcut" -- rcut only needs to be generous enough that the candidate pool below is
      // full of *genuinely* nearest neighbors; anything past that must be discarded, or PTM's
      // own graph matching sees a different point set and gives different (rcut-dependent)
      // results. Same convention as OVITO's own PTM code (PTMAlgorithm.cpp): always gather up to
      // PTM_MAX_INPUT_POINTS-1 nearest neighbors regardless of which structures are enabled --
      // ptm_index() itself already gates each individual structure check on whether enough points
      // were supplied (see the `num_points >= PTM_NUM_POINTS_*` checks in ptm_index.cpp), and its
      // combinatorial matching is designed to pick the right subset/correspondence out of a
      // candidate pool larger than any one template needs. Partial selection sort finds the true
      // n globally-nearest neighbors (inner loop scans the full jnum, not just the first n) --
      // buf.copy(src,dst) is a one-way overwrite (not a swap), so exchange the two slots by hand.
      const int n = std::min( jnum, static_cast<int>(PTM_MAX_INPUT_POINTS) - 1 );
      for(int i=0;i<n;i++)
      {
        int m = i;
        for(int j=i+1;j<jnum;j++) { if( buf.d2[j] < buf.d2[m] ) { m = j; } }
        if( m != i )
        {
          std::swap( buf.drx[i], buf.drx[m] );
          std::swap( buf.dry[i], buf.dry[m] );
          std::swap( buf.drz[i], buf.drz[m] );
          std::swap( buf.d2[i], buf.d2[m] );
        }
      }

      double struct_type = static_cast<double>( PTM_MATCH_NONE );
      double quat[4] = { 1.0, 0.0, 0.0, 0.0 }; // identity
      double rmsd = -1.0; // negative = no match attempted or found

      if( n >= PTM_NUM_NBRS_SC ) // SC (6) is the smallest requirement among the enabled checks
      {
        ptm_atomicenv_t env;
        env.num = n + 1;
        env.points[0][0] = 0.0; env.points[0][1] = 0.0; env.points[0][2] = 0.0;
        env.atom_indices[0] = 0;
        env.numbers[0] = 1;
        for(int i=0;i<n;i++)
        {
          env.points[i+1][0] = buf.drx[i];
          env.points[i+1][1] = buf.dry[i];
          env.points[i+1][2] = buf.drz[i];
          env.atom_indices[i+1] = static_cast<size_t>(i+1);
          env.numbers[i+1] = 1;
        }

        const int tid = omp_get_thread_num();
        ptm_result_t result;
        ptm_atomicenv_t output_env;
        ptm_index( m_thread_handles[tid], 0, ptm_get_neighbours_from_prebuilt_env, &env,
                   m_flags, false, &result, &output_env );

        if( result.structure_type != PTM_MATCH_NONE )
        {
          // PTM always returns its best combinatorial match, however poor -- it never refuses on
          // its own. Without a cutoff, a genuinely different local structure (e.g. FCC atoms
          // against the BCC-only template checked here) still gets *some* correspondence, just a
          // bad one (rmsd far from 0). Reject those explicitly, same convention OVITO's own PTM
          // modifier uses (rmsdCutoff, default 0.1) -- rmsd/orientation are kept either way, only
          // the discrete structure_type classification is gated, so a rejected fit is still
          // visible for diagnostics.
          rmsd = result.rmsd;
          quat[0] = result.orientation[0];
          quat[1] = result.orientation[1];
          quat[2] = result.orientation[2];
          quat[3] = result.orientation[3];
          if( m_rmsd_cutoff <= 0.0 || rmsd <= m_rmsd_cutoff )
          {
            struct_type = static_cast<double>( result.structure_type );
          }
        }
      }

      const size_t i = m_cell_particle_offset[buf.cell] + buf.part;
      m_struct_type_out[i] = struct_type;
      m_orientation_out[4*i+0] = quat[0];
      m_orientation_out[4*i+1] = quat[1];
      m_orientation_out[4*i+2] = quat[2];
      m_orientation_out[4*i+3] = quat[3];
      m_rmsd_out[i] = rmsd;
    }
  };

  template<class GridT>
  class ComputePTM : public OperatorNode
  {
    ADD_SLOT( GridT               , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain              , domain          , INPUT , REQUIRED );
    ADD_SLOT( double              , rcut            , INPUT , REQUIRED , DocString{"Neighbor search cutoff -- must be generous enough to include at least BCC's 14 first+second-shell neighbors (the largest requirement among the checked structures: FCC/HCP/ICO 12, BCC 14, SC 6)"} );
    ADD_SLOT( double              , rmsd_cutoff     , INPUT , 0.1 , DocString{"Reject a PTM match whose RMSD exceeds this (structure_type reported as PTM_MATCH_NONE, rmsd/orientation still written for diagnostics). Same default as OVITO's own PTM modifier. <=0 disables the cutoff."} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double              , rcut_max        , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_struct_type , OUTPUT , DocString{"Flat per-particle structure type buffer (0=none, 1=FCC, 2=HCP, 3=BCC, 4=ICO, 5=SC -- see PTM_MATCH_* in ptm_constants.h): ptm_struct_type[ cell_particle_offset[cell] + particle ], see grid->cell_particle_offset_data()"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_orientation , OUTPUT , DocString{"Flat per-particle lattice orientation quaternion buffer (w,x,y,z, PTM's own fit result): ptm_orientation[ 4*(cell_particle_offset[cell]+particle) + component ]"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_rmsd , OUTPUT , DocString{"Flat per-particle fit-quality RMSD buffer (negative if no match attempted/found), same indexing as ptm_struct_type"} );
    ADD_SLOT( long , ptm_n_matched , OUTPUT , DocString{"Number of particles PTM found a match for (out of grid->number_of_particles())"} );

    std::vector<ptm_local_handle_t> m_thread_handles;

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( grid->number_of_cells() == 0 ) return;

      if( chunk_neighbors->m_chunk_size != 1 )
      {
        fatal_error() << "compute_ptm: requires chunk_neighbors built with chunk_size=1 (set random_access: true"
                       << " on the chunk_neighbors config), got chunk_size=" << chunk_neighbors->m_chunk_size << std::endl;
      }

      ptm_initialize_global();
      const size_t nt = static_cast<size_t>( omp_get_max_threads() );
      if( nt > m_thread_handles.size() )
      {
        const size_t old_nt = m_thread_handles.size();
        m_thread_handles.resize(nt);
        for(size_t i=old_nt;i<nt;i++) { m_thread_handles[i] = ptm_initialize_local(); }
      }

      const size_t total_particles = grid->number_of_particles();
      ptm_struct_type->resize( total_particles );
      ptm_orientation->resize( total_particles * 4 );
      ptm_rmsd->resize( total_particles );

      using ComputeBuffer = ComputePairBuffer2<false,false>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      PTMFunctor compute_op = { PTMFunctor::ENABLED_FLAGS, *rmsd_cutoff, m_thread_handles.data(), grid->cell_particle_offset_data()
                               , ptm_struct_type->data(), ptm_orientation->data(), ptm_rmsd->data() };

      ComputePairNullWeightIterator cp_weight{};
      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr DefaultPositionFields posfields = {};
      static constexpr ComputeParticlePairOpts<false,true,false> cp_opts = {}; // Symmetric=false, PreferComputeBuffer=true

      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();
      const auto cells = grid->cells();
      const double rcut2 = (*rcut) * (*rcut);

#     pragma omp parallel for collapse(3) schedule(dynamic)
      for(ssize_t k=gl;k<dims.k-gl;k++)
      for(ssize_t j=gl;j<dims.j-gl;j++)
      for(ssize_t i=gl;i<dims.i-gl;i++)
      {
        const IJK cell_a_loc{i,j,k};
        const size_t cell_a = static_cast<size_t>( grid_ijk_to_index(dims,cell_a_loc) );
        compute_cell_particle_pairs_cell( cells, dims, cell_a_loc, cell_a, rcut2
                                         , compute_buf, optional, compute_op
                                         , compute_field_set, onika::UIntConst<1>{}, cp_opts
                                         , posfields, std::index_sequence<>{} );
      }

      size_t n_matched = 0;
      double rmsd_sum = 0.0;
      double angle_sum = 0.0; // deviation from identity orientation, degrees
      for(size_t p=0;p<total_particles;p++)
      {
        if( (*ptm_struct_type)[p] != static_cast<double>(PTM_MATCH_NONE) )
        {
          ++n_matched;
          rmsd_sum += (*ptm_rmsd)[p];
          const double qw = std::min( 1.0, std::abs( (*ptm_orientation)[4*p+0] ) );
          angle_sum += 2.0 * std::acos(qw) * (180.0/M_PI);
        }
      }
      *ptm_n_matched = static_cast<long>(n_matched);
      lout << "compute_ptm: " << n_matched << " / " << total_particles << " particles matched"
           << ( n_matched>0 ? ( ", mean rmsd=" + std::to_string(rmsd_sum/n_matched)
                                 + ", mean |angle from identity|=" + std::to_string(angle_sum/n_matched) + " deg" )
                             : std::string() ) << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Per-particle structure identification and local lattice orientation via PTM (Polyhedral Template
Matching), vendored from OVITO. Checks FCC, HCP, BCC, ICO and SC (not DCUB/DHEX/GRAPHENE, which
need a neighbour query this operator doesn't implement -- see compute_ptm.cu). Gathers each
particle's nearest neighbors (sorted by true distance, up to PTM_MAX_INPUT_POINTS-1=18 within
rcut), feeds them to PTM's ptm_index(), and writes three flat per-particle buffers, indexed via
grid->cell_particle_offset_data() (same convention as compute_bispectrum's output): ptm_struct_type
(0=none, 1=FCC, 2=HCP, 3=BCC, 4=ICO, 5=SC), ptm_orientation
(PTM's own fitted lattice orientation quaternion, w,x,y,z, 4 doubles/particle), and ptm_rmsd
(fit quality, negative if no match attempted/found). A match whose RMSD exceeds rmsd_cutoff is
reported as PTM_MATCH_NONE (rmsd/orientation are still written, only the classification is
gated) -- PTM itself always returns its best combinatorial match, however poor, so this cutoff is
what actually rejects a genuinely different local structure. Host/OpenMP only -- PTM is not
GPU-portable.

Usage example:

compute_ptm:
  rcut: 3.6 ang   # generous enough for BCC's 8 first + 6 second shell neighbors
  rmsd_cutoff: 0.1

)EOF";
    }
  };

  // Materializes compute_ptm's flat CudaMMVector buffers into named per-particle grid fields
  // (field::mk_generic_real/mk_generic_mat3), so any consumer that expects a real grid field
  // (write_xyz, write_grid_vtk, write_delaunay_vtk's color_field, ...) can read PTM's results
  // directly. Purely pointwise, no neighbor search needed -- unlike compute_ptm.cu itself, this
  // has no chunk-neighbor/CS complication to work around, so it uses exanb's plain per-particle
  // compute_cell_particles (compute_cell_particles.h) instead of the pair-compute machinery.
  struct PTMFieldsFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const double * const __restrict__ m_struct_type_in = nullptr;
    const double * const __restrict__ m_orientation_in = nullptr; // quaternion, 4/particle
    const double * const __restrict__ m_rmsd_in = nullptr;

    inline void operator () ( size_t cell, unsigned int part, double& struct_type_out, Mat3d& orient_out, double& rmsd_out ) const
    {
      const size_t i = m_cell_particle_offset[cell] + part;
      struct_type_out = m_struct_type_in[i];
      rmsd_out = m_rmsd_in[i];
      const double q0=m_orientation_in[4*i+0], q1=m_orientation_in[4*i+1], q2=m_orientation_in[4*i+2], q3=m_orientation_in[4*i+3];
      orient_out = Mat3d{
        1.-2.*(q2*q2+q3*q3),   2.*(q1*q2-q0*q3),     2.*(q1*q3+q0*q2),
        2.*(q1*q2+q0*q3),     1.-2.*(q1*q1+q3*q3),   2.*(q2*q3-q0*q1),
        2.*(q1*q3-q0*q2),     2.*(q2*q3+q0*q1),     1.-2.*(q1*q1+q2*q2)
      };
    }
  };

  template<class GridT>
  class PTMFields : public OperatorNode
  {
    ADD_SLOT( GridT , grid , INPUT_OUTPUT );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_struct_type , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_orientation , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_rmsd        , INPUT , REQUIRED );
    ADD_SLOT( std::string , struct_field , INPUT , std::string("ptm_type")        , DocString{"Name of the resulting per-particle structure-type grid field"} );
    ADD_SLOT( std::string , orient_field , INPUT , std::string("ptm_orientation") , DocString{"Name of the resulting per-particle lattice-orientation rotation-tensor grid field"} );
    ADD_SLOT( std::string , rmsd_field   , INPUT , std::string("ptm_rmsd")        , DocString{"Name of the resulting per-particle RMSD grid field"} );

  public:
    inline void execute () override final
    {
      auto struct_acc = grid->field_accessor( field::mk_generic_real( *struct_field ) );
      auto orient_acc = grid->field_accessor( field::mk_generic_mat3( *orient_field ) );
      auto rmsd_acc = grid->field_accessor( field::mk_generic_real( *rmsd_field ) );

      PTMFieldsFunctor func = { grid->cell_particle_offset_data(), ptm_struct_type->data(), ptm_orientation->data(), ptm_rmsd->data() };
      compute_cell_particles( *grid, false, func, onika::make_flat_tuple( struct_acc, orient_acc, rmsd_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Copies compute_ptm's flat per-particle output buffers (ptm_struct_type/ptm_orientation/ptm_rmsd)
into named per-particle grid fields, so any consumer expecting a real grid field (write_xyz,
write_grid_vtk, write_delaunay_vtk's color_field, ...) can read PTM's results directly. Purely
pointwise (no neighbor search) -- run any time after compute_ptm. Orientation is converted from
PTM's quaternion (as stored in ptm_orientation) to a rotation tensor field. Only touches owned
particles -- ghost copies of these fields are not synchronized across MPI ranks.

Usage example:

compute_ptm: { rcut: 3.6 ang }
ptm_fields: { struct_field: ptm_type, orient_field: ptm_orientation, rmsd_field: ptm_rmsd }

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_ptm)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_ptm", make_grid_variant_operator< ComputePTM > );
    OperatorNodeFactory::instance()->register_factory( "ptm_fields", make_grid_variant_operator< PTMFields > );
  }

}

namespace exanb
{
  template<>
  struct ComputePairTraits< exaStamp::PTMFunctor >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible    = false;
    static inline constexpr bool CudaCompatible          = false; // PTM is host-only
  };
}
