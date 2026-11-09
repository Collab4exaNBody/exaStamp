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
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exaStamp/analysis_particle/bounded_pair_buffer.h>

#include <ptm_functions.h>
#include <ptm_constants.h>
#include <ptm_initialize_data.h>

#include <omp.h>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <algorithm>

// Polyhedral Template Matching (Larsen, Schmidt, Schiotz, MSMSE 24 (2016) 055007, vendored library in
// lib/ and include/): structure type (FCC, HCP, BCC, ICO, SC), lattice orientation and RMSD per
// particle. Host only, one PTM handle per OpenMP thread. Needs chunk_size 1 neighbor lists.
namespace exaStamp
{
  using namespace exanb;

  // single-shell structures query their neighbors once: return the prepared environment
  static int ptm_get_neighbours_from_prebuilt_env(void* vdata, size_t /*unused*/, size_t /*atom_index*/, int /*num*/, ptm_atomicenv_t* env)
  {
    *env = *reinterpret_cast<const ptm_atomicenv_t*>(vdata);
    return 0;
  }

  // PTM uses at most the PTM_MAX_INPUT_POINTS-1 nearest neighbors
  using PTMComputeBuffer = NearestNeighborsComputeBuffer< PTM_MAX_INPUT_POINTS - 1 >;

  struct PTMFunctor
  {
    // DCUB/DHEX/GRAPHENE need a two-shell neighbor query, not implemented
    static constexpr int32_t ENABLED_FLAGS = PTM_CHECK_SC | PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC;
    const int32_t m_flags = ENABLED_FLAGS;
    const double m_rmsd_cutoff = 0.1; // <= 0 disables the cutoff
    ptm_local_handle_t * const m_thread_handles = nullptr; // one per OpenMP thread, sized by caller
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    double * const __restrict__ m_struct_type_out = nullptr; // 1 per particle
    double * const __restrict__ m_orientation_out = nullptr; // 4 per particle (quaternion w,x,y,z)
    double * const __restrict__ m_rmsd_out = nullptr;        // 1 per particle

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () ( int jnum, ComputeBufferT& buf, CellParticlesT /*cells*/) const
    {
      // the 18 nearest neighbors, sorted by distance (NearestNeighborsComputeBuffer), whatever rcut is
      const int n = std::min( jnum, static_cast<int>(PTM_MAX_INPUT_POINTS) - 1 );

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
          // PTM always returns its best match: reject poor fits, keep rmsd and orientation for diagnostics
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
    ADD_SLOT( MPI_Comm            , mpi             , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT               , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain              , domain          , INPUT , REQUIRED );
    ADD_SLOT( double              , rcut            , INPUT , REQUIRED , DocString{"Neighbor search cutoff, must include 14 neighbors for BCC (12 FCC/HCP/ICO, 6 SC). Only the 18 nearest are used"} );
    ADD_SLOT( double              , rmsd_cutoff     , INPUT , 0.1 , DocString{"Matches with a larger RMSD are reported as no structure (RMSD and orientation still written). <= 0 disables it"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double              , rcut_max        , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_struct_type , OUTPUT , DocString{"Per-particle structure type (0 none, 1 FCC, 2 HCP, 3 BCC, 4 ICO, 5 SC), indexed by cell_particle_offset[cell] + particle"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_orientation , OUTPUT , DocString{"Per-particle orientation quaternion (w,x,y,z), 4 values per particle, same indexing"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_rmsd , OUTPUT , DocString{"Per-particle RMSD (negative if no match), same indexing"} );
    ADD_SLOT( long , ptm_n_matched , OUTPUT , DocString{"Number of matched particles on this rank"} );

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
      // ghost entries are not computed: give them the "no match" values instead of stale ones
      std::fill( ptm_struct_type->begin(), ptm_struct_type->end(), static_cast<double>(PTM_MATCH_NONE) );
      std::fill( ptm_rmsd->begin(), ptm_rmsd->end(), -1.0 );
      for(size_t p=0;p<total_particles;p++)
      {
        double * q = ptm_orientation->data() + 4*p;
        q[0] = 1.0; q[1] = 0.0; q[2] = 0.0; q[3] = 0.0;
      }

      using ComputeBuffer = PTMComputeBuffer;
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
                                         , compute_field_set, cp_opts
                                         , posfields, std::index_sequence<>{} );
      }

      size_t n_matched = 0;
      size_t n_owned = 0;
      double rmsd_sum = 0.0;
      double angle_sum = 0.0; // deviation from identity orientation, degrees
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      for(ssize_t k=gl;k<dims.k-gl;k++)
      for(ssize_t j=gl;j<dims.j-gl;j++)
      for(ssize_t i=gl;i<dims.i-gl;i++)
      {
        const size_t cell = static_cast<size_t>( grid_ijk_to_index(dims,IJK{i,j,k}) );
        const size_t np = cells[cell].size();
        for(size_t q=0;q<np;q++)
        {
          const size_t p = cell_particle_offset[cell] + q;
          ++n_owned;
          if( (*ptm_struct_type)[p] != static_cast<double>(PTM_MATCH_NONE) )
          {
            ++n_matched;
            rmsd_sum += (*ptm_rmsd)[p];
            const double qw = std::min( 1.0, std::abs( (*ptm_orientation)[4*p+0] ) );
            angle_sum += 2.0 * std::acos(qw) * (180.0/M_PI);
          }
        }
      }
      *ptm_n_matched = static_cast<long>(n_matched);

      // log: totals over all ranks
      double sums[4] = { double(n_matched), double(n_owned), rmsd_sum, angle_sum };
      MPI_Allreduce( MPI_IN_PLACE, sums, 4, MPI_DOUBLE, MPI_SUM, *mpi );
      const long g_matched = static_cast<long>( sums[0] );
      lout << "compute_ptm: " << g_matched << " / " << static_cast<long>(sums[1]) << " particles matched"
           << ( g_matched>0 ? ( ", mean rmsd=" + std::to_string(sums[2]/sums[0])
                                 + ", mean |angle from identity|=" + std::to_string(sums[3]/sums[0]) + " deg" )
                             : std::string() ) << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Polyhedral Template Matching: structure type (0 none, 1 FCC, 2 HCP, 3 BCC, 4 ICO, 5 SC), lattice
orientation quaternion (w,x,y,z) and RMSD per particle, from the 18 nearest neighbors within rcut.
Writes the flat buffers ptm_struct_type, ptm_orientation and ptm_rmsd; use ptm_fields to copy them
into grid fields. A match with an RMSD above rmsd_cutoff is reported as no structure. Host only.

Usage example:

compute_ptm:
  rcut: 3.6 ang
  rmsd_cutoff: 0.1

)EOF";
    }
  };

  // compute_ptm buffers -> grid fields (orientation as a rotation matrix)
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

Copies compute_ptm's buffers into per-particle grid fields, orientation as a rotation matrix. Owned
particles only: run ghost_update_opt on the fields if ghosts need them.

Usage example:

compute_ptm: { rcut: 3.6 ang }
ptm_fields: { struct_field: ptm_type, orient_field: ptm_orientation, rmsd_field: ptm_rmsd }

)EOF";
    }
  };

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
