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
#include <onika/string_utils.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

namespace exaStamp
{
  using namespace exanb;

  // maximum number of nearest neighbors (nnn) : they live in the per-particle compute context
  // (a per-thread local), so this bounds the context size (~1 KB)
  static inline constexpr unsigned int LOCAL_CSP_MAX_NNN = 32;

  // per-particle context : the (at most nnn) nearest neighbors seen so far, sorted by increasing distance
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) LocalCspExtStorage
  {
    double m_d2[LOCAL_CSP_MAX_NNN];
    double m_rx[LOCAL_CSP_MAX_NNN];
    double m_ry[LOCAL_CSP_MAX_NNN];
    double m_rz[LOCAL_CSP_MAX_NNN];
    unsigned int m_count = 0;
  };

  // Per-atom centro-symmetry parameter (Kelchner et al., Phys. Rev. B 58, 1998), same algorithm as
  // LAMMPS compute centro/atom and the legacy compute_local_csp_old operator :
  //   CSP_i = sum of the nnn/2 smallest |r_ij + r_ik|^2 over all pairs (j,k) of the nnn nearest neighbors of i
  // The nnn nearest neighbors are selected on the fly by insertion (no per-particle neighbor buffer),
  // and the nnn/2 smallest pair values are also kept by insertion, so no heap allocation nor sort is needed.
  template<class CspFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) LocalCspFunctor
  {
    const unsigned int m_nnn = 12;
    const unsigned int m_nhalf = 6;
    CspFieldT m_csp_field = {};

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, LocalCellsT, size_t, size_t, exanb::ComputePairParticleContextStart) const
    {
      ctx.ext.m_count = 0;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
        ComputeBufferT& ctx
      , const Vec3d& dr , double d2
      , LocalCellsT , size_t , size_t
      , double ) const
    {
      auto & ext = ctx.ext;
      unsigned int pos = 0;
      if( ext.m_count < m_nnn ) { pos = ext.m_count ++ ; }
      else if( d2 < ext.m_d2[m_nnn-1] ) { pos = m_nnn - 1; } // replaces the farthest kept neighbor
      else return;                                            // ties keep the first neighbor seen
      // insertion sort step : shift farther neighbors up (strict comparison keeps the order of equal distances)
      while( pos > 0 && ext.m_d2[pos-1] > d2 )
      {
        ext.m_d2[pos] = ext.m_d2[pos-1];
        ext.m_rx[pos] = ext.m_rx[pos-1];
        ext.m_ry[pos] = ext.m_ry[pos-1];
        ext.m_rz[pos] = ext.m_rz[pos-1];
        -- pos;
      }
      ext.m_d2[pos] = d2;
      ext.m_rx[pos] = dr.x;
      ext.m_ry[pos] = dr.y;
      ext.m_rz[pos] = dr.z;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop) const
    {
      const auto & ext = ctx.ext;
      double csp = 0.0;
      if( ext.m_count >= m_nnn ) // fewer than nnn neighbors within rcut : csp = 0 (as LAMMPS and the legacy operator)
      {
        // nhalf smallest |r_j + r_k|^2 , sorted by increasing value
        double smallest[LOCAL_CSP_MAX_NNN/2];
        unsigned int n_smallest = 0;
        for(unsigned int j=0;j<m_nnn;j++)
        {
          for(unsigned int k=j+1;k<m_nnn;k++)
          {
            const double dx = ext.m_rx[j] + ext.m_rx[k];
            const double dy = ext.m_ry[j] + ext.m_ry[k];
            const double dz = ext.m_rz[j] + ext.m_rz[k];
            const double v = dx*dx + dy*dy + dz*dz;
            unsigned int pos = 0;
            if( n_smallest < m_nhalf ) { pos = n_smallest ++ ; }
            else if( v < smallest[m_nhalf-1] ) { pos = m_nhalf - 1; }
            else continue;
            while( pos > 0 && smallest[pos-1] > v ) { smallest[pos] = smallest[pos-1]; -- pos; }
            smallest[pos] = v;
          }
        }
        for(unsigned int i=0;i<m_nhalf;i++) csp += smallest[i]; // increasing order, as the legacy operator
      }
      cells[cell_a][m_csp_field][p_a] = csp;
    }
  };

  template<class GridT>
  class ComputeLocalCsp : public OperatorNode
  {
    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT , DocString{"Local sub-domain particles grid"} );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED , DocString{"Simulation domain"} );
    ADD_SLOT( double                    , rcut            , INPUT        , REQUIRED , DocString{"Neighbor search radius, must enclose the nnn nearest neighbors"} );
    ADD_SLOT( size_t                    , nnn             , INPUT        , 12 , DocString{"Number of nearest neighbors (even, in [2;32]) : 12 for FCC, 8 for BCC"} );
    ADD_SLOT( std::string               , name            , INPUT        , std::string("csp") , DocString{"Name of the resulting per-particle scalar field"} );
    ADD_SLOT( bool                      , ghost           , INPUT        , false , DocString{"Also compute for ghost particles"} );
    ADD_SLOT( bool                      , verbose         , INPUT        , false , DocString{"Print parameters"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT        , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double                    , rcut_max        , INPUT        , 0.0 , DocString{"Neighbor list cutoff (global rcut_max), must be >= rcut"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( chunk_neighbors->number_of_cells() != grid->number_of_cells() )
      {
        fatal_error() << "compute_local_csp: neighbor list does not match the grid (missing chunk_neighbors ?)" << std::endl;
      }
      if( ! ( *rcut > 0.0 ) ) { fatal_error() << "compute_local_csp: rcut must be > 0 (got " << *rcut << ")" << std::endl; }
      // the neighbor list is sized from the global rcut_max at initialization (nbh_dist), before analysis operators run :
      // neighbors between rcut_max and rcut would be silently missing
      if( *rcut > *rcut_max )
      {
        fatal_error() << "compute_local_csp: rcut=" << *rcut << " exceeds the neighbor list cutoff rcut_max=" << *rcut_max
                      << ", neighbors would be missing. Set 'rcut_max' >= rcut in the global section of the input file" << std::endl;
      }
      if( *nnn < 2 || *nnn > LOCAL_CSP_MAX_NNN || ( *nnn % 2 ) != 0 )
      {
        fatal_error() << "compute_local_csp: nnn=" << *nnn << " must be even and in [2;" << LOCAL_CSP_MAX_NNN << "]" << std::endl;
      }

      if( *verbose )
      {
        lout << onika::format_string("compute_local_csp: rcut=%.5f nnn=%d\n" , *rcut , int(*nnn) );
      }

      auto csp_acc = grid->field_accessor( field::mk_generic_real( *name ) );

      using ComputeBuffer = ComputePairBuffer2<false,false,LocalCspExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      const unsigned int n = static_cast<unsigned int>( *nnn );
      LocalCspFunctor<decltype(csp_acc)> compute_op = { n , n / 2 , csp_acc };

      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr std::integral_constant<bool,true> force_use_cells_accessor = {};
      static constexpr DefaultPositionFields posfields = {};
      compute_cell_particle_pairs2( *grid, *rcut, *ghost, optional, compute_buf, compute_op, compute_field_set
                                   , posfields, parallel_execution_context(), force_use_cells_accessor );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the per-atom centro-symmetry parameter (Kelchner et al., Phys. Rev. B 58, 1998,
same algorithm as LAMMPS compute centro/atom) as a per-particle scalar field,
GPU-compatibility to be tested:

  CSP_i = sum of the nnn/2 smallest |r_ij + r_ik|^2 over all pairs (j,k) of the nnn
          nearest neighbors of atom i

CSP is 0 for a perfect centro-symmetric lattice, and grows near defects and surfaces.
Atoms with fewer than nnn neighbors within rcut get CSP = 0.

Requirements / limits:
  - nnn even, in [2;32] : 12 for FCC, 8 for BCC
  - rcut must enclose the nnn nearest neighbors (e.g. between the 1st and 2nd shell)
  - rcut must not exceed the neighbor list cutoff: set 'rcut_max: <value>' (>= rcut) in the
    global section when rcut is larger than the potential cutoff

Usage example:

compute_local_csp:
  rcut: 3.0 ang
  nnn: 12
  name: csp

)EOF";
    }
  };

  template<class GridT> using ComputeLocalCspTmpl = ComputeLocalCsp<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_local_csp)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_local_csp", make_grid_variant_operator< ComputeLocalCspTmpl > );
  }

}

namespace exanb
{
  // specialize functor traits to allow Cuda execution space
  template<class CspFieldT>
  struct ComputePairTraits< exaStamp::LocalCspFunctor<CspFieldT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };
}
