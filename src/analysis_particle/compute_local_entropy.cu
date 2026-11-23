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

#include <exaStamp/thermo_state/thermodynamic_state.h>

#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  // maximum number of radial bins of the smoothed per-atom g(r) : bins live in the per-particle
  // compute context (a per-thread local), so this bounds the context size (~1 KB)
  static inline constexpr size_t LOCAL_ENTROPY_MAX_BINS = 128;

  // per-particle accumulators : unnormalized gaussian-smeared radial histogram and neighbor count
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) LocalEntropyExtStorage
  {
    double m_g[LOCAL_ENTROPY_MAX_BINS];
    unsigned int m_count = 0;
  };

  // Per-atom two-body excess entropy (Piaget & Pablo, J. Chem. Phys. 147, 2017), same
  // discretization as the legacy compute_local_entropy_old operator. Bin radii and the
  // normalization prefactor are evaluated on the fly instead of read from host arrays,
  // and the prefactor/density factor is applied once per bin at the end instead of per neighbor.
  template<class EntropyFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) LocalEntropyFunctor
  {
    static inline constexpr double gtol = 1.0e-10; // below this g(r), the integrand g.ln(g)-g+1 is taken as 1
    static inline constexpr double eps = 1.0e-8;   // guards the bin index of a distance falling exactly on a bin edge

    const size_t m_nbins = 0;
    const size_t m_deltabin = 0;        // half width of the gaussian window, in bins
    const double m_drm = 0.0;           // bin spacing
    const double m_inv_drm = 0.0;
    const double m_inv_sigmasq2 = 0.0;  // 1 / (2 sigma^2)
    const double m_prefac_const = 0.0;  // 1 / (4 pi rho_global sqrt(2 pi) sigma)
    const double m_rho_global = 0.0;
    const double m_local_volume = 0.0;  // 4/3 pi rcut^3
    const bool m_local = false;
    EntropyFieldT m_entropy_field = {};

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, LocalCellsT, size_t, size_t, exanb::ComputePairParticleContextStart) const
    {
      for(size_t k=0;k<m_nbins;k++) ctx.ext.m_g[k] = 0.0;
      ctx.ext.m_count = 0;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
        ComputeBufferT& ctx
      , const Vec3d& , double d2
      , LocalCellsT , size_t , size_t
      , double ) const
    {
      const double rij = sqrt( d2 );
      const size_t bin = static_cast<size_t>( rij * m_inv_drm + eps );
      const size_t minbin = ( bin > m_deltabin ) ? ( bin - m_deltabin ) : 0;
      const size_t maxbin = ( m_nbins < bin + m_deltabin ) ? m_nbins : ( bin + m_deltabin );
      for(size_t k=minbin;k<maxbin;k++)
      {
        const double dr = k * m_drm - rij;
        ctx.ext.m_g[k] += exp( - dr * dr * m_inv_sigmasq2 );
      }
      ++ ctx.ext.m_count;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop) const
    {
      const double rho = m_local ? ( static_cast<double>(ctx.ext.m_count) / m_local_volume ) : m_rho_global;
      double entropy = 0.0;
      if( rho > 0.0 )
      {
        const double rhofac = m_rho_global / rho;
        double integral = 0.0;
        for(size_t k=0;k<m_nbins;k++)
        {
          const double rksq = k * k * m_drm * m_drm;
          const double r1sq = ( k == 0 ) ? ( m_drm * m_drm ) : rksq;   // bin 0 reuses the bin 1 prefactor (avoids r=0 division)
          const double g = ctx.ext.m_g[k] * ( m_prefac_const / r1sq ) * rhofac;
          const double f = ( g >= gtol ) ? ( g * log(g) - g + 1.0 ) : 1.0;
          const double w = ( k == 0 || k == m_nbins - 1 ) ? 0.5 : 1.0;   // trapezoid rule
          integral += f * rksq * w;
        }
        entropy = -2.0 * M_PI * rho * integral * m_drm;
      }
      cells[cell_a][m_entropy_field][p_a] = entropy;
    }
  };

  template<class GridT>
  class ComputeLocalEntropy : public OperatorNode
  {
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT , DocString{"Local sub-domain particles grid"} );
    ADD_SLOT( Domain                    , domain              , INPUT        , REQUIRED , DocString{"Simulation domain"} );
    ADD_SLOT( ThermodynamicState        , thermodynamic_state , INPUT        , REQUIRED , DocString{"Thermodynamic state, provides the global density (particle count / volume). Must be up to date"} );
    ADD_SLOT( double                    , rcut                , INPUT        , REQUIRED , DocString{"Cutoff radius of the local pair distribution function"} );
    ADD_SLOT( double                    , sigma               , INPUT        , REQUIRED , DocString{"Width of the gaussian smearing of neighbor distances"} );
    ADD_SLOT( size_t                    , nbins               , INPUT        , 0 , DocString{"Number of radial bins, 0 means rcut/sigma+1 (must end in [2;128])"} );
    ADD_SLOT( bool                      , local               , INPUT        , false , DocString{"Use the local density (neighbor count / sphere volume) instead of the global one"} );
    ADD_SLOT( std::string               , name                , INPUT        , std::string("entropy") , DocString{"Name of the resulting per-particle scalar field"} );
    ADD_SLOT( bool                      , ghost               , INPUT        , false , DocString{"Also compute for ghost particles"} );
    ADD_SLOT( bool                      , verbose             , INPUT        , false , DocString{"Print parameters"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors     , INPUT        , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double                    , rcut_max            , INPUT        , 0.0 , DocString{"Neighbor list cutoff (global rcut_max), must be >= rcut"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( chunk_neighbors->number_of_cells() != grid->number_of_cells() )
      {
        fatal_error() << "compute_local_entropy: neighbor list does not match the grid (missing chunk_neighbors ?)" << std::endl;
      }
      if( ! ( *rcut > 0.0 ) ) { fatal_error() << "compute_local_entropy: rcut must be > 0 (got " << *rcut << ")" << std::endl; }
      if( ! ( *sigma > 0.0 ) ) { fatal_error() << "compute_local_entropy: sigma must be > 0 (got " << *sigma << ")" << std::endl; }
      // the neighbor list is sized from the global rcut_max at initialization (nbh_dist), before analysis operators run :
      // raising rcut_max here would come too late, so neighbors between rcut_max and rcut would be silently missing
      if( *rcut > *rcut_max )
      {
        fatal_error() << "compute_local_entropy: rcut=" << *rcut << " exceeds the neighbor list cutoff rcut_max=" << *rcut_max
                      << ", neighbors would be missing. Set 'rcut_max' >= rcut in the global section of the input file" << std::endl;
      }

      const size_t n_bins = ( *nbins > 0 ) ? *nbins : ( static_cast<size_t>( (*rcut) / (*sigma) ) + 1 );
      if( n_bins < 2 || n_bins > LOCAL_ENTROPY_MAX_BINS )
      {
        fatal_error() << "compute_local_entropy: nbins=" << n_bins << " out of range [2;" << LOCAL_ENTROPY_MAX_BINS
                      << "] (set nbins explicitly, or use a larger sigma / smaller rcut)" << std::endl;
      }

      const double volume = thermodynamic_state->volume();
      const size_t particle_count = thermodynamic_state->particle_count();
      if( ! ( volume > 0.0 ) || particle_count == 0 )
      {
        fatal_error() << "compute_local_entropy: invalid thermodynamic_state (volume=" << volume << ", particle_count=" << particle_count
                      << "), compute it before this operator" << std::endl;
      }

      const double drm = (*rcut) / ( n_bins - 1 );
      const double rho_global = static_cast<double>(particle_count) / volume;
      const double inv_sigmasq2 = 1.0 / ( 2.0 * (*sigma) * (*sigma) );

      // half width of the gaussian window, in bins : first bin where exp(-r^2/2sigma^2) < 1e-5
      size_t deltabin = 0;
      while( deltabin < n_bins && std::exp( - ( deltabin * deltabin * drm * drm ) * inv_sigmasq2 ) >= 1.0e-5 ) ++ deltabin;

      if( *verbose )
      {
        lout << onika::format_string("compute_local_entropy: rcut=%.5f sigma=%.5f nbins=%d deltabin=%d local=%s\n"
                                    , *rcut , *sigma , int(n_bins) , int(deltabin) , (*local) ? "true" : "false" );
      }

      auto entropy_acc = grid->field_accessor( field::mk_generic_real( *name ) );

      using ComputeBuffer = ComputePairBuffer2<false,false,LocalEntropyExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      LocalEntropyFunctor<decltype(entropy_acc)> compute_op =
        { n_bins , deltabin , drm , 1.0 / drm , inv_sigmasq2
        , 1.0 / ( 4.0 * M_PI * rho_global * std::sqrt( 2.0 * M_PI ) * (*sigma) )
        , rho_global , (4.0/3.0) * M_PI * (*rcut) * (*rcut) * (*rcut) , *local , entropy_acc };

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

Computes the per-atom two-body excess entropy (Piaget & Pablo, J. Chem. Phys. 147, 2017,
same definition as LAMMPS compute entropy/atom without averaging) as a per-particle scalar
field, GPU-compatibility to be tested:

  s_i = -2 pi rho integral_0^rcut [ g_i(r) ln g_i(r) - g_i(r) + 1 ] r^2 dr

where g_i(r) is the gaussian-smeared (width sigma) radial distribution of the neighbors of
atom i within rcut, and rho the global density (or the local one, neighbor count over the
rcut sphere volume, if local is true). The integral uses the trapezoid rule over nbins bins.

Requirements / limits:
  - thermodynamic_state must be up to date (global density = particle count / volume)
  - nbins (given, or rcut/sigma+1 when 0) must be in [2;128]
  - rcut must not exceed the neighbor list cutoff: set 'rcut_max: <value>' (>= rcut) in the
    global section when rcut is larger than the potential cutoff

Usage example:

compute_local_entropy:
  rcut: 5.0 ang
  sigma: 0.25 ang
  local: false
  name: entropy

)EOF";
    }
  };

  template<class GridT> using ComputeLocalEntropyTmpl = ComputeLocalEntropy<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_local_entropy)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_local_entropy", make_grid_variant_operator< ComputeLocalEntropyTmpl > );
  }

}

namespace exanb
{
  // specialize functor traits to allow Cuda execution space
  template<class EntropyFieldT>
  struct ComputePairTraits< exaStamp::LocalEntropyFunctor<EntropyFieldT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };
}
