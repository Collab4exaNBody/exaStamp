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
#include <exaStamp/analysis_particle/strain_functional/sfd_math.h>

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// Strain Functional Descriptors (Kober, Tavenner, Adams, Mathew, arXiv:2402.04191):
// rotational invariants of the Gaussian-weighted neighbor position moments up to order 4.
// The math lives in sfd_math.h.
namespace exaStamp
{
  using namespace exanb;

  // per-particle running sums of w * x^p y^q z^r (p+q+r <= 4), positions in units of sigma
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) StrainFunctionalExtStorage
  {
    double m_raw[ sfd::N_MONOMIALS ] = {};

    ONIKA_HOST_DEVICE_FUNC
    inline void reset()
    {
      for(int i=0;i<sfd::N_MONOMIALS;i++) m_raw[i] = 0.0;
      m_raw[0] = 1.0; // central particle: weight 1 at the origin
    }
  };

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) StrainFunctionalFunctor
  {
    const double m_rcut2 = 0.0;
    const double m_inv_sigma = 1.0;
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    double * const __restrict__ m_invariants = nullptr; // N_INVARIANTS per particle
    double * const __restrict__ m_spherical = nullptr;  // N_SPHERICAL per particle, or nullptr

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, LocalCellsT, size_t, size_t, exanb::ComputePairParticleContextStart) const
    {
      ctx.ext.reset();
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
        ComputeBufferT& ctx
      , const Vec3d& dr, double d2
      , LocalCellsT, size_t, size_t
      , double ) const
    {
      if( d2 >= m_rcut2 ) return;
      // r_ij = r_i - r_j = -dr, in units of sigma
      const double x = -dr.x * m_inv_sigma, y = -dr.y * m_inv_sigma, z = -dr.z * m_inv_sigma;
      const double w = exp( -0.5 * ( x*x + y*y + z*z ) );
      sfd::accumulate_monomials( ctx.ext.m_raw , w , x , y , z );
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, LocalCellsT, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop) const
    {
      const size_t i = m_cell_particle_offset[cell_a] + p_a;
      const sfd::SFDSpherical s = sfd::spherical_decomposition( ctx.ext.m_raw );
      sfd::invariants( s , m_invariants + i * sfd::N_INVARIANTS );
      if( m_spherical != nullptr ) sfd::spherical_values( s , m_spherical + i * sfd::N_SPHERICAL );
    }
  };

  // one column of a flat per-particle buffer -> grid field
  struct StrainFunctionalColumnFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const double * const __restrict__ m_values = nullptr;
    const size_t m_stride = 0;
    const size_t m_column = 0;

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t cell, unsigned int part, double& out ) const
    {
      out = m_values[ ( m_cell_particle_offset[cell] + part ) * m_stride + m_column ];
    }
  };

}

namespace exanb
{
  template<>
  struct ComputePairTraits< exaStamp::StrainFunctionalFunctor >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };

  template<>
  struct ComputeCellParticlesTraits< exaStamp::StrainFunctionalColumnFunctor >
  {
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeStrainFunctional : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;

    ADD_SLOT( GridT                     , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain            , INPUT , REQUIRED );
    ADD_SLOT( double                    , rcut              , INPUT , REQUIRED , DocString{"Neighbor cutoff distance"} );
    ADD_SLOT( double                    , sigma             , INPUT , REQUIRED , DocString{"Width of the Gaussian neighbor weight exp(-r^2/(2 sigma^2)); moments are expressed in units of sigma"} );
    ADD_SLOT( bool                      , spherical_tensors , INPUT , false , DocString{"Also store the 35 lab-frame spherical tensor components in sfd_spherical_tensors"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double                    , rcut_max          , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( DoubleVector              , sfd_invariants    , OUTPUT , DocString{"Per-particle invariants and orientation factors, 36 values per particle at (cell_particle_offset[cell] + particle) * 36, see strain_functional_fields for the column names"} );
    ADD_SLOT( DoubleVector              , sfd_spherical_tensors , OUTPUT , DocString{"Per-particle spherical tensor components, 35 values per particle, same indexing; empty unless spherical_tensors is true"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( *sigma <= 0.0 ) { fatal_error() << "compute_strain_functional: sigma must be > 0" << std::endl; }

      const size_t n_particles = grid->number_of_particles();
      // ghost entries are not computed and stay zero
      sfd_invariants->assign( n_particles * sfd::N_INVARIANTS , 0.0 );
      if( *spherical_tensors ) sfd_spherical_tensors->assign( n_particles * sfd::N_SPHERICAL , 0.0 );
      else sfd_spherical_tensors->clear();
      if( grid->number_of_cells() == 0 ) return;

      using ComputeBuffer = ComputePairBuffer2<false,false,StrainFunctionalExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      StrainFunctionalFunctor compute_op = { (*rcut)*(*rcut) , 1.0/(*sigma) , grid->cell_particle_offset_data()
                                           , sfd_invariants->data() , (*spherical_tensors) ? sfd_spherical_tensors->data() : nullptr };

      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr std::integral_constant<bool,true> force_use_cells_accessor = {};
      static constexpr DefaultPositionFields posfields = {};
      compute_cell_particle_pairs2( *grid, *rcut, false, optional, compute_buf, compute_op, compute_field_set
                                   , posfields, parallel_execution_context(), force_use_cells_accessor );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Strain Functional Descriptors (E. M. Kober, J. P. Tavenner, C. M. Adams, N. Mathew,
arXiv:2402.04191). For each particle, the neighbors within rcut (and the particle itself)
are weighted by exp(-r^2/(2 sigma^2)). The weighted position moments up to order 4,
taken about the weighted centroid and in units of sigma, are split into SO(3)
irreducible spherical tensors. Their rotational invariants (1, 3, 7 and 12 values for
orders 1 to 4) and lab-frame orientation factors (3 per order) are written to the flat
per-particle buffer sfd_invariants (36 values per particle, the first one being the
total weight w). With spherical_tensors: true, the 35 spherical tensor components are
also written to sfd_spherical_tensors. GPU compatible.

Use strain_functional_fields to copy columns into per-particle grid fields.

Usage example:

compute_strain_functional:
  rcut: 5.1 ang
  sigma: 0.91 ang
  spherical_tensors: false

)EOF";
    }
  };

  template<class GridT>
  class StrainFunctionalFields : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;
    using StringVector = std::vector<std::string>;

    ADD_SLOT( GridT        , grid                  , INPUT_OUTPUT );
    ADD_SLOT( DoubleVector , sfd_invariants        , INPUT , REQUIRED );
    ADD_SLOT( DoubleVector , sfd_spherical_tensors , INPUT , DoubleVector{} );
    ADD_SLOT( StringVector , fields                , INPUT , StringVector{} , DocString{"Columns to copy (e.g. [ w, P2_I0, P4_I0, P22_0 ]), empty for all 36 invariant columns"} );
    ADD_SLOT( std::string  , prefix                , INPUT , std::string("sfd_") , DocString{"Grid field name = prefix + column name"} );

  public:
    inline void execute () override final
    {
      const char* const* inv_names = sfd::sfd_invariant_names();
      const char* const* st_names = sfd::sfd_spherical_names();
      StringVector columns = *fields;
      if( columns.empty() ) columns.assign( inv_names , inv_names + sfd::N_INVARIANTS );

      const size_t n_particles = grid->number_of_particles();
      for(const auto& name : columns)
      {
        const DoubleVector* values = nullptr;
        size_t stride = 0, column = 0;
        const auto inv_it = std::find_if( inv_names , inv_names + sfd::N_INVARIANTS , [&name](const char* s){ return name == s; } );
        const auto st_it = std::find_if( st_names , st_names + sfd::N_SPHERICAL , [&name](const char* s){ return name == s; } );
        if( inv_it != inv_names + sfd::N_INVARIANTS ) { values = &(*sfd_invariants); stride = sfd::N_INVARIANTS; column = inv_it - inv_names; }
        else if( st_it != st_names + sfd::N_SPHERICAL ) { values = &(*sfd_spherical_tensors); stride = sfd::N_SPHERICAL; column = st_it - st_names; }
        else { fatal_error() << "strain_functional_fields: unknown column '" << name << "'" << std::endl; }
        if( values->size() != n_particles * stride )
        {
          fatal_error() << "strain_functional_fields: column '" << name << "' not available"
                        << ( stride == sfd::N_SPHERICAL ? " (set spherical_tensors: true in compute_strain_functional)" : "" ) << std::endl;
        }

        auto acc = grid->field_accessor( field::mk_generic_real( (*prefix) + name ) );
        StrainFunctionalColumnFunctor func = { grid->cell_particle_offset_data() , values->data() , stride , column };
        compute_cell_particles( *grid , false , func , onika::make_flat_tuple( acc ) , parallel_execution_context() );
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Copies columns of compute_strain_functional's buffers into per-particle grid fields named
prefix + column (default prefix sfd_). Owned particles only: run ghost_update_opt on the
fields if ghosts need them.

Invariant columns (sfd_invariants): w, P1_I0, P2_I0..P2_I2, P3_I0..P3_I6, P4_I0..P4_I11,
O1_I0..O1_I2, O2_I0..O2_I2, O3_I0..O3_I2, O4_I0..O4_I2.
Spherical tensor columns (sfd_spherical_tensors, needs spherical_tensors: true):
P00_0, P11_0, P11_1r, P11_1i, P20_0, P22_0, P22_1r, P22_1i, P22_2r, P22_2i,
P31a_0, P31a_1r, P31a_1i, P33a_0, P33a_1r ... P33a_3i, P40_0, P42a_0 ... P42a_2i,
P44a_0, P44a_1r ... P44a_4i.

Usage example:

compute_strain_functional: { rcut: 5.1 ang, sigma: 0.91 ang }
strain_functional_fields: { fields: [ w, P2_I0, P3_I0, P4_I0 ] }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_strain_functional)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_strain_functional", make_grid_variant_operator< ComputeStrainFunctional > );
    OperatorNodeFactory::instance()->register_factory( "strain_functional_fields", make_grid_variant_operator< StrainFunctionalFields > );
  }

}
