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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <cmath>

// Velocity gradient tensor L, computed the same way as compute_deformation_gradient_tensor
// computes F (weighted least-squares tensor fit over neighbors, polynomial distance weight),
// but entirely within the current configuration -- no reference grid/xform needed, since L is
// a rate (velocity difference vs. position difference), not a gradient of a coordinate map.
namespace exaStamp
{
  using namespace exanb;

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) VelocityGradientExtStorage
  {
    Vec3d m_v0 = {};
    Mat3d m_tensorA = {};
    Mat3d m_tensorB = {};

    ONIKA_HOST_DEVICE_FUNC
    inline void reset()
    {
      m_tensorA = Mat3d{};
      m_tensorB = Mat3d{};
    }
  };

  template<class VelGradFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) VelocityGradientFunctor
  {
    const double m_rcut_sq = 0.0;
    const double a0 = 1.0;
    const double a1 = 0.0;
    const double a2 = 0.0;
    const double a3 = 0.0;
    VelGradFieldT m_velgrad_field = {};

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart) const
    {
      ctx.ext.reset();
      ctx.ext.m_v0 = Vec3d{ cells[cell_a][field::vx][p_a] , cells[cell_a][field::vy][p_a] , cells[cell_a][field::vz][p_a] };
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop) const
    {
      Mat3d L = AikBkj( ctx.ext.m_tensorB, inverse(ctx.ext.m_tensorA) );
      // L is a rate tensor: its "nothing happening" value is zero, not identity (unlike F)
      const bool has_nan = ! ( (L.m11==L.m11) && (L.m12==L.m12) && (L.m13==L.m13)
                            && (L.m21==L.m21) && (L.m22==L.m22) && (L.m23==L.m23)
                            && (L.m31==L.m31) && (L.m32==L.m32) && (L.m33==L.m33) );
      if( has_nan ) { L = Mat3d{}; }
      cells[cell_a][m_velgrad_field][p_a] = L;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
       ComputeBufferT& ctx
      , const Vec3d& dr, double d2
      , LocalCellsT cells, size_t cell_b, size_t p_b
      , double /*scale*/) const
    {
      if( d2 <= m_rcut_sq )
      {
        double w = a0 + a2*d2;
        if( a1!=0.0 || a3!=0.0 ) { const double d = sqrt(d2); w += a1*d + a3*d2*d; }
        const Vec3d vb { cells[cell_b][field::vx][p_b] , cells[cell_b][field::vy][p_b] , cells[cell_b][field::vz][p_b] };
        const Vec3d dv = vb - ctx.ext.m_v0;
        ctx.ext.m_tensorA += tensor(dr,dr) * w;
        ctx.ext.m_tensorB += tensor(dv,dr) * w;
      }
    }
  };
}

namespace exanb
{
  template<class VelGradFieldT>
  struct ComputePairTraits< exaStamp::VelocityGradientFunctor<VelGradFieldT> >
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
  class ComputeVelocityGradientTensor : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;

    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT , DocString{"Local sub-domain particles grid"} );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED , DocString{"Simulation domain"} );
    ADD_SLOT( double                    , rcut            , INPUT        , REQUIRED , DocString{"Cutoff distance for the neighbors contributing to the local velocity gradient"} );
    ADD_SLOT( DoubleVector              , weight_function , INPUT        , DoubleVector{ {1.0} } , DocString{"List of [a0,...,an] coefficients for the polynomial distance weighting function : a0*x^0 + a1*x^1 + ... +an*x^n"} );
    ADD_SLOT( std::string               , velgrad_field   , INPUT        , std::string("velgrad") , DocString{"Name of the resulting velocity gradient tensor field"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT        , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double                    , rcut_max        , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( grid->number_of_cells() == 0 ) return;
      if( weight_function->size() > 4 )
      {
        fatal_error() << "weighting function polynomial has a maximum degree of 3 (maximum 4 coefficients)" << std::endl;
      }

      double poly_coefs[4] = { 1.0 , 0.0 , 0.0 , 0.0 };
      for(size_t i=0;i<weight_function->size() && i<4;i++) poly_coefs[i] = weight_function->at(i);

      auto velgrad_acc = grid->field_accessor( field::mk_generic_mat3( *velgrad_field ) );

      using ComputeBuffer = ComputePairBuffer2<false,false,VelocityGradientExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      VelocityGradientFunctor<decltype(velgrad_acc)> compute_op =
        { (*rcut)*(*rcut) , poly_coefs[0] , poly_coefs[1] , poly_coefs[2] , poly_coefs[3] , velgrad_acc };

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

Computes the velocity gradient tensor L per particle, directly as a native field,
GPU-compatible, from a weighted least-squares fit over neighbors in the current
configuration only (no reference grid needed, unlike compute_deformation_gradient_tensor,
since L is a rate: neighbor position differences vs. neighbor velocity differences).
The neighbor weighting is a polynomial of the current-frame neighbor distance, same
convention as average_neighbors_scalar / compute_deformation_gradient_tensor.

IMPORTANT: needs ghost particle velocities to be up to date. The default per-step
pipeline's fast path (update_particles_fast_body, config_move_particles.msp) only
calls ghost_update_r (positions only) -- ghost vx/vy/vz can be stale on any step that
didn't just do a full ghost update. Call ghost_update_r_v right before this operator
to be sure (verified: without it, results silently differ from the correct value by
~10% locally, no error or warning of any kind).

Usage example:

ghost_update_r_v
compute_velocity_gradient_tensor:
  rcut: 8.0 ang
  weight_function: [ 1.0 , 0.0 , -0.01 ]
  velgrad_field: velgrad

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_velocity_gradient_tensor)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_velocity_gradient_tensor", make_grid_variant_operator< ComputeVelocityGradientTensor > );
  }

}
