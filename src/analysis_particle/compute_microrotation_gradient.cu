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
#include <onika/math/basic_types_yaml.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/position_long_term_backup.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <cmath>

// Spatial gradient of the microrotation vector mu across neighbors, matched between
// a reference configuration (grid_t0/backup_r_lt) and the current one -- first step
// of the dislocation-detection chain (Zimmerman et al.-style analysis, ported from
// compute_local_analysis_particle_metrics.cpp's RefGradientComputeOp). Unlike
// compute_deformation_gradient_tensor, this only ever works in the reference frame:
// mu itself is a current-configuration quantity already, we're just computing how it
// varies spatially over the *reference* neighborhood -- no current-frame position is
// needed at all, so there's no second (current-box) periodic-image consistency wrap.
namespace exaStamp
{
  using namespace exanb;

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) MicrorotationGradientExtStorage
  {
    Vec3d m_pos0 = {};
    Vec3d m_mu0 = {};
    Mat3d m_tensorA = {}; // sum of w * (dPosRef (x) dPosRef)
    Mat3d m_tensorB = {}; // sum of w * (dMu (x) dPosRef)

    ONIKA_HOST_DEVICE_FUNC
    inline void reset()
    {
      m_tensorA = Mat3d{};
      m_tensorB = Mat3d{};
    }
  };

  ONIKA_HOST_DEVICE_FUNC
  inline bool mat3d_has_nan_mg( const Mat3d& m )
  {
    return ! ( (m.m11==m.m11) && (m.m12==m.m12) && (m.m13==m.m13)
            && (m.m21==m.m21) && (m.m22==m.m22) && (m.m23==m.m23)
            && (m.m31==m.m31) && (m.m32==m.m32) && (m.m33==m.m33) );
  }

  template<class GridT, class MuFieldT, class VecGradFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) MicrorotationGradientFunctor
  {
    using CellsT = decltype( GridT{}.cells() );

    const double m_rcut_sq = 0.0;
    const double a0 = 1.0;
    const double a1 = 0.0;
    const double a2 = 0.0;
    const double a3 = 0.0;
    const Mat3d m_xform_t0 {};
    const Mat3d m_hh0 {}; // m_xform_t0 * lattice, reference configuration box
    const CellsT m_cells_t0 = nullptr;
    MuFieldT m_mu_field = {};
    VecGradFieldT m_vecgrad_field = {};

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart) const
    {
      ctx.ext.reset();
      Vec3d r0 { m_cells_t0[cell_a][field::rx][p_a] , m_cells_t0[cell_a][field::ry][p_a] , m_cells_t0[cell_a][field::rz][p_a] };
      ctx.ext.m_pos0 = m_xform_t0 * r0;
      ctx.ext.m_mu0 = cells[cell_a][m_mu_field][p_a];
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop) const
    {
      Mat3d vecgrad = AikBkj( ctx.ext.m_tensorB, inverse(ctx.ext.m_tensorA) );
      // spatial gradient: "no gradient" baseline is zero, not identity (same reasoning as L)
      if( mat3d_has_nan_mg(vecgrad) ) { vecgrad = Mat3d{}; }
      cells[cell_a][m_vecgrad_field][p_a] = vecgrad;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
       ComputeBufferT& ctx
      , const Vec3d& /*dr*/, double /*d2*/
      , LocalCellsT cells, size_t cell_b, size_t p_b
      , double /*scale*/) const
    {
      Vec3d r0b { m_cells_t0[cell_b][field::rx][p_b] , m_cells_t0[cell_b][field::ry][p_b] , m_cells_t0[cell_b][field::rz][p_b] };
      Vec3d deltaPosInit = ( m_xform_t0 * r0b ) - ctx.ext.m_pos0;

      // wrap to minimum image in the reference box (only frame that matters here)
      Vec3d deltaPosInitCris = inverse(m_hh0) * deltaPosInit;
      if( deltaPosInitCris.x > 0.5 ) { deltaPosInitCris.x -= 1.0; } else if( deltaPosInitCris.x < -0.5 ) { deltaPosInitCris.x += 1.0; }
      if( deltaPosInitCris.y > 0.5 ) { deltaPosInitCris.y -= 1.0; } else if( deltaPosInitCris.y < -0.5 ) { deltaPosInitCris.y += 1.0; }
      if( deltaPosInitCris.z > 0.5 ) { deltaPosInitCris.z -= 1.0; } else if( deltaPosInitCris.z < -0.5 ) { deltaPosInitCris.z += 1.0; }
      deltaPosInit = m_hh0 * deltaPosInitCris;

      const double rrInit2 = norm2(deltaPosInit);
      if( rrInit2 <= m_rcut_sq )
      {
        double w = a0 + a2*rrInit2;
        if( a1!=0.0 || a3!=0.0 ) { const double d = sqrt(rrInit2); w += a1*d + a3*rrInit2*d; }
        const Vec3d mu_b = cells[cell_b][m_mu_field][p_b];
        const Vec3d deltaMu = mu_b - ctx.ext.m_mu0;
        ctx.ext.m_tensorA += tensor(deltaPosInit,deltaPosInit) * w;
        ctx.ext.m_tensorB += tensor(deltaMu,deltaPosInit) * w;
      }
    }
  };
}

namespace exanb
{
  template<class GridT, class MuFieldT, class VecGradFieldT>
  struct ComputePairTraits< exaStamp::MicrorotationGradientFunctor<GridT,MuFieldT,VecGradFieldT> >
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
  class ComputeMicrorotationGradient : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;

    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT , DocString{"Local sub-domain particles grid"} );
    ADD_SLOT( GridT                     , grid_t0         , INPUT        , REQUIRED , DocString{"Reference configuration grid (same particle ordering as grid)"} );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED , DocString{"Simulation domain"} );
    ADD_SLOT( PositionLongTermBackup    , backup_r_lt     , INPUT        , REQUIRED , DocString{"Reference configuration position/xform backup (see backup_r_lt), used here only for its m_xform"} );
    ADD_SLOT( double                    , rcut            , INPUT        , REQUIRED , DocString{"Cutoff distance, in the reference configuration, for the neighbors contributing to the local microrotation gradient"} );
    ADD_SLOT( DoubleVector               , weight_function , INPUT        , DoubleVector{ {1.0} } , DocString{"List of [a0,...,an] coefficients for the polynomial distance weighting function : a0*x^0 + a1*x^1 + ... +an*x^n, applied to the reference-frame neighbor distance"} );
    ADD_SLOT( std::string               , microrot_field  , INPUT        , std::string("microrotation") , DocString{"Name of the input microrotation vector field (see compute_microrotation). Must have fresh ghost data -- run ghost_update_opt on it first."} );
    ADD_SLOT( std::string               , vecgrad_field   , INPUT        , std::string("vecgrad") , DocString{"Name of the resulting spatial gradient of the microrotation vector field"} );
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
      if( ! grid->has_allocated_field( field::mk_generic_vec3( *microrot_field ) ) )
      {
        fatal_error() << "compute_microrotation_gradient: input field '" << *microrot_field << "' does not exist (run compute_microrotation first, or check microrot_field)" << std::endl;
      }

      double poly_coefs[4] = { 1.0 , 0.0 , 0.0 , 0.0 };
      for(size_t i=0;i<weight_function->size() && i<4;i++) poly_coefs[i] = weight_function->at(i);

      const Mat3d xform_t0 = backup_r_lt->m_xform;
      const Mat3d lattice = diag_matrix( domain->extent() - domain->origin() );
      const Mat3d hh0 = xform_t0 * lattice;

      auto mu_acc = grid->field_accessor( field::mk_generic_vec3( *microrot_field ) );
      auto vecgrad_acc = grid->field_accessor( field::mk_generic_mat3( *vecgrad_field ) );

      using ComputeBuffer = ComputePairBuffer2<false,false,MicrorotationGradientExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      MicrorotationGradientFunctor<GridT,decltype(mu_acc),decltype(vecgrad_acc)> compute_op =
        { (*rcut)*(*rcut) , poly_coefs[0] , poly_coefs[1] , poly_coefs[2] , poly_coefs[3] , xform_t0 , hh0 , grid_t0->cells() , mu_acc , vecgrad_acc };

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

Computes the spatial gradient of the microrotation vector field (mu) per particle,
directly as a native field, GPU-compatible. Matches neighbors between a reference
configuration (grid_t0/backup_r_lt) and the current one, same as
compute_deformation_gradient_tensor, but the gradient itself is computed entirely in
the reference frame (mu is already a current-configuration quantity, we're only
measuring how it varies spatially -- no current-frame periodic-image consistency
step is needed). First step of the dislocation-detection chain (see
compute_dislocation_indicators for the second step). Requires the microrotation
input field's ghost data to be fresh -- run ghost_update_opt on it first, since
neighbor mu values (potentially on ghost particles) are read directly.

Usage example:

ghost_update_opt: { opt_fields: [ "microrotation" ] }
compute_microrotation_gradient:
  grid_t0: <reference configuration grid>
  backup_r_lt: <PositionLongTermBackup used to restore grid_t0>
  rcut: 8.0 ang
  weight_function: [ 1.0 , 0.0 , -0.01 ]
  microrot_field: microrotation
  vecgrad_field: vecgrad

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_microrotation_gradient)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_microrotation_gradient", make_grid_variant_operator< ComputeMicrorotationGradient > );
  }

}
