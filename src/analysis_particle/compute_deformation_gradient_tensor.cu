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

namespace exaStamp
{
  using namespace exanb;

  // Running per-particle accumulators : reference position of the central particle,
  // the two tensors whose ratio (BF . AF^-1) is the deformation gradient tensor, and
  // the (unweighted) sum/count of "slipped" neighbors for the slip vector (Zimmerman
  // et al., Phys. Rev. Lett. 87, 165507 (2001), Eq. 1).
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) DeformationGradientExtStorage
  {
    Vec3d m_pos0 = {};
    Mat3d m_tensorAF = {};
    Mat3d m_tensorBF = {};
    Vec3d m_slip_sum = {};
    long m_slip_count = 0;

    ONIKA_HOST_DEVICE_FUNC
    inline void reset()
    {
      m_tensorAF = Mat3d{};
      m_tensorBF = Mat3d{};
      m_slip_sum = Vec3d{};
      m_slip_count = 0;
    }
  };

  ONIKA_HOST_DEVICE_FUNC
  inline bool mat3d_has_nan(const Mat3d& m)
  {
    return ! ( (m.m11==m.m11) && (m.m12==m.m12) && (m.m13==m.m13)
            && (m.m21==m.m21) && (m.m22==m.m22) && (m.m23==m.m23)
            && (m.m31==m.m31) && (m.m32==m.m32) && (m.m33==m.m33) );
  }

  // computes, per particle, the deformation gradient tensor F between a reference
  // configuration (m_cells_t0/m_xform_t0) and the current one, as a weighted least
  // squares fit over neighbors matched between both configurations (same idea as
  // exanb's average_neighbors_scalar, but accumulating outer-product tensors instead
  // of a scalar sum, and with the neighbor weight evaluated on the reference distance
  // so the smoothing radius is unaffected by the current deformation).
  template<class GridT, class DefGradFieldT, class SlipFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) DeformationGradientFunctor
  {
    using CellsT = decltype( GridT{}.cells() );

    const double m_rcut_sq = 0.0;
    const double a0 = 1.0;
    const double a1 = 0.0;
    const double a2 = 0.0;
    const double a3 = 0.0;
    const Mat3d m_xform_t0 {};
    const Mat3d m_hh0 {}; // m_xform_t0 * lattice, reference configuration box
    const Mat3d m_hht {}; // m_xform    * lattice, current configuration box
    const CellsT m_cells_t0 = nullptr;
    DefGradFieldT m_defgrad_field = {};
    SlipFieldT m_slip_field = {};

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, LocalCellsT, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart) const
    {
      ctx.ext.reset();
      Vec3d r0 { m_cells_t0[cell_a][field::rx][p_a] , m_cells_t0[cell_a][field::ry][p_a] , m_cells_t0[cell_a][field::rz][p_a] };
      ctx.ext.m_pos0 = m_xform_t0 * r0;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, LocalCellsT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop) const
    {
      Mat3d F = AikBkj( ctx.ext.m_tensorBF, inverse(ctx.ext.m_tensorAF) );
      if( mat3d_has_nan(F) ) { F = make_identity_matrix(); }
      cells[cell_a][m_defgrad_field][p_a] = F;

      // slip vector (Zimmerman et al. PRL 87, 165507 (2001), Eq. 1) :
      // s = -(1/n_slipped) * sum over slipped neighbors of (x_cur - x_ref)
      const Vec3d slip = ( ctx.ext.m_slip_count > 0 ) ? ( ctx.ext.m_slip_sum * ( -1.0 / double(ctx.ext.m_slip_count) ) ) : Vec3d{};
      cells[cell_a][m_slip_field][p_a] = slip;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
       ComputeBufferT& ctx
      , const Vec3d& dr, double /*d2*/
      , LocalCellsT cells, size_t cell_b, size_t p_b
      , double /*scale*/) const
    {
      Vec3d r0b { m_cells_t0[cell_b][field::rx][p_b] , m_cells_t0[cell_b][field::ry][p_b] , m_cells_t0[cell_b][field::rz][p_b] };
      Vec3d deltaPosInit = ( m_xform_t0 * r0b ) - ctx.ext.m_pos0;
      Vec3d deltaPosCour = dr;

      // pick the same periodic image in both configurations : wrap deltaPosInit to
      // its minimum image in the reference box, then re-apply the same lattice
      // translation to deltaPosCour before wrapping it in the current (deformed) box
      Vec3d deltaPosInitCris = inverse(m_hh0) * deltaPosInit;
      Vec3d dtr {};
      if( deltaPosInitCris.x > 0.5 ) { deltaPosInitCris.x -= 1.0; dtr.x -= 1.0; } else if( deltaPosInitCris.x < -0.5 ) { deltaPosInitCris.x += 1.0; dtr.x += 1.0; }
      if( deltaPosInitCris.y > 0.5 ) { deltaPosInitCris.y -= 1.0; dtr.y -= 1.0; } else if( deltaPosInitCris.y < -0.5 ) { deltaPosInitCris.y += 1.0; dtr.y += 1.0; }
      if( deltaPosInitCris.z > 0.5 ) { deltaPosInitCris.z -= 1.0; dtr.z -= 1.0; } else if( deltaPosInitCris.z < -0.5 ) { deltaPosInitCris.z += 1.0; dtr.z += 1.0; }
      deltaPosInit = m_hh0 * deltaPosInitCris;

      Vec3d deltaPosCourCris = inverse(m_hht) * deltaPosCour + dtr;
      if( deltaPosCourCris.x > 0.5 ) { deltaPosCourCris.x -= 1.0; } else if( deltaPosCourCris.x < -0.5 ) { deltaPosCourCris.x += 1.0; }
      if( deltaPosCourCris.y > 0.5 ) { deltaPosCourCris.y -= 1.0; } else if( deltaPosCourCris.y < -0.5 ) { deltaPosCourCris.y += 1.0; }
      if( deltaPosCourCris.z > 0.5 ) { deltaPosCourCris.z -= 1.0; } else if( deltaPosCourCris.z < -0.5 ) { deltaPosCourCris.z += 1.0; }
      deltaPosCour = m_hht * deltaPosCourCris;

      const double rrInit2 = norm2(deltaPosInit);
      if( rrInit2 <= m_rcut_sq )
      {
        double w = a0 + a2*rrInit2;
        if( a1!=0.0 || a3!=0.0 ) { const double d = sqrt(rrInit2); w += a1*d + a3*rrInit2*d; }
        // scaling by 1e20 is a numerical conditioning trick only (cancels out in
        // AikBkj(tensorBF,inverse(tensorAF)), same as the original per-atom implementation)
        ctx.ext.m_tensorAF += tensor(deltaPosInit,deltaPosInit) * w * 1.0e20;
        ctx.ext.m_tensorBF += tensor(deltaPosCour,deltaPosInit) * w * 1.0e20;

        // slip vector accumulation is unweighted (no w factor), only counting
        // neighbors whose relative position changed by more than 1/10 of rcut
        // (filters out thermal-noise-level "slip")
        const Vec3d delta_cour_init = deltaPosCour - deltaPosInit;
        if( norm2(delta_cour_init) >= ( m_rcut_sq / 100.0 ) )
        {
          ctx.ext.m_slip_sum += delta_cour_init;
          ++ ctx.ext.m_slip_count;
        }
      }
    }
  };

  template<class GridT>
  class ComputeDeformationGradientTensor : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;

    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT , DocString{"Local sub-domain particles grid"} );
    ADD_SLOT( GridT                     , grid_t0         , INPUT        , REQUIRED , DocString{"Reference configuration grid (same particle ordering as grid)"} );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED , DocString{"Simulation domain"} );
    ADD_SLOT( PositionLongTermBackup    , backup_r_lt     , INPUT        , REQUIRED , DocString{"Reference configuration position/xform backup (see backup_r_lt), used here only for its m_xform"} );
    ADD_SLOT( double                    , rcut            , INPUT        , REQUIRED , DocString{"Cutoff distance, in the reference configuration, for the neighbors contributing to the local deformation gradient"} );
    ADD_SLOT( DoubleVector               , weight_function , INPUT        , DoubleVector{ {1.0} } , DocString{"List of [a0,...,an] coefficients for the polynomial distance weighting function : a0*x^0 + a1*x^1 + ... +an*x^n, applied to the reference-frame neighbor distance"} );
    ADD_SLOT( std::string               , defgrad_field   , INPUT        , std::string("defgrad") , DocString{"Name of the resulting per-particle deformation gradient tensor field"} );
    ADD_SLOT( std::string               , slip_field      , INPUT        , std::string("slip") , DocString{"Name of the resulting per-particle slip vector field (Zimmerman et al. PRL 87, 165507 (2001))"} );
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

      const Mat3d xform_t0 = backup_r_lt->m_xform;
      const Mat3d xform = domain->xform();
      const Mat3d lattice = diag_matrix( domain->extent() - domain->origin() );
      const Mat3d hh0 = xform_t0 * lattice;
      const Mat3d hht = xform * lattice;

      auto defgrad_acc = grid->field_accessor( field::mk_generic_mat3( *defgrad_field ) );
      auto slip_acc = grid->field_accessor( field::mk_generic_vec3( *slip_field ) );

      using ComputeBuffer = ComputePairBuffer2<false,false,DeformationGradientExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      DeformationGradientFunctor<GridT,decltype(defgrad_acc),decltype(slip_acc)> compute_op =
        { (*rcut)*(*rcut) , poly_coefs[0] , poly_coefs[1] , poly_coefs[2] , poly_coefs[3] , xform_t0 , hh0 , hht , grid_t0->cells() , defgrad_acc , slip_acc };

      LinearXForm cp_xform { xform };
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

Computes the local deformation gradient tensor (Mat3d) per particle, directly as a
native per-particle field, GPU-compatible. For each particle, matches its neighbors
between a reference configuration (grid_t0/backup_r_lt) and the current one, and fits F
as a weighted least-squares tensor ratio over that neighborhood. The neighbor weight
is a polynomial of the reference-frame distance, same convention as average_neighbors_scalar.
The reference xform is read from backup_r_lt->m_xform, so no separate xform backup
(e.g. backup_xform) is needed as long as grid_t0 was itself restored from that same backup.

Also computes the slip vector s (Zimmerman, Kelchner, Klein, Hamilton, Foiles,
Phys. Rev. Lett. 87, 165507 (2001), Eq. 1): s = -(1/n_slipped) * sum over slipped
neighbors of (x_current - x_reference), unweighted, only counting neighbors whose
relative position changed by at least 1/10 of rcut between the two configurations.
Reuses the exact same reference/current neighbor matching already computed for F,
at near-zero extra cost.

Usage example:

compute_deformation_gradient_tensor:
  grid_t0: <reference configuration grid>
  backup_r_lt: <PositionLongTermBackup used to restore grid_t0>
  rcut: 8.0 ang
  weight_function: [ 1.0 , 0.0 , -0.01 ] # => 1 + 0.0 r - 0.01 r^2, r being the reference-frame neighbor distance
  defgrad_field: defgrad
  slip_field: slip

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_deformation_gradient_tensor)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_deformation_gradient_tensor", make_grid_variant_operator< ComputeDeformationGradientTensor > );
  }

}

namespace exanb
{
  // specialize functor traits to allow Cuda execution space
  template<class GridT, class DefGradFieldT, class SlipFieldT>
  struct ComputePairTraits< exaStamp::DeformationGradientFunctor<GridT,DefGradFieldT,SlipFieldT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };
}
