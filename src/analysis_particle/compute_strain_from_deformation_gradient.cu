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
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particles.h>

#include <exaStamp/analysis_particle/basic_algebra_for_mechanics.h>

// Pointwise tensors derived from an already-computed deformation gradient field F
// (see compute_deformation_gradient_tensor). No neighbor list involved: F is read
// and the derived tensor(s) written back, one particle at a time, so this is plain
// per-particle work and trivially GPU-compatible. To add another derived quantity,
// clone one of the two functors/operators below with the relevant closed-form
// expression of F.
namespace exaStamp
{
  using namespace exanb;

  struct GreenLagrangeStrainFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& F, Mat3d& E ) const
    {
      E = 0.5 * ( AikBkj( transpose(F), F ) - make_identity_matrix() );
    }
  };

  struct PolarDecompositionFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& F, Mat3d& R, Mat3d& U ) const
    {
      RU_decomposition( F, R, U );
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::GreenLagrangeStrainFunctor> { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::PolarDecompositionFunctor>  { static inline constexpr bool CudaCompatible = true; };
}

namespace exaStamp
{
  template<class GridT>
  class ComputeGreenLagrangeStrain : public OperatorNode
  {
    ADD_SLOT( GridT       , grid          , INPUT_OUTPUT );
    ADD_SLOT( std::string , defgrad_field , INPUT , std::string("defgrad")        , DocString{"Name of the input deformation gradient tensor field"} );
    ADD_SLOT( std::string , strain_field  , INPUT , std::string("green_lagrange") , DocString{"Name of the resulting Green-Lagrange strain tensor field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *defgrad_field ) ) )
      {
        fatal_error() << "compute_green_lagrange_strain: input field '" << *defgrad_field << "' does not exist (run compute_deformation_gradient_tensor first, or check defgrad_field)" << std::endl;
      }
      auto F_acc = grid->field_const_accessor( field::mk_generic_mat3( *defgrad_field ) );
      auto E_acc = grid->field_accessor( field::mk_generic_mat3( *strain_field ) );
      compute_cell_particles( *grid, false, GreenLagrangeStrainFunctor{}, onika::make_flat_tuple( F_acc, E_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the Green-Lagrange strain tensor E = 1/2 (F^T F - I) per particle, from an
already-computed deformation gradient tensor field (see compute_deformation_gradient_tensor).
Pointwise, no neighbor list needed.

Usage example:

compute_green_lagrange_strain:
  defgrad_field: defgrad
  strain_field: green_lagrange

)EOF";
    }
  };

  template<class GridT>
  class ComputePolarDecomposition : public OperatorNode
  {
    ADD_SLOT( GridT       , grid          , INPUT_OUTPUT );
    ADD_SLOT( std::string , defgrad_field , INPUT , std::string("defgrad")  , DocString{"Name of the input deformation gradient tensor field"} );
    ADD_SLOT( std::string , rot_field     , INPUT , std::string("rotation") , DocString{"Name of the resulting pure rotation tensor field"} );
    ADD_SLOT( std::string , stretch_field , INPUT , std::string("stretch")  , DocString{"Name of the resulting pure (right) stretch tensor field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *defgrad_field ) ) )
      {
        fatal_error() << "compute_polar_decomposition: input field '" << *defgrad_field << "' does not exist (run compute_deformation_gradient_tensor first, or check defgrad_field)" << std::endl;
      }
      auto F_acc = grid->field_const_accessor( field::mk_generic_mat3( *defgrad_field ) );
      auto R_acc = grid->field_accessor( field::mk_generic_mat3( *rot_field ) );
      auto U_acc = grid->field_accessor( field::mk_generic_mat3( *stretch_field ) );
      compute_cell_particles( *grid, false, PolarDecompositionFunctor{}, onika::make_flat_tuple( F_acc, R_acc, U_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the polar decomposition F = R U (pure rotation R, pure right stretch tensor U)
per particle, from an already-computed deformation gradient tensor field (see
compute_deformation_gradient_tensor). Pointwise, no neighbor list needed.

Usage example:

compute_polar_decomposition:
  defgrad_field: defgrad
  rot_field: rotation
  stretch_field: stretch

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_strain_from_deformation_gradient)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_green_lagrange_strain", make_grid_variant_operator< ComputeGreenLagrangeStrain > );
    OperatorNodeFactory::instance()->register_factory( "compute_polar_decomposition", make_grid_variant_operator< ComputePolarDecomposition > );
  }

}
