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

  // axial vector of the skew part of a rotation tensor (e.g. R from compute_polar_decomposition)
  struct MicrorotationFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& R, Vec3d& mu ) const
    {
      const Mat3d skewR = 0.5 * ( R - transpose(R) );
      mu.x = 0.5 * ( skewR.m32 - skewR.m23 );
      mu.y = 0.5 * ( skewR.m13 - skewR.m31 );
      mu.z = 0.5 * ( skewR.m21 - skewR.m12 );
    }
  };

  ONIKA_HOST_DEVICE_FUNC inline double slip_tripod_sign( double x ) { return ( x >= 0.0 ) ? 1.0 : -1.0; }

  // local slip-plane basis tripod (Burgers-parallel l, Burgers-orthogonal m, glide-plane normal n)
  // from A = F - I, same construction as compute_local_mechanical_metrics.cpp's dislocation analysis
  struct SlipTripodFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& F, Vec3d& l, Vec3d& m, Vec3d& n ) const
    {
      const Mat3d A = F - make_identity_matrix();
      const Mat3d Ltens = A * transpose(A);
      const Mat3d Ntens = transpose(A) * A;
      constexpr double tol = 1.e-4;

      Vec3d lvec;
      if( Ltens.m11 > Ltens.m22 && Ltens.m11 > Ltens.m33 )
      {
        lvec.x = sqrt(Ltens.m11);
        lvec.y = slip_tripod_sign(Ltens.m12) * sqrt(Ltens.m22);
        lvec.z = slip_tripod_sign(Ltens.m13) * sqrt(Ltens.m33);
      }
      else if( Ltens.m22 > Ltens.m11 && Ltens.m22 > Ltens.m33 )
      {
        lvec.x = slip_tripod_sign(Ltens.m12) * sqrt(Ltens.m11);
        lvec.y = sqrt(Ltens.m22);
        lvec.z = slip_tripod_sign(Ltens.m23) * sqrt(Ltens.m33);
      }
      else
      {
        lvec.x = slip_tripod_sign(Ltens.m13) * sqrt(Ltens.m11);
        lvec.y = slip_tripod_sign(Ltens.m23) * sqrt(Ltens.m22);
        lvec.z = sqrt(Ltens.m33);
      }
      if( fabs(lvec.z) > tol )      { lvec = lvec * slip_tripod_sign(lvec.z); }
      else if( fabs(lvec.y) > tol ) { lvec = lvec * slip_tripod_sign(lvec.y); }

      Vec3d nvec;
      if( Ntens.m11 > Ntens.m22 && Ntens.m11 > Ntens.m33 )
      {
        nvec.x = sqrt(Ntens.m11);
        nvec.y = slip_tripod_sign(Ntens.m12) * sqrt(Ntens.m22);
        nvec.z = slip_tripod_sign(Ntens.m13) * sqrt(Ntens.m33);
      }
      else if( Ntens.m22 > Ntens.m11 && Ntens.m22 > Ntens.m33 )
      {
        nvec.x = slip_tripod_sign(Ntens.m12) * sqrt(Ntens.m11);
        nvec.y = sqrt(Ntens.m22);
        nvec.z = slip_tripod_sign(Ntens.m23) * sqrt(Ntens.m33);
      }
      else
      {
        nvec.x = slip_tripod_sign(Ntens.m13) * sqrt(Ntens.m11);
        nvec.y = slip_tripod_sign(Ntens.m23) * sqrt(Ntens.m22);
        nvec.z = sqrt(Ntens.m33);
      }
      if( fabs(nvec.z) > tol )      { nvec = nvec * slip_tripod_sign(nvec.z); }
      else if( fabs(nvec.y) > tol ) { nvec = nvec * slip_tripod_sign(nvec.y); }

      l = lvec;
      n = nvec;
      m = cross( nvec, lvec );
    }
  };

  // relative local volume change J = det(F) : J=1 no change, J>1 expansion, J<1 compaction
  struct JacobianFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& F, double& J ) const
    {
      J = determinant(F);
    }
  };

  // principal invariants of a symmetric tensor (e.g. the Green-Lagrange strain E),
  // basis-independent: I1=tr(T), I2=sum of principal 2x2 minors, I3=det(T)
  struct TensorInvariantsFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& T, double& I1, double& I2, double& I3 ) const
    {
      I1 = T.m11 + T.m22 + T.m33;
      I2 = ( T.m11*T.m22 - T.m12*T.m21 ) + ( T.m22*T.m33 - T.m23*T.m32 ) + ( T.m11*T.m33 - T.m13*T.m31 );
      I3 = determinant(T);
    }
  };

  // von Mises equivalent of a symmetric tensor, not itself one of the 3 principal
  // invariants (it's a combination of them via the deviatoric part) but the usual
  // single scalar for visualizing shear/distortion intensity. Same formula already
  // used for stress in exaStamp/thermo_state/thermodynamic_state.h's vonmises_scal(),
  // just without that function's volume normalization (T here is already a per-particle
  // tensor, not a virial).
  struct VonMisesFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& T, double& vm ) const
    {
      vm = sqrt( 0.5 * (
                    ( T.m11 - T.m22 ) * ( T.m11 - T.m22 )
                  + ( T.m22 - T.m33 ) * ( T.m22 - T.m33 )
                  + ( T.m33 - T.m11 ) * ( T.m33 - T.m11 )
                  + 6.0 * ( T.m12*T.m12 + T.m13*T.m13 + T.m23*T.m23 )
                  ) );
    }
  };

  // OVITO's shear strain measure (atomic strain / Polyhedral Template Matching
  // convention): proportional to von Mises above by a constant factor sqrt(3)
  // (vm = sqrt(3) * shear_strain), kept as its own operator since it's the
  // convention people comparing against OVITO output will expect bit-for-bit.
  struct ShearStrainFunctor
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Mat3d& T, double& shear ) const
    {
      shear = sqrt(
                    T.m12*T.m12 + T.m13*T.m13 + T.m23*T.m23
                  + (1.0/6.0) * (
                      ( T.m11 - T.m22 ) * ( T.m11 - T.m22 )
                    + ( T.m22 - T.m33 ) * ( T.m22 - T.m33 )
                    + ( T.m33 - T.m11 ) * ( T.m33 - T.m11 )
                    )
                  );
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::GreenLagrangeStrainFunctor> { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::PolarDecompositionFunctor>  { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::MicrorotationFunctor>       { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::SlipTripodFunctor>          { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::JacobianFunctor>            { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::TensorInvariantsFunctor>    { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::VonMisesFunctor>            { static inline constexpr bool CudaCompatible = true; };
  template<> struct ComputeCellParticlesTraits<exaStamp::ShearStrainFunctor>         { static inline constexpr bool CudaCompatible = true; };
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

  template<class GridT>
  class ComputeMicrorotation : public OperatorNode
  {
    ADD_SLOT( GridT       , grid           , INPUT_OUTPUT );
    ADD_SLOT( std::string , rot_field      , INPUT , std::string("rotation")      , DocString{"Name of the input pure rotation tensor field (see compute_polar_decomposition)"} );
    ADD_SLOT( std::string , microrot_field , INPUT , std::string("microrotation") , DocString{"Name of the resulting microrotation vector field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *rot_field ) ) )
      {
        fatal_error() << "compute_microrotation: input field '" << *rot_field << "' does not exist (run compute_polar_decomposition first, or check rot_field)" << std::endl;
      }
      auto R_acc  = grid->field_const_accessor( field::mk_generic_mat3( *rot_field ) );
      auto mu_acc = grid->field_accessor( field::mk_generic_vec3( *microrot_field ) );
      compute_cell_particles( *grid, false, MicrorotationFunctor{}, onika::make_flat_tuple( R_acc, mu_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the microrotation vector mu (axial vector of the skew part of the pure
rotation tensor R) per particle, from an already-computed rotation tensor field
(see compute_polar_decomposition). Pointwise, no neighbor list needed.

Usage example:

compute_microrotation:
  rot_field: rotation
  microrot_field: microrotation

)EOF";
    }
  };

  template<class GridT>
  class ComputeSlipTripod : public OperatorNode
  {
    ADD_SLOT( GridT       , grid              , INPUT_OUTPUT );
    ADD_SLOT( std::string , defgrad_field     , INPUT , std::string("defgrad")     , DocString{"Name of the input deformation gradient tensor field"} );
    ADD_SLOT( std::string , burgerpar_field   , INPUT , std::string("burgerpar")   , DocString{"Name of the resulting Burgers-parallel slip basis vector field"} );
    ADD_SLOT( std::string , burgerortho_field , INPUT , std::string("burgerortho") , DocString{"Name of the resulting Burgers-orthogonal slip basis vector field"} );
    ADD_SLOT( std::string , glide_field       , INPUT , std::string("glide")       , DocString{"Name of the resulting glide-plane basis vector field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *defgrad_field ) ) )
      {
        fatal_error() << "compute_slip_tripod: input field '" << *defgrad_field << "' does not exist (run compute_deformation_gradient_tensor first, or check defgrad_field)" << std::endl;
      }
      auto F_acc = grid->field_const_accessor( field::mk_generic_mat3( *defgrad_field ) );
      auto l_acc = grid->field_accessor( field::mk_generic_vec3( *burgerpar_field ) );
      auto m_acc = grid->field_accessor( field::mk_generic_vec3( *burgerortho_field ) );
      auto n_acc = grid->field_accessor( field::mk_generic_vec3( *glide_field ) );
      compute_cell_particles( *grid, false, SlipTripodFunctor{}, onika::make_flat_tuple( F_acc, l_acc, m_acc, n_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the local slip-plane basis tripod (Burgers-parallel l, Burgers-orthogonal m,
glide-plane normal n) per particle, from an already-computed deformation gradient
tensor field (see compute_deformation_gradient_tensor). Pointwise, no neighbor list
needed. Same construction as compute_local_mechanical_metrics.cpp's dislocation
analysis tripod (l,m,n), built from A = F - I.

Usage example:

compute_slip_tripod:
  defgrad_field: defgrad
  burgerpar_field: burgerpar
  burgerortho_field: burgerortho
  glide_field: glide

)EOF";
    }
  };

  template<class GridT>
  class ComputeJacobian : public OperatorNode
  {
    ADD_SLOT( GridT       , grid           , INPUT_OUTPUT );
    ADD_SLOT( std::string , defgrad_field  , INPUT , std::string("defgrad") , DocString{"Name of the input deformation gradient tensor field"} );
    ADD_SLOT( std::string , jacobian_field , INPUT , std::string("jacobian"), DocString{"Name of the resulting Jacobian (J=det(F)) scalar field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *defgrad_field ) ) )
      {
        fatal_error() << "compute_jacobian: input field '" << *defgrad_field << "' does not exist (run compute_deformation_gradient_tensor first, or check defgrad_field)" << std::endl;
      }
      auto F_acc = grid->field_const_accessor( field::mk_generic_mat3( *defgrad_field ) );
      auto J_acc = grid->field_accessor( field::mk_generic_real( *jacobian_field ) );
      compute_cell_particles( *grid, false, JacobianFunctor{}, onika::make_flat_tuple( F_acc, J_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the Jacobian J = det(F) per particle: the local relative volume change
(J=1 no change, J>1 expansion, J<1 compaction), from an already-computed deformation
gradient tensor field (see compute_deformation_gradient_tensor). Pointwise, no
neighbor list needed.

Usage example:

compute_jacobian:
  defgrad_field: defgrad
  jacobian_field: jacobian

)EOF";
    }
  };

  template<class GridT>
  class ComputeTensorInvariants : public OperatorNode
  {
    ADD_SLOT( GridT       , grid         , INPUT_OUTPUT );
    ADD_SLOT( std::string , tensor_field , INPUT , std::string("green_lagrange") , DocString{"Name of the input symmetric tensor field (e.g. see compute_green_lagrange_strain)"} );
    ADD_SLOT( std::string , i1_field     , INPUT , std::string("strain_i1") , DocString{"Name of the resulting 1st invariant (trace) field"} );
    ADD_SLOT( std::string , i2_field     , INPUT , std::string("strain_i2") , DocString{"Name of the resulting 2nd invariant (sum of principal minors) field"} );
    ADD_SLOT( std::string , i3_field     , INPUT , std::string("strain_i3") , DocString{"Name of the resulting 3rd invariant (determinant) field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *tensor_field ) ) )
      {
        fatal_error() << "compute_strain_invariants: input field '" << *tensor_field << "' does not exist (check tensor_field)" << std::endl;
      }
      auto T_acc  = grid->field_const_accessor( field::mk_generic_mat3( *tensor_field ) );
      auto I1_acc = grid->field_accessor( field::mk_generic_real( *i1_field ) );
      auto I2_acc = grid->field_accessor( field::mk_generic_real( *i2_field ) );
      auto I3_acc = grid->field_accessor( field::mk_generic_real( *i3_field ) );
      compute_cell_particles( *grid, false, TensorInvariantsFunctor{}, onika::make_flat_tuple( T_acc, I1_acc, I2_acc, I3_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the three principal (basis-independent) invariants of a symmetric tensor
field per particle: I1=tr(T), I2=sum of principal 2x2 minors, I3=det(T). Typically
applied to the Green-Lagrange strain tensor (see compute_green_lagrange_strain), but
works on any symmetric Mat3d field. Pointwise, no neighbor list needed.

Usage example:

compute_strain_invariants:
  tensor_field: green_lagrange
  i1_field: strain_i1
  i2_field: strain_i2
  i3_field: strain_i3

)EOF";
    }
  };

  template<class GridT>
  class ComputeVonMisesStrain : public OperatorNode
  {
    ADD_SLOT( GridT       , grid          , INPUT_OUTPUT );
    ADD_SLOT( std::string , tensor_field  , INPUT , std::string("green_lagrange")   , DocString{"Name of the input symmetric tensor field (e.g. see compute_green_lagrange_strain)"} );
    ADD_SLOT( std::string , vonmises_field, INPUT , std::string("von_mises") , DocString{"Name of the resulting von Mises equivalent scalar field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *tensor_field ) ) )
      {
        fatal_error() << "compute_von_mises_strain: input field '" << *tensor_field << "' does not exist (check tensor_field)" << std::endl;
      }
      auto T_acc  = grid->field_const_accessor( field::mk_generic_mat3( *tensor_field ) );
      auto vm_acc = grid->field_accessor( field::mk_generic_real( *vonmises_field ) );
      compute_cell_particles( *grid, false, VonMisesFunctor{}, onika::make_flat_tuple( T_acc, vm_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes the von Mises equivalent of a symmetric tensor field per particle: the usual
scalar summary of shear/distortion intensity, sqrt(1/2 [(T11-T22)^2+(T22-T33)^2+(T33-T11)^2
+ 6(T12^2+T13^2+T23^2)]). Typically applied to the Green-Lagrange strain tensor (see
compute_green_lagrange_strain), same formula convention already used for stress in
thermodynamic_state.h's vonmises_scal(). Pointwise, no neighbor list needed.

Usage example:

compute_von_mises_strain:
  tensor_field: green_lagrange
  vonmises_field: von_mises

)EOF";
    }
  };

  template<class GridT>
  class ComputeShearStrain : public OperatorNode
  {
    ADD_SLOT( GridT       , grid              , INPUT_OUTPUT );
    ADD_SLOT( std::string , tensor_field      , INPUT , std::string("green_lagrange") , DocString{"Name of the input symmetric tensor field (e.g. see compute_green_lagrange_strain)"} );
    ADD_SLOT( std::string , shear_strain_field, INPUT , std::string("shear_strain")   , DocString{"Name of the resulting shear strain scalar field"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( ! grid->has_allocated_field( field::mk_generic_mat3( *tensor_field ) ) )
      {
        fatal_error() << "compute_shear_strain: input field '" << *tensor_field << "' does not exist (check tensor_field)" << std::endl;
      }
      auto T_acc     = grid->field_const_accessor( field::mk_generic_mat3( *tensor_field ) );
      auto shear_acc = grid->field_accessor( field::mk_generic_real( *shear_strain_field ) );
      compute_cell_particles( *grid, false, ShearStrainFunctor{}, onika::make_flat_tuple( T_acc, shear_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes OVITO's shear strain measure of a symmetric tensor field per particle:
sqrt( T12^2+T13^2+T23^2 + 1/6 [(T11-T22)^2+(T22-T33)^2+(T33-T11)^2] ). Typically
applied to the Green-Lagrange strain tensor (see compute_green_lagrange_strain).
Proportional to compute_von_mises_strain's output by a constant factor sqrt(3)
(von Mises = sqrt(3) * shear_strain) but kept as a separate operator matching
OVITO's exact convention bit-for-bit. Pointwise, no neighbor list needed.

Usage example:

compute_shear_strain:
  tensor_field: green_lagrange
  shear_strain_field: shear_strain

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_strain_from_deformation_gradient)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_green_lagrange_strain", make_grid_variant_operator< ComputeGreenLagrangeStrain > );
    OperatorNodeFactory::instance()->register_factory( "compute_polar_decomposition", make_grid_variant_operator< ComputePolarDecomposition > );
    OperatorNodeFactory::instance()->register_factory( "compute_microrotation", make_grid_variant_operator< ComputeMicrorotation > );
    OperatorNodeFactory::instance()->register_factory( "compute_slip_tripod", make_grid_variant_operator< ComputeSlipTripod > );
    OperatorNodeFactory::instance()->register_factory( "compute_jacobian", make_grid_variant_operator< ComputeJacobian > );
    OperatorNodeFactory::instance()->register_factory( "compute_strain_invariants", make_grid_variant_operator< ComputeTensorInvariants > );
    OperatorNodeFactory::instance()->register_factory( "compute_von_mises_strain", make_grid_variant_operator< ComputeVonMisesStrain > );
    OperatorNodeFactory::instance()->register_factory( "compute_shear_strain", make_grid_variant_operator< ComputeShearStrain > );
  }

}
