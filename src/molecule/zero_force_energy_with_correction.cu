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

//  // DO NOT REMOVE THIS LINE

//  // DO NOT REMOVE THIS LINE

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types_operators.h>

#include <onika/cuda/cuda.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exaStamp/molecule/molecule_compute_param.h>
#include <exaStamp/particle_species/particle_specie.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  struct InitForceEnergyWithCorrection
  {
    const double * __restrict__ m_type_energy_correction = nullptr;
    const Mat3d * __restrict__ m_type_virial_correction = nullptr;
    double m_inv_volume = 1.0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep, Mat3d& virial ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type] * m_inv_volume;
      virial = m_type_virial_correction[type] * m_inv_volume;
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type] * m_inv_volume;
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep, Vec3d& couple, Mat3d& virial ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type] * m_inv_volume;
      couple = Vec3d { 0.0 , 0.0 , 0.0 };
      virial = m_type_virial_correction[type] * m_inv_volume;
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep, Vec3d& couple) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type] * m_inv_volume;
      couple = Vec3d { 0.0 , 0.0 , 0.0 };
    }

  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::InitForceEnergyWithCorrection>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT>
  class InitVirialEnergyCorrection : public OperatorNode
  {  
    ADD_SLOT( GridT , grid  , INPUT_OUTPUT );
    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( ParticleSpecies              , species           , INPUT , REQUIRED );
    ADD_SLOT( Domain , domain , INPUT , REQUIRED , DocString{"Simulation domain, for the box volume"} );
    ADD_SLOT( bool  , ghost  , INPUT , false , DocString{"Also process ghost atoms"} );

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(
Resets forces before the force computation, and initializes the potential energy and the virial of each atom with the long
range correction of its type (computed by intramolecular_setup, zero without long_range_correction), divided by the
current box volume : the correction follows volume changes.
)EOF";
    }

    inline void execute () override final
    {
      molecule_compute_parameters->m_energy_correction.resize( species->size() , 0.0 );
      molecule_compute_parameters->m_virial_correction.resize( species->size() , Mat3d{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0} );

      auto type = grid->field_accessor( field::type );
      auto fx = grid->field_accessor( field::fx );
      auto fy = grid->field_accessor( field::fy );
      auto fz = grid->field_accessor( field::fz );
      auto ep = grid->field_accessor( field::ep );

      // the corrections are stored times the box volume
      double volume = 1.0;
      if( ! domain->xform_is_identity() )
      {
        const Mat3d mat = domain->xform();
        const Vec3d a { mat.m11, mat.m21, mat.m31 };
        const Vec3d b { mat.m12, mat.m22, mat.m32 };
        const Vec3d c { mat.m13, mat.m23, mat.m33 };
        volume = dot( cross(a,b) , c );
      }
      volume *= bounds_volume( domain->bounds() );

      InitForceEnergyWithCorrection func = { molecule_compute_parameters->m_energy_correction.data() , molecule_compute_parameters->m_virial_correction.data() , 1.0 / volume };

      if( grid->has_allocated_field( field::couple ) )
      {
        auto couple = grid->field_accessor( field::couple );        
        if( grid->has_allocated_field( field::virial ) )
        {
          auto virial = grid->field_accessor( field::virial );        
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(type,fx,fy,fz,ep,couple,virial) , parallel_execution_context() );
        }
        else
        {
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(type,fx,fy,fz,ep,couple) , parallel_execution_context() );
        }
      }
      else
      {
        if( grid->has_allocated_field( field::virial ) )
        {
          auto virial = grid->field_accessor( field::virial );        
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(type,fx,fy,fz,ep,virial) , parallel_execution_context() );
        }
        else
        {
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(type,fx,fy,fz,ep) , parallel_execution_context() );
        }
      }
    }

  };
    
 // === register factories ===  
  ONIKA_AUTORUN_INIT(zero_force_energy_with_correction)
  {
   OperatorNodeFactory::instance()->register_factory( "zero_force_energy_with_correction", make_grid_variant_operator< InitVirialEnergyCorrection > );
  }

}

