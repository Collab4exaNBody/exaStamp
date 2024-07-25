// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

// #pragma xstamp_grid_variant // DO NOT REMOVE THIS LINE

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>

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
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep, Mat3d& virial ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type];
      virial = m_type_virial_correction[type];
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type];
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep, Vec3d& couple, Mat3d& virial ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type];
      couple = Vec3d { 0.0 , 0.0 , 0.0 };
      virial = m_type_virial_correction[type];
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type, double& fx, double& fy, double& fz, double& ep, Vec3d& couple) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = m_type_energy_correction[type];
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
    ADD_SLOT( bool  , ghost  , INPUT , false );

  public:

    inline void execute () override final
    {
      molecule_compute_parameters->m_energy_correction.resize( species->size() , 0.0 );
      molecule_compute_parameters->m_virial_correction.resize( species->size() , Mat3d{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0} );

      auto type = grid->field_accessor( field::type );
      auto fx = grid->field_accessor( field::fx );
      auto fy = grid->field_accessor( field::fy );
      auto fz = grid->field_accessor( field::fz );
      auto ep = grid->field_accessor( field::ep );

      InitForceEnergyWithCorrection func = { molecule_compute_parameters->m_energy_correction.data() , molecule_compute_parameters->m_virial_correction.data() };

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
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "zero_force_energy_with_correction", make_grid_variant_operator< InitVirialEnergyCorrection > );
  }

}

