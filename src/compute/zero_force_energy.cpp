// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

// #pragma xstamp_grid_variant // DO NOT REMOVE THIS LINE

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>

#include <onika/cuda/cuda.h>
#include <exanb/compute/compute_cell_particles.h>


namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  struct ZeroCellParticleFields
  {
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double& fx, double& fy, double& fz, double& ep, Mat3d& virial ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = 0.0;
      virial = Mat3d { 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 };
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double& fx, double& fy, double& fz, double& ep ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = 0.0;
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double& fx, double& fy, double& fz, double& ep, Vec3d& couple, Mat3d& virial ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = 0.0;
      couple = Vec3d { 0.0 , 0.0 , 0.0 };
      virial = Mat3d { 0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0 };
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double& fx, double& fy, double& fz, double& ep, Vec3d& couple) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      ep = 0.0;
      couple = Vec3d { 0.0 , 0.0 , 0.0 };
    }

  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::ZeroCellParticleFields>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT, class _FX, class _FY, class _FZ, class _EP, class _VIRIAL, class _COUPLE >
  class ZeroForceEnergyNode : public OperatorNode
  {  
    ADD_SLOT( GridT , grid  , INPUT_OUTPUT );
    ADD_SLOT( bool  , ghost  , INPUT , false );

  public:

    inline void execute () override final
    {
      auto fx = grid->field_accessor( onika::soatl::FieldId<_FX> {} );
      auto fy = grid->field_accessor( onika::soatl::FieldId<_FY> {} );
      auto fz = grid->field_accessor( onika::soatl::FieldId<_FZ> {} );
      auto ep = grid->field_accessor( onika::soatl::FieldId<_EP> {} );

      onika::soatl::FieldId<_VIRIAL> virial_field = {};
      onika::soatl::FieldId<_COUPLE> couple_field = {};

      ZeroCellParticleFields func = {};

      if( grid->has_allocated_field(couple_field) )
      {
        auto couple = grid->field_accessor( couple_field );        
        if( grid->has_allocated_field(virial_field) )
        {
          auto virial = grid->field_accessor( virial_field );        
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(fx,fy,fz,ep,couple,virial) , parallel_execution_context() );
        }
        else
        {
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(fx,fy,fz,ep,couple) , parallel_execution_context() );
        }
      }
      else
      {
        if( grid->has_allocated_field(virial_field) )
        {
          auto virial = grid->field_accessor( virial_field );        
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(fx,fy,fz,ep,virial) , parallel_execution_context() );
        }
        else
        {
          compute_cell_particles( *grid , *ghost , func, onika::make_flat_tuple(fx,fy,fz,ep) , parallel_execution_context() );
        }
      }
    }

  };
  
  template<class GridT> using ZeroForceEnergy = ZeroForceEnergyNode<GridT,field::_fx,field::_fy,field::_fz, field::_ep, field::_virial, field::_couple >;
  template<class GridT> using ZeroForceEnergyFlat = ZeroForceEnergyNode<GridT,field::_flat_fx,field::_flat_fy,field::_flat_fz, field::_flat_ep, field::_virial, field::_couple>;
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(zero_force_energy)
  {
   OperatorNodeFactory::instance()->register_factory( "zero_force_energy", make_grid_variant_operator< ZeroForceEnergy > );
   OperatorNodeFactory::instance()->register_factory( "zero_force_energy_flat", make_grid_variant_operator< ZeroForceEnergyFlat > );
  }

}

