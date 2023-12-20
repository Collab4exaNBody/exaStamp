#pragma xstamp_cuda_enable

#pragma xstamp_grid_variant

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
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

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_ep >
    >
  class ZeroForceEnergyNode : public OperatorNode
  {
    // compile time constant indicating if grid has type virial
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr bool has_virial_field = has_virial_field_t::value;

    using has_couple_field_t = typename GridT::CellParticles::template HasField < field::_couple > ;
    static constexpr bool has_couple_field = has_couple_field_t::value;

    using zero_field_set_t = std::conditional_t< has_couple_field ,
                                    std::conditional_t< has_virial_field ,
                                                 FieldSet<field::_fx, field::_fy, field::_fz,field::_ep,field::_couple,field::_virial> ,
                                                 FieldSet<field::_fx, field::_fy, field::_fz,field::_ep,field::_couple> > ,
                                    std::conditional_t< has_virial_field ,
                                                 FieldSet<field::_fx, field::_fy, field::_fz,field::_ep,field::_virial> ,
                                                 FieldSet<field::_fx, field::_fy, field::_fz,field::_ep> > > ;
    static constexpr zero_field_set_t zero_field_set{};
  
    ADD_SLOT( GridT , grid  , INPUT_OUTPUT );
    ADD_SLOT( bool  , ghost  , INPUT , false );

  public:

    inline void execute () override final
    {
      compute_cell_particles( *grid , *ghost , ZeroCellParticleFields{} , zero_field_set , parallel_execution_context() );
    }

  };
  
  template<class GridT> using ZeroForceEnergyNodeTmpl = ZeroForceEnergyNode<GridT>;
  
 // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "zero_force_energy", make_grid_variant_operator< ZeroForceEnergyNodeTmpl > );
  }

}

