// #pragma xstamp_cuda_enable  // DO NOT REMOVE THIS LINE !!

// #pragma xstamp_grid_variant  // DO NOT REMOVE THIS LINE !!

#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/particle_species/particle_specie_yaml.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/quantity.h>
#include <exanb/compute/compute_cell_particles.h>

namespace exaStamp
{
  using namespace exanb;

  struct ForceToAccelFunctor
  {
    ParticleSpecie* m_species = nullptr;
    double m_defaultMass = 1.0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double& fx, double & fy, double & fz, unsigned int type) const
    {       
      const double m = m_species[type].m_mass;
      fx /= m;
      fy /= m;
      fz /= m;
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double& fx, double & fy, double & fz) const
    {
      fx /= m_defaultMass;
      fy /= m_defaultMass;
      fz /= m_defaultMass;
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::ForceToAccelFunctor>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  template<class GridT, class _FX, class _FY, class _FZ, class _Type >
  class ForceToAcceleration : public OperatorNode
  {
    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      ForceToAccelFunctor func = { species->data() , species->at(0).m_mass };
      auto fx = grid->field_accessor( onika::soatl::FieldId<_FX>{} );
      auto fy = grid->field_accessor( onika::soatl::FieldId<_FY>{} );
      auto fz = grid->field_accessor( onika::soatl::FieldId<_FZ>{} );

      onika::soatl::FieldId<_Type> type_field = {};
      if( grid->has_allocated_field( type_field ) )
      {    
        auto type = grid->field_accessor( type_field );
        auto cp_fields = onika::make_flat_tuple( fx, fy, fz, type );
        compute_cell_particles( *grid , false , func , cp_fields , parallel_execution_context() );
      }
      else
      {
        auto cp_fields = onika::make_flat_tuple( fx, fy, fz );
        compute_cell_particles( *grid , false , func , cp_fields , parallel_execution_context() );
      }
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Converts force field to acceleration field, dividing it by particle mass
)EOF";
    }

  };

  template<class GridT> using ForceToAccelerationTmpl = ForceToAcceleration<GridT, field::_fx, field::_fy, field::_fz, field::_type >;
  template<class GridT> using ForceToAccelerationFlat = ForceToAcceleration<GridT, field::_flat_fx, field::_flat_fy, field::_flat_fz, field::_flat_type >;

  // === register factories ===
  ONIKA_AUTORUN_INIT(force_to_accel)
  {
   OperatorNodeFactory::instance()->register_factory( "force_to_accel", make_grid_variant_operator< ForceToAccelerationTmpl > );
   OperatorNodeFactory::instance()->register_factory( "force_to_accel_flat", make_grid_variant_operator< ForceToAccelerationFlat > );
  }

}
