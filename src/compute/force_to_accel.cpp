// #pragma xstamp_cuda_enable  // DO NOT REMOVE THIS LINE !!

// #pragma xstamp_grid_variant  // DO NOT REMOVE THIS LINE !!

#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/particle_species/particle_specie_yaml.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/fields.h>
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

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz >
    >
  class ForceToAcceleration : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;
    using PointerTuple = onika::soatl::FieldPointerTuple< GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize , field::_ax , field::_ay , field::_az >;
    using f2a_field_set_t = std::conditional_t< has_type_field ,
                                                 FieldSet<field::_fx, field::_fy, field::_fz,field::_type> ,
                                                 FieldSet<field::_fx, field::_fy, field::_fz> >;
    static constexpr f2a_field_set_t f2a_field_set {};

    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      ForceToAccelFunctor func = { species->data() , species->at(0).m_mass };
      compute_cell_particles( *grid , false , func , f2a_field_set , parallel_execution_context() );
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

  template<class GridT> using ForceToAccelerationTmpl = ForceToAcceleration<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "force_to_accel", make_grid_variant_operator< ForceToAccelerationTmpl > );
  }

}
