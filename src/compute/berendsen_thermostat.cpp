#include <memory>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/quantity.h>
#include <exanb/core/physics_constants.h>
#include <exanb/core/unityConverterHelper.h>
#include <onika/memory/allocator.h>
#include <exanb/core/parallel_random.h>
#include <exanb/grid_cell_particles/particle_region.h>
#include <exaStamp/compute/thermodynamic_state.h>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  class BerendsenThermostatNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
    static constexpr bool has_id_field = has_id_field_t::value;

    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );
    ADD_SLOT( double         , tau     , INPUT , 0.1 );
    ADD_SLOT( double         , T       , INPUT , REQUIRED );
    ADD_SLOT( double         , dt      , INPUT , REQUIRED );
    ADD_SLOT( ThermodynamicState      , thermodynamic_state , INPUT );

    ADD_SLOT( ParticleRegions   , particle_regions , INPUT , OPTIONAL );
    ADD_SLOT( ParticleRegionCSG , region           , INPUT , OPTIONAL );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      static const double k = UnityConverterHelper::convert(legacy_constant::boltzmann, "J/K");

      if( grid->number_of_cells() == 0 ) return;

      GridT& grid              = *(this->grid);
      const double tau         = *(this->tau);
      const double T           = *(this->T);
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      double dt                = *(this->dt);
      ParticleSpecies& species = *(this->species);
      long natoms = sim_info.particle_count();
      static constexpr double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann;      
      double T_cur = sim_info.temperature_scal() / natoms * conv_temperature;

      double lambda = sqrt( 1. + dt * ( T/T_cur - 1 ) / tau );
      ldbg << "berendsen: tau="<<tau<<", T="<<T<<", dt="<<dt<<std::endl;

      size_t nSpecies = species.size();
      double masses[ nSpecies ];
      for(size_t i=0;i<nSpecies;i++)
      {
        masses[i] = species[i].m_mass;
      }
      if( !has_type_field && nSpecies>1 )
      {
        fatal_error() << "There are "<<nSpecies<<" atom species but grid has no type field"<<std::endl;
      }

      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();      

      ParticleRegionCSGShallowCopy prcsg;
      if( region.has_value() && !particle_regions.has_value() )
      {
        fatal_error() << "region is defined, but particle_regions has no value" << std::endl;
      }        
      if( region.has_value() && region->m_nb_operands==0 )
      {
        region->build_from_expression_string( particle_regions->data() , particle_regions->size() );
      }
      if( region.has_value() ) prcsg = *region;

#     pragma omp parallel
      {

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );

          const auto* __restrict__ rx = cells[i][field::rx];
          const auto* __restrict__ ry = cells[i][field::ry];
          const auto* __restrict__ rz = cells[i][field::rz];
          [[maybe_unused]] const auto* __restrict__ ids = cells[i].field_pointer_or_null(field::id);

          auto* __restrict__ fx = cells[i][field::fx]; ONIKA_ASSUME_ALIGNED(fx);
          auto* __restrict__ fy = cells[i][field::fy]; ONIKA_ASSUME_ALIGNED(fy);
          auto* __restrict__ fz = cells[i][field::fz]; ONIKA_ASSUME_ALIGNED(fz);

          const auto* __restrict__ vx = cells[i][field::vx];
          const auto* __restrict__ vy = cells[i][field::vy];
          const auto* __restrict__ vz = cells[i][field::vz];

          const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);
          const unsigned int n = cells[i].size();

          for(unsigned int j=0;j<n;j++)
          {
            double mass = masses[0];
            if constexpr ( has_type_field ) { mass = masses[ atom_type[j] ]; }
            uint64_t p_id = 0;
            if constexpr (has_id_field) { p_id = ids[j]; }
            if( prcsg.contains( Vec3d{rx[j],ry[j],rz[j]} , p_id ) )
            {
              fx[j] *= lambda ;
              fy[j] *= lambda ;
              fz[j] *= lambda ;              
            }
          }

        }
        GRID_OMP_FOR_END
	    }

    }

    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsScalar() )
      {
        tmp["T"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize( tmp );
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Apply a berendsen thermostat on particles.
if a single value is given as the node description, it is understood as the T parameter (temperature)

Uses formulation as in LAMMPS documentation :
=============================================
Apply a Berendsen thermostat as described in (Schneider) to a group of atoms which models an interaction with a background implicit solvent. Used with fix nve, this command performs Brownian dynamics (BD), since the total force on each atom will have the form:
F = Fc + Ff + Fr
Ff = - (m / damp) v
Fr is proportional to sqrt(Kb T m / (dt damp))
Fc is the conservative force computed via the usual inter-particle interactions (pair_style, bond_style, etc).
The Ff and Fr terms are added by this fix on a per-particle basis. See the pair_style dpd/tstat command for a thermostatting option that adds similar terms on a pairwise basis to pairs of interacting particles.
Ff is a frictional drag or viscous damping term proportional to the particle’s velocity. The proportionality constant for each atom is computed as m/damp, where m is the mass of the particle and damp is the damping factor specified by the user.
Fr is a force due to solvent atoms at a temperature T randomly bumping into the particle. As derived from the fluctuation/dissipation theorem, its magnitude as shown above is proportional to sqrt(Kb T m / dt damp), where Kb is the Boltzmann constant, T is the desired temperature, m is the mass of the particle, dt is the timestep size, and damp is the damping factor. Random numbers are used to randomize the direction and magnitude of this force as described in (Dunweg), where a uniform random number is used (instead of a Gaussian random number) for speed.

exemple 1:
==========
berendsen_thermostat: 500 K

exemple 2:
==========
berendsen_thermostat:
  T: 800 K
  tau: 0.1 ps

Note: do not process particles in ghost layers.
)EOF";
    }

  };

  template<class GridT> using BerendsenThermostatNodeTmpl = BerendsenThermostatNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory(
    "berendsen_thermostat",
    make_grid_variant_operator< BerendsenThermostatNodeTmpl >
    );
  }

}
