#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/quantity.h>
#include <exanb/core/physics_constants.h>
#include <exanb/core/unityConverterHelper.h>
#include <onika/parallel/random.h>
#include <onika/physics/units.h>
#include <onika/math/basic_types_yaml.h>
#include <exaStamp/particle_species/particle_specie_yaml.h>

#include <yaml-cpp/yaml.h>
#include <string>

namespace exaStamp
{
  struct MaterialLangevin
  {
    double gamma = 0.1;
    double T = 0.0;
    double mass = 0.0;
    std::string type;
  };

}

namespace YAML
{
  template<> struct convert< exaStamp::MaterialLangevin >
  {
    static inline bool decode(const Node& node, exaStamp::MaterialLangevin& v)
    {
      using exanb::Quantity;
      
      if( ! node.IsMap() ) return false;
      if( ! node["gamma"] ) return false;
      if( ! node["T"] ) return false;
      if( ! node["type"] ) return false;

      v.gamma = node["gamma"].as<Quantity>().convert();
      v.T = node["T"].as<Quantity>().convert();
      v.mass = 0.0;
      v.type = node["type"].as<std::string>();      
      return true;
    }
  };
}


namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_vx, field::_vy, field::_vz, field::_type >
    >
  class MaterialLangevinThermostat : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );
    ADD_SLOT( double         , dt      , INPUT , REQUIRED );
    ADD_SLOT( std::vector<MaterialLangevin> , matlangevin , INPUT , REQUIRED );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      static const double k = UnityConverterHelper::convert(legacy_constant::boltzmann, "J/K");

      ldbg << "per material langevin: dt="<<*dt<<std::endl;

      std::map< std::string ,  MaterialLangevin > map_langevin;
      for(const auto& l : *matlangevin)
      {
        map_langevin[ l.type ] = MaterialLangevin{ l.gamma, l.T };
      }

      const size_t nSpecies = species->size();
      std::vector<MaterialLangevin> langevin_params( nSpecies );
      for(size_t i=0;i<nSpecies;i++)
      {
        auto it = map_langevin.find( species->at(i).name() );
        if( it != map_langevin.end() )
        {
          langevin_params[i] = it->second;
          langevin_params[i].mass = species->at(i).m_mass;
          ldbg<<"\t"<<species->at(i).name()<<" : gamma="<<langevin_params[i].gamma<<", T="<<langevin_params[i].T<<", mass="<<langevin_params[i].mass<<std::endl;
        }
      }

      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();
      const double dt = *(this->dt);    

#     pragma omp parallel
      {
        auto& re = rand::random_engine();
        std::normal_distribution<double> f_rand(0.0,1.0);

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );

          auto* __restrict__ fx = cells[i][field::fx]; ONIKA_ASSUME_ALIGNED(fx);
          auto* __restrict__ fy = cells[i][field::fy]; ONIKA_ASSUME_ALIGNED(fy);
          auto* __restrict__ fz = cells[i][field::fz]; ONIKA_ASSUME_ALIGNED(fz);

          const auto* __restrict__ vx = cells[i][field::vx]; ONIKA_ASSUME_ALIGNED(vx);
          const auto* __restrict__ vy = cells[i][field::vy]; ONIKA_ASSUME_ALIGNED(vy);
          const auto* __restrict__ vz = cells[i][field::vz]; ONIKA_ASSUME_ALIGNED(vz);

          const auto* __restrict__ atom_type = cells[i][field::type]; ONIKA_ASSUME_ALIGNED(atom_type);
          const unsigned int n = cells[i].size();

          for(unsigned int j=0;j<n;j++)
          {
            const unsigned int type = atom_type[j];
            const double mass = langevin_params[ type ].mass;
            const double gamma = langevin_params[ type ].gamma;
            const double T = langevin_params[ type ].T;
            
            fx[j] +=  /* Ff */ - ( mass / gamma ) * vx[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
            fy[j] +=  /* Ff */ - ( mass / gamma ) * vy[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
            fz[j] +=  /* Ff */ - ( mass / gamma ) * vz[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
          }

        }
        GRID_OMP_FOR_END
	    }
    }

    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsSequence() )
      {
        tmp["matlangevin"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize( tmp );
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Apply a langevin thermostat on particles.
if a single value is given as the node description, it is understood as the T parameter (temperature)

Uses formulation as in LAMMPS documentation :
=============================================
Apply a Langevin thermostat as described in (Schneider) to a group of atoms which models an interaction with a background implicit solvent. Used with fix nve, this command performs Brownian dynamics (BD), since the total force on each atom will have the form:
F = Fc + Ff + Fr
Ff = - (m / damp) v
Fr is proportional to sqrt(Kb T m / (dt damp))
Fc is the conservative force computed via the usual inter-particle interactions (pair_style, bond_style, etc).
The Ff and Fr terms are added by this fix on a per-particle basis. See the pair_style dpd/tstat command for a thermostatting option that adds similar terms on a pairwise basis to pairs of interacting particles.
Ff is a frictional drag or viscous damping term proportional to the particle’s velocity. The proportionality constant for each atom is computed as m/damp, where m is the mass of the particle and damp is the damping factor specified by the user.
Fr is a force due to solvent atoms at a temperature T randomly bumping into the particle. As derived from the fluctuation/dissipation theorem, its magnitude as shown above is proportional to sqrt(Kb T m / dt damp), where Kb is the Boltzmann constant, T is the desired temperature, m is the mass of the particle, dt is the timestep size, and damp is the damping factor. Random numbers are used to randomize the direction and magnitude of this force as described in (Dunweg), where a uniform random number is used (instead of a Gaussian random number) for speed.

exemple 1:
==========
langevin_thermostat: 500 K

exemple 2:
==========
langevin_thermostat:
  T: 800 K
  gamma: 0.25

Note: do not process particles in ghost layers.
)EOF";
    }

  };

  template<class GridT> using MaterialLangevinThermostatTmpl = MaterialLangevinThermostat<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(material_thermostat)
  {
   OperatorNodeFactory::instance()->register_factory(
    "material_langevin_thermostat",
    make_grid_variant_operator< MaterialLangevinThermostatTmpl >
    );
  }

}
