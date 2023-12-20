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

namespace exaStamp
{

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  class PiecewiseLangevinThermostatNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    using TimeVec = std::vector<double>;
    using TempVec = std::vector<double>;
    
    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );
    ADD_SLOT( std::string    , specy   , INPUT );
    ADD_SLOT( double         , gamma   , INPUT , 0.1 );
    ADD_SLOT( TimeVec        , time_serie  , INPUT , REQUIRED );
    ADD_SLOT( TempVec        , temperature_serie  , INPUT , REQUIRED );    
    ADD_SLOT( double         , dt      , INPUT , REQUIRED );
    ADD_SLOT( int64_t        , timestep     , INPUT , REQUIRED );
    ADD_SLOT( double         , time_nve      , INPUT );

    template<typename YFunc>
    static inline double interpolate( const std::vector<double>& X , double ix, YFunc yfunc )
    {
      assert( std::is_sorted( X.begin() , X.end() ) );
      size_t N = X.size();
      
      if( N == 0 )
      {
        return 0.0;
      }
      if( N == 1 )
      {
        return yfunc(0);
      }
      if( ix < X[0] )
      {
        return yfunc(0);
      }
      
      size_t k = 0;
      while( ix > X[k+1] && k<(N-2) ) { ++k; }
      
      if( ix > X[k+1] )
      {
        return yfunc(k+1);
      }

      double t = (ix - X[k]) / ( X[k+1] - X[k] );
      assert( t>=0 && t<=1.0 );

      return yfunc(k)*(1.-t) + yfunc(k+1)*t;
    }
    
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      if( temperature_serie->size() != time_serie->size() )
      {
        lerr << "number of time values does not match number of temperature." << std::endl;
        std::abort();
      }
      
      static const double k = UnityConverterHelper::convert(legacy_constant::boltzmann, "J/K");

      GridT& grid              = *(this->grid);
      const double gamma       = *(this->gamma);
      double dt                = *(this->dt);
      ParticleSpecies& species = *(this->species);
      double curtime = dt * (*timestep);
      const TempVec& tempvec = *temperature_serie;
      
      const double T = interpolate( *time_serie , curtime , [&tempvec](size_t i)->double { return tempvec[i]; } );

      ldbg << "langevin: gamma="<<gamma<<", T="<<T<<", dt="<<dt<<std::endl;

      size_t nSpecies = species.size();
      m_masses.resize( nSpecies );
      unsigned int specy_index = 0;
      for(size_t i=0;i<nSpecies;i++)
	{
	  m_masses[i] = species[i].m_mass;
	  if( specy.has_value() && species[i].m_name == *specy )
	    {
	      specy_index = i;
	    }
	}

      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();      

      if (time_nve.has_value()) {
	double nvetime = (*time_nve);

	if(curtime < nvetime) {
	
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

	      const auto* __restrict__ vx = cells[i][field::vx];
	      const auto* __restrict__ vy = cells[i][field::vy];
	      const auto* __restrict__ vz = cells[i][field::vz];

	      const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);
	      const unsigned int n = cells[i].size();

	      if( atom_type == nullptr )
		{
		  const double mass = m_masses[specy_index];
		  for(unsigned int j=0;j<n;j++)
		    {
		      fx[j] +=  /* Ff */ - ( mass / gamma ) * vx[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fy[j] +=  /* Ff */ - ( mass / gamma ) * vy[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fz[j] +=  /* Ff */ - ( mass / gamma ) * vz[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		    }
		}
	      else
		{
		  for(unsigned int j=0;j<n;j++)
		    {
		      const double mass = m_masses[ atom_type[j] ];
		      fx[j] +=  /* Ff */ - ( mass / gamma ) * vx[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fy[j] +=  /* Ff */ - ( mass / gamma ) * vy[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fz[j] +=  /* Ff */ - ( mass / gamma ) * vz[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		    }
		}

	    }
	  GRID_OMP_FOR_END
	    }
	}
      } else {
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

	      const auto* __restrict__ vx = cells[i][field::vx];
	      const auto* __restrict__ vy = cells[i][field::vy];
	      const auto* __restrict__ vz = cells[i][field::vz];

	      const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);
	      const unsigned int n = cells[i].size();

	      if( atom_type == nullptr )
		{
		  const double mass = m_masses[specy_index];
		  for(unsigned int j=0;j<n;j++)
		    {
		      fx[j] +=  /* Ff */ - ( mass / gamma ) * vx[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fy[j] +=  /* Ff */ - ( mass / gamma ) * vy[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fz[j] +=  /* Ff */ - ( mass / gamma ) * vz[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		    }
		}
	      else
		{
		  for(unsigned int j=0;j<n;j++)
		    {
		      const double mass = m_masses[ atom_type[j] ];
		      fx[j] +=  /* Ff */ - ( mass / gamma ) * vx[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fy[j] +=  /* Ff */ - ( mass / gamma ) * vy[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		      fz[j] +=  /* Ff */ - ( mass / gamma ) * vz[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( dt * gamma ) ) ;
		    }
		}

	    }
	  GRID_OMP_FOR_END
	    }
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
Apply a langevin thermostat on particles.
if a single value is given as the node description, it is understood as the T parameter (temperature)

Uses formulation as in LAMMPS documentation :
=============================================
Apply a PiecewiseLangevin thermostat as described in (Schneider) to a group of atoms which models an interaction with a background implicit solvent. Used with fix nve, this command performs Brownian dynamics (BD), since the total force on each atom will have the form:
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

  private:
    std::vector<double> m_masses;
  };

  template<class GridT> using PiecewiseLangevinThermostatNodeTmpl = PiecewiseLangevinThermostatNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory(
    "piecewise_langevin_thermostat",
    make_grid_variant_operator< PiecewiseLangevinThermostatNodeTmpl >
    );
  }

}
