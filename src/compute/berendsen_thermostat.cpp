#include <memory>

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
#include <onika/memory/allocator.h>
#include <onika/parallel/random.h>
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

    // compile time constant indicating if grid has id field
    using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
    static constexpr bool has_id_field = has_id_field_t::value;

    using TimeVec = std::vector<double>;
    using TempVec = std::vector<double>;
    
    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );
    ADD_SLOT( double         , dt      , INPUT , REQUIRED );
    ADD_SLOT( long           , timestep     , INPUT , REQUIRED );
    ADD_SLOT( long           , simulation_end_iteration , INPUT , REQUIRED );
    ADD_SLOT( long           , simulation_start_iteration , INPUT , 0 );
    ADD_SLOT( ThermodynamicState      , thermodynamic_state , INPUT );
    ADD_SLOT( ParticleRegions   , particle_regions , INPUT , OPTIONAL );
    ADD_SLOT( ParticleRegionCSG , region           , INPUT , OPTIONAL );
    ADD_SLOT( double         , tau     , INPUT , 0.1 );
    ADD_SLOT( double         , T       , INPUT , OPTIONAL );
    ADD_SLOT( double         , Tstart  , INPUT , OPTIONAL );
    ADD_SLOT( double         , Tstop   , INPUT , OPTIONAL );
    ADD_SLOT( TimeVec        , tserie  , INPUT , OPTIONAL );
    ADD_SLOT( TempVec        , Tserie  , INPUT , OPTIONAL );    

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

      if( grid->number_of_cells() == 0 ) return;

      GridT& grid              = *(this->grid);
      const double tau         = *(this->tau);
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      double dt                = *(this->dt);
      ParticleSpecies& species = *(this->species);

      // Getting current temperature
      static constexpr double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann;
      double Tcurrent = sim_info.temperature_scal() / sim_info.particle_count() * conv_temperature;

      // Checking definition of target temperature
      bool constant_T = T.has_value();
      bool linear_T = Tstart.has_value() && Tstop.has_value();
      bool interpolated_T = tserie.has_value() && Tserie.has_value();
      if ( (constant_T && linear_T && interpolated_T) || (constant_T && linear_T) || (constant_T && interpolated_T) || (linear_T && interpolated_T)) {
        lerr << "Multiple definition of target temperature are provided" << std::endl;
        lout << "You must define the target temperature using one of the three following solutions:" << std::endl;
        lout << "berendsen_thermostat:" << std::endl;
        lout << "  T: 300 K" << std::endl;
        lout << "  tau: 0.1 ps" << std::endl;
        lout << "berendsen_thermostat:" << std::endl;
        lout << "  Tstart: 300 K" << std::endl;
        lout << "  Tstop: 300 K" << std::endl;        
        lout << "  tau: 0.1 ps" << std::endl;
        lout << "berendsen_thermostat:" << std::endl;
        lout << "  tserie: [0, 1]" << std::endl;
        lout << "  Tserie: [300, 300]" << std::endl;        
        lout << "  tau: 0.1 ps" << std::endl;
        std::abort();
      }

      // Getting target temperature      
      double Ttarget = 0.;
      if (constant_T) {
        Ttarget = *T;
      } else if (linear_T) {
        double delta = (*timestep)-(*simulation_start_iteration);
        double frac = delta/((*simulation_end_iteration) - (*simulation_start_iteration));
        Ttarget = *Tstart + frac * ( (*Tstop)-(*Tstart) );
      } else if (interpolated_T) {
        double tcurrent = dt * (*timestep);
        const TempVec& tempvec = *Tserie;        
        Ttarget = interpolate( *tserie , tcurrent , [&tempvec](size_t i)->double { return tempvec[i]; } );
      }
      const double lambda = sqrt( 1. + dt/tau*(Ttarget/Tcurrent - 1.0));

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

          auto* __restrict__ vx = cells[i][field::vx];
          auto* __restrict__ vy = cells[i][field::vy];
          auto* __restrict__ vz = cells[i][field::vz];

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
              vx[j] *= lambda ;
              vy[j] *= lambda ;
              vz[j] *= lambda ;              
            }
          }

        }
        GRID_OMP_FOR_END
	    }

    }

  };

  template<class GridT> using BerendsenThermostatNodeTmpl = BerendsenThermostatNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(berendsen_thermostat)
  {
   OperatorNodeFactory::instance()->register_factory(
    "berendsen_thermostat",
    make_grid_variant_operator< BerendsenThermostatNodeTmpl >
    );
  }

}
