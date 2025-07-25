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

#include <memory>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <onika/parallel/random.h>
#include <exanb/grid_cell_particles/particle_region.h>
#include <exaStamp/unit_system.h>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  class LangevinThermostatNode : public OperatorNode
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
    ADD_SLOT( ParticleRegions   , particle_regions , INPUT , OPTIONAL );
    ADD_SLOT( ParticleRegionCSG , region           , INPUT , OPTIONAL );
    ADD_SLOT( double         , gamma   , INPUT , 0.1 );
    ADD_SLOT( size_t         , seed    , INPUT , OPTIONAL );
    ADD_SLOT( double         , T       , INPUT , OPTIONAL );
    ADD_SLOT( double         , Tstart  , INPUT , OPTIONAL );
    ADD_SLOT( double         , Tstop   , INPUT , OPTIONAL );
    ADD_SLOT( TimeVec        , tserie  , INPUT , OPTIONAL );
    ADD_SLOT( TempVec        , Tserie  , INPUT , OPTIONAL );
    ADD_SLOT( bool             , deterministic_noise , INPUT , false );
    
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
      static constexpr double k = EXASTAMP_CONST_QUANTITY( onika::physics::boltzmann * J / K );

      if( grid->number_of_cells() == 0 ) return;

      GridT& grid              = *(this->grid);
      const double gamma       = *(this->gamma);
      double dt                = *(this->dt);
      ParticleSpecies& species = *(this->species);

      // Checking definition of target temperature
      bool constant_T = T.has_value();
      bool linear_T = Tstart.has_value() && Tstop.has_value();
      bool interpolated_T = tserie.has_value() && Tserie.has_value();
      if ( (constant_T && linear_T && interpolated_T) || (constant_T && linear_T) || (constant_T && interpolated_T) || (linear_T && interpolated_T)) {
        lerr << "Multiple definition of target temperature are provided" << std::endl;
        lout << "You must define the target temperature using one of the three following solutions:" << std::endl;
        lout << "langevin_thermostat:" << std::endl;
        lout << "  T: 300 K" << std::endl;
        lout << "  gamma: 0.1 ps" << std::endl;
        lout << "langevin_thermostat:" << std::endl;
        lout << "  Tstart: 300 K" << std::endl;
        lout << "  Tstop: 300 K" << std::endl;        
        lout << "  gamma: 0.1 ps" << std::endl;
        lout << "langevin_thermostat:" << std::endl;
        lout << "  tserie: [0, 1]" << std::endl;
        lout << "  Tserie: [300, 300]" << std::endl;        
        lout << "  gamma: 0.1 ps" << std::endl;
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
      ldbg << "langevin: gamma="<<gamma<<", T="<<Ttarget<<", dt="<<dt<<std::endl;
      
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
      
      const int nthreads = *deterministic_noise ? 1 : omp_get_max_threads();

#     pragma omp parallel num_threads(nthreads)
      {
        std::mt19937_64 det_re;
        std::mt19937_64 & re = *deterministic_noise ? det_re : onika::parallel::random_engine() ;
        std::normal_distribution<double> f_rand(0.,1.);

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
            det_re.seed( ( ( (i*1023) ^ j ) * 1023 ) ^ (*timestep) );
            double mass = masses[0];
            if constexpr ( has_type_field ) { mass = masses[ atom_type[j] ]; }
            uint64_t p_id = 0;
            if constexpr (has_id_field) { p_id = ids[j]; }
            if( prcsg.contains( Vec3d{rx[j],ry[j],rz[j]} , p_id ) )
            {
              fx[j] +=  /* Ff */ - ( mass / gamma ) * vx[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * Ttarget * mass / ( dt * gamma ) ) ;
              fy[j] +=  /* Ff */ - ( mass / gamma ) * vy[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * Ttarget * mass / ( dt * gamma ) ) ;
              fz[j] +=  /* Ff */ - ( mass / gamma ) * vz[j]  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * Ttarget * mass / ( dt * gamma ) ) ;
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

  template<class GridT> using LangevinThermostatNodeTmpl = LangevinThermostatNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(langevin_thermostat)
  {
   OperatorNodeFactory::instance()->register_factory(
    "langevin_thermostat",
    make_grid_variant_operator< LangevinThermostatNodeTmpl >
    );
  }

}
