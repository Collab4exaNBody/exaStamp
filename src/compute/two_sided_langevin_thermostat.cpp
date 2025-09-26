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
#include <exanb/core/domain.h>

#include <algorithm>

namespace exaStamp
{

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_vx, field::_vy, field::_vz >
    >
  class TwoSidedLangevinThermostatNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;
    using PlanesVector = std::vector<Plane3d>;

    ADD_SLOT( GridT          , grid      , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species   , INPUT , REQUIRED );
    ADD_SLOT( std::string    , specy     , INPUT );
    ADD_SLOT( double         , gamma_a   , INPUT , 0.1 );
    ADD_SLOT( double         , T_a       , INPUT , REQUIRED );
    ADD_SLOT( double         , gamma_b   , INPUT , 0.1 );
    ADD_SLOT( double         , T_b       , INPUT , REQUIRED );
    ADD_SLOT( double         , dt        , INPUT , REQUIRED );

    ADD_SLOT( PlanesVector   , planes    , INPUT , PlanesVector( { Plane3d{ Vec3d{1.0,0.0,0.0} , 0.0 } } ) );

    ADD_SLOT( double         , thickness , INPUT , 1.0 );
    ADD_SLOT( Domain         , domain    , INPUT , REQUIRED );
    ADD_SLOT( bool             , deterministic_noise , INPUT , false );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      const double k = onika::physics::make_quantity( onika::physics::boltzmann, "J/K" ).convert();

      ldbg << "langevin: gamma_a="<<*gamma_a<<", Ta="<<*T_a<<", dt="<<*dt
           << ", gamma_b="<<*gamma_b<<", Tb="<<*T_b<<", planes:";
      for(const auto& p:*planes) ldbg<<" (N="<<p.N <<" D="<<p.D<<")";
      ldbg<<std::endl;

      double masses[ species->size() ];

      size_t nSpecies = species->size();
      unsigned int specy_index = 0;
      for(size_t i=0;i<nSpecies;i++)
      {
        masses[i] = species->at(i).m_mass;
        if( specy.has_value() && species->at(i).m_name == *specy )
        {
          specy_index = i;
        }
      }
      if(specy_index>=nSpecies) { lerr<<"No atom species found\n"; }
      //specy_index=specy_index;

      auto cells = grid->cells();
      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();
      IJK gstart { gl, gl, gl };
      IJK gend = dims - IJK{ gl, gl, gl };
      IJK gdims = gend - gstart;
      const auto dom_dims = domain->grid_dimension();
      const auto dom_start = grid->offset();
      
      const Mat3d xform = domain->xform();
//      const Plane3d P = *plane;
//      const double D = - (*offset);
      const double Ta = *T_a;
      const double Tb = *T_b;
      const double Ga = *gamma_a;
      const double Gb = *gamma_b;
      const double delta_t = *dt;

      /*
      size_t count_a = 0;
      size_t count_b = 0;
      Vec3d Ke_a = { 0 , 0 , 0 };
      Vec3d Ke_b = { 0 , 0 , 0 };
      Vec3d mom_a = { 0 , 0 , 0 };
      Vec3d mom_b = { 0 , 0 , 0 };
      double mass_a = 0.0;
      double mass_b = 0.0;
      */

      const int nthreads = *deterministic_noise ? 1 : omp_get_max_threads();
#     pragma omp parallel num_threads(nthreads)
      {
        std::mt19937_64 det_re;
        std::mt19937_64 & re = *deterministic_noise ? det_re : onika::parallel::random_engine() ;
        std::normal_distribution<double> f_rand(0.0,1.0);
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          const auto i = grid_ijk_to_index( dims , loc + gstart );
          const auto domain_cell_idx = grid_ijk_to_index( dom_dims , loc + gstart + dom_start );
          
          auto* __restrict__ fx = cells[i][field::fx]; ONIKA_ASSUME_ALIGNED(fx);
          auto* __restrict__ fy = cells[i][field::fy]; ONIKA_ASSUME_ALIGNED(fy);
          auto* __restrict__ fz = cells[i][field::fz]; ONIKA_ASSUME_ALIGNED(fz);

          const auto* __restrict__ vx = cells[i][field::vx];
          const auto* __restrict__ vy = cells[i][field::vy];
          const auto* __restrict__ vz = cells[i][field::vz];

          const auto* __restrict__ rx = cells[i][field::rx];
          const auto* __restrict__ ry = cells[i][field::ry];
          const auto* __restrict__ rz = cells[i][field::rz];

          const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type); if(atom_type==nullptr){};
          const unsigned int n = cells[i].size();

          det_re.seed( domain_cell_idx * 1023 );
          
          for(unsigned int j=0;j<n;j++)
          {
            double mass = 0.0;
            if constexpr ( has_type_field ) { mass = masses[ atom_type[j] ]; }
            else { mass = masses[specy_index]; }
            
            const Vec3d r = xform * Vec3d{ rx[j], ry[j], rz[j] };
            const Vec3d v = { vx[j], vy[j], vz[j] };
     
            double plane_dist = std::numeric_limits<double>::max();
            for(const auto& p:*planes)
            {
              plane_dist = std::min( plane_dist , dot( r , p.N ) + p.D );
            }
            
            double t = std::clamp( plane_dist / (*thickness) , -1.0 , 1.0 ) *0.5 + 0.5;
            
            double gamma = (1.0-t) * Ga + t * Gb; 
            double T = (1.0-t) * Ta + t * Tb;

/*
            if( plane_dist < -(*thickness) )
            {
             ++ count_a;
             mom_a += v * mass;
             Ke_a += v * v * mass;
             mass_a += mass;
            }
            else if( plane_dist > (*thickness) )
            {
             ++ count_b;
             mom_b += v * mass;
             Ke_b += v * v * mass;
             mass_b += mass;
            }
*/
            
            fx[j] +=  /* Ff */ - ( mass / gamma ) * v.x  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( delta_t * gamma ) ) ;
            fy[j] +=  /* Ff */ - ( mass / gamma ) * v.y  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( delta_t * gamma ) ) ;
            fz[j] +=  /* Ff */ - ( mass / gamma ) * v.z  +  /* Fr */ f_rand(re) * std::sqrt( 2.0 * k * T * mass / ( delta_t * gamma ) ) ;
          }

        }
        GRID_OMP_FOR_END
      }
      
/*    
      Vec3d temp_a = 2. * ( Ke_a*0.5 - 0.5 * mom_a * mom_a / mass_a );
      double temp_a_scal = ( conv_temperature * ( temp_a.x + temp_a.y + temp_a.z ) / 3. ) / count_a;

      Vec3d temp_b = 2. * ( Ke_b*0.5 - 0.5 * mom_b * mom_b / mass_b );
      double temp_b_scal = ( conv_temperature * ( temp_b.x + temp_b.y + temp_b.z ) / 3. ) / count_b;

      ldbg<<"count_a="<<count_a<<", count_b="<<count_b<<", temp_a="<<temp_a_scal<<", temp_b="<<temp_b_scal<<std::endl;           
*/      
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Apply a langevin thermostat on particles, with two different set of parameters.
a plane separation tells wich parameter set to use on which particles.

one set of parameters (gamma_a and T_a) is used on particles lying on the negative side of the plane,
the other set (gamma_b and T_a) is used on the positive side of the plane.

the plane is described by a normal and a distance to the origin (parameter offset), such that plane's equation is :
normal . Ri + offset = 0
where Ri is the coordinate of the Ith particle.

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

exemple :
=========
two_sided_langevin_thermostat:
  T_a: 800 K
  gamma_a: 0.25
  T_b: 500 K
  gamma_b: 0.25
  normal: [ 1.0 , 0.0 , 0.0 ]
  offset: 0.5 ang

Note: do not process particles in ghost layers.
)EOF";
    }

  };

  template<class GridT> using TwoSidedLangevinThermostatNodeTmpl = TwoSidedLangevinThermostatNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(two_sided_langevin_thermostat)
  {
   OperatorNodeFactory::instance()->register_factory(
    "two_sided_langevin_thermostat",
    make_grid_variant_operator< TwoSidedLangevinThermostatNodeTmpl >
    );
  }

}
