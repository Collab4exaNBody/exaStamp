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

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  class FixedAtomsNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );
    ADD_SLOT( std::string    , specy   , INPUT );
    ADD_SLOT( Vec3d          , normal , INPUT , Vec3d{1.0,0.0,0.0} );
    ADD_SLOT( double         , min_limit       , INPUT , REQUIRED );
    ADD_SLOT( double         , max_limit       , INPUT , REQUIRED );
    ADD_SLOT( Domain         , domain    , INPUT , REQUIRED );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      //static const double k = UnityConverterHelper::convert(onika::physics::boltzmann, "J/K");

      GridT& grid              = *(this->grid);
      ParticleSpecies& species = *(this->species);
      Vec3d N = *(this->normal);
      double min_limit = *(this->min_limit);
      double max_limit = *(this->max_limit);
      const Mat3d xform = domain->xform();

      size_t nSpecies = species.size();
      m_masses.resize( nSpecies );
      /*
      unsigned int specy_index = 0;
      for(size_t i=0;i<nSpecies;i++)
      {
        m_masses[i] = species[i].m_mass;
        if( specy.has_value() && species[i].m_name == *specy )
        {
          specy_index = i;
        }
      }
      */
      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();      

#     pragma omp parallel
      {

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );

          auto* __restrict__ fx = cells[i][field::fx]; ONIKA_ASSUME_ALIGNED(fx);
          auto* __restrict__ fy = cells[i][field::fy]; ONIKA_ASSUME_ALIGNED(fy);
          auto* __restrict__ fz = cells[i][field::fz]; ONIKA_ASSUME_ALIGNED(fz);

          const auto* __restrict__ rx = cells[i][field::rx];
          const auto* __restrict__ ry = cells[i][field::ry];
          const auto* __restrict__ rz = cells[i][field::rz];

          //const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);
          const unsigned int n = cells[i].size();

	  for(unsigned int j=0;j<n;j++)
            {
	      
	      const Vec3d r = xform * Vec3d{ rx[j], ry[j], rz[j] };
	      
	      double dist_to_min = dot( r , N ) - min_limit;
	      double dist_to_max = dot( r , N ) - max_limit;

	      if ( (dist_to_min >= 0.) && (dist_to_max <= 0.) ) {
		fx[j] = 0.;
		fy[j] = 0.;
		fz[j] = 0.;
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

  private:
    std::vector<double> m_masses;
  };

  template<class GridT> using FixedAtomsNodeTmpl = FixedAtomsNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(fixed_atoms)
  {
   OperatorNodeFactory::instance()->register_factory(
    "fixed_atoms",
    make_grid_variant_operator< FixedAtomsNodeTmpl >
    );
  }

}
