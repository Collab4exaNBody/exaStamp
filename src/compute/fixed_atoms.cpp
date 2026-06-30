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
TODO
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
