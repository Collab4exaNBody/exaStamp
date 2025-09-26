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
//#include "exanb/memory.h"
#include <onika/physics/constants.h>
#include <exanb/core/domain.h>

namespace exaStamp
{

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_rxf, field::_ryf, field::_rzf >
    >
  class PositionsFilteringNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    ADD_SLOT( GridT              , grid    , INPUT_OUTPUT);
    ADD_SLOT( double             , alpha   , INPUT , 0.01 );
    ADD_SLOT( double             , dt      , INPUT , REQUIRED );
    ADD_SLOT( uint64_t           , timestep            , INPUT );
    ADD_SLOT( double             , time_start, INPUT);
    ADD_SLOT( Domain , domain    , INPUT , REQUIRED );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {

      GridT& grid              = *(this->grid);
      const double alpha       = *(this->alpha);
      double dt                = *(this->dt);
      double starttime         = *(this->time_start);
      double curtime           = dt * (*timestep);
      
      ldbg << "atomic positions filtering: alpha="<<alpha<<std::endl;

      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();      
      double pi = M_PI;

      double xmin = domain->bounds().bmin.x;
      double xmax = domain->bounds().bmax.x;
      double ymin = domain->bounds().bmin.y;
      double ymax = domain->bounds().bmax.y;
      double zmin = domain->bounds().bmin.z;
      double zmax = domain->bounds().bmax.z;
#     pragma omp parallel
      {

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );

          auto* __restrict__ rxf = cells[i][field::rxf];
          auto* __restrict__ ryf = cells[i][field::ryf];
          auto* __restrict__ rzf = cells[i][field::rzf];

          const auto* __restrict__ rx = cells[i][field::rx];
          const auto* __restrict__ ry = cells[i][field::ry];
          const auto* __restrict__ rz = cells[i][field::rz];

          const unsigned int n = cells[i].size();
          double dxsi, dzeta, theta,thetaf;
	  
	  if (curtime <= starttime) {
            for(unsigned int j=0;j<n;j++)
	      {
		rxf[j] = rx[j];
		ryf[j] = ry[j];
		rzf[j] = rz[j];
	      }
	  } else {
            for(unsigned int j=0;j<n;j++)
	      {
		// operation for x component
		theta = (rx[j]-xmin)/(xmax-xmin) * 2.0 * pi;
		thetaf = (rxf[j]-xmin)/(xmax-xmin) * 2.0 * pi;		
		dxsi = alpha * cos(theta) + (1. - alpha) * cos(thetaf);
		dzeta = alpha * sin(theta) + (1. - alpha) * sin(thetaf);		
		theta = atan2(-dzeta,-dxsi)+pi;
		rxf[j] = theta/(2.0 * pi)*(xmax-xmin)+xmin;

		// operation for y component
		theta = (ry[j]-ymin)/(ymax-ymin) * 2.0 * pi;
		thetaf = (ryf[j]-ymin)/(ymax-ymin) * 2.0 * pi;		
		dxsi = alpha * cos(theta) + (1. - alpha) * cos(thetaf);
		dzeta = alpha * sin(theta) + (1. - alpha) * sin(thetaf);		
		theta = atan2(-dzeta,-dxsi)+pi;
		ryf[j] = theta/(2.0 * pi)*(ymax-ymin)+ymin;

		// operation for z component
		theta = (rz[j]-zmin)/(zmax-zmin) * 2.0 * pi;
		thetaf = (rzf[j]-zmin)/(zmax-zmin) * 2.0 * pi;		
		dxsi = alpha * cos(theta) + (1. - alpha) * cos(thetaf);
		dzeta = alpha * sin(theta) + (1. - alpha) * sin(thetaf);		
		theta = atan2(-dzeta,-dxsi)+pi;
		rzf[j] = theta/(2.0 * pi)*(zmax-zmin)+zmin;
		
	      }
	  }
	  
        }
        GRID_OMP_FOR_END
      }
    }

  };

  template<class GridT> using PositionsFilteringNodeTmpl = PositionsFilteringNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(positions_filtering)
  {
   OperatorNodeFactory::instance()->register_factory(
    "positions_filtering",
    make_grid_variant_operator< PositionsFilteringNodeTmpl >
    );
  }

}
