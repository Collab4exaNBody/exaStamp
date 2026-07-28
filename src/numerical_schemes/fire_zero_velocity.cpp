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

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  // Zeroes vx,vy,vz for every atom: FIRE's initial v(t)=0, and the "set velocity to zero" step
  // of the uphill-motion correction (Bitzek et al., F3/F4; Guenole et al. Algorithm 2 line 29).
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz >
    >
  struct FIREZeroVelocity : public OperatorNode
  {
    ADD_SLOT( GridT , grid , INPUT_OUTPUT );

    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();
      auto cells = grid.cells();

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          auto* __restrict__ vx = cells[i][field::vx];
          auto* __restrict__ vy = cells[i][field::vy];
          auto* __restrict__ vz = cells[i][field::vz];
#         pragma omp simd
          for(size_t j=0;j<n;j++)
          {
            vx[j] = 0.0;
            vy[j] = 0.0;
            vz[j] = 0.0;
          }
        }
        GRID_OMP_FOR_END
      }
    }
  };

  template<class GridT> using FIREZeroVelocityTmpl = FIREZeroVelocity<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(fire_zero_velocity)
  {
   OperatorNodeFactory::instance()->register_factory( "fire_zero_velocity", make_grid_variant_operator< FIREZeroVelocityTmpl > );
  }

}
