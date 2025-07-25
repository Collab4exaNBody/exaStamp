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
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/log.h>
#include <onika/parallel/random.h>
#include <onika/math/quaternion_operators.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_orient >
    >
  class RandomOrient : public OperatorNode
  {
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );

  public:
    inline void execute () override final
    {
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);

#     pragma omp parallel
      {
        auto& re = onika::parallel::random_engine();
        std::normal_distribution<double> f_rand(-1.0 , 1.0) ;
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, schedule(dynamic) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);
          auto* __restrict__ orient = cells[cell_i][field::orient];
          size_t n = cells[cell_i].size();
          for(size_t j=0;j<n;j++)
          {
            orient[j] = normalize( Quaternion{ f_rand(re) , f_rand(re) , f_rand(re) , f_rand(re) } );
          }
        }
        GRID_OMP_FOR_END
      }
    }
  };

  template<class GridT> using RandomOrientTmpl = RandomOrient<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(random_orient)
  {
    OperatorNodeFactory::instance()->register_factory("random_orient", make_grid_variant_operator< RandomOrientTmpl >);
  }

}
