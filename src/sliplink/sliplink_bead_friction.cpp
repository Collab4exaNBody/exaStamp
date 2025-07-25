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
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <onika/memory/allocator.h>
#include <onika/parallel/random.h>
#include <exaStamp/sliplink/sliplink.h>

#include <memory>


namespace exaStamp
{

  template< class GridT >
  struct SliplinkForceOverdamped : public OperatorNode
  {      
    ADD_SLOT( GridT , grid ,INPUT_OUTPUT);
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT, REQUIRED );

    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      
      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();

      const double sigma1 = sliplink_config->sigma1;

#     pragma omp parallel
      {
        auto& re = onika::parallel::random_engine();
        std::normal_distribution<double> gaussian(0.0,sigma1);

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc)
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          
          auto* __restrict__ rx = cells[i][ field::rx ]; ONIKA_ASSUME_ALIGNED(rx);
          auto* __restrict__ ry = cells[i][ field::ry ]; ONIKA_ASSUME_ALIGNED(ry);
          auto* __restrict__ rz = cells[i][ field::rz ]; ONIKA_ASSUME_ALIGNED(rz);

          for(size_t k=0;k<n;k++)
          {
			rx[k] += gaussian(re);
			ry[k] += gaussian(re);
			rz[k] += gaussian(re);
          }
        }
        GRID_OMP_FOR_END
      }
    }

  };

  template<class GridT> using SliplinkForceOverdampedTmpl = SliplinkForceOverdamped<GridT>;
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(sliplink_bead_friction)
  {
   OperatorNodeFactory::instance()->register_factory( "sliplink_bead_friction", make_grid_variant_operator< SliplinkForceOverdampedTmpl > );
  }

}

