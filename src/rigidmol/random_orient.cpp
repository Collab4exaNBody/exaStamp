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
#include <exanb/core/domain.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_orient >
    >
  class RandomOrient : public OperatorNode
  {
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( Domain          , domain         , INPUT );
    ADD_SLOT( bool            , deterministic_noise , INPUT , false );

  public:
    inline void execute () override final
    {
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();      
      IJK gstart { gl, gl, gl };
      IJK gend = dims - IJK{ gl, gl, gl };
      IJK gdims = gend - gstart;
      const auto dom_dims = domain->grid_dimension();
      const auto dom_start = grid->offset();

      const int nthreads = *deterministic_noise ? 1 : omp_get_max_threads();
#     pragma omp parallel num_threads(nthreads)
      {
        std::mt19937_64 det_re;
        std::mt19937_64 & re = *deterministic_noise ? det_re : onika::parallel::random_engine() ;
        std::normal_distribution<double> f_rand(-1.,1.);
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          const auto i = grid_ijk_to_index( dims , loc + gstart );
          const auto domain_cell_idx = grid_ijk_to_index( dom_dims , loc + gstart + dom_start );
          
          auto* __restrict__ orient = cells[i][field::orient];
          size_t n = cells[i].size();
          det_re.seed( domain_cell_idx * 1023 );          
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
