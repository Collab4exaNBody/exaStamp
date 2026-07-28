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
#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  // FIRE velocity mixing: v <- (1-alpha)*v + alpha*F*(|v|/|F|), the same "mixing factor" ratio
  // (|v|/|F|, a single global scalar from fire_norms) applied to every atom's vector - this is
  // an Euler discretization of the bias term in Bitzek et al.'s equation of motion, steering v
  // towards the force direction without changing its magnitude.
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_vx, field::_vy, field::_vz >
    >
  struct FIREVelocityMix : public OperatorNode
  {
    ADD_SLOT( GridT  , grid         , INPUT_OUTPUT );
    ADD_SLOT( double , alpha        , INPUT , REQUIRED );
    ADD_SLOT( double , force_sqnorm , INPUT , REQUIRED );
    ADD_SLOT( double , veloc_sqnorm , INPUT , REQUIRED );

    inline void execute () override final
    {
      if( *force_sqnorm <= 0.0 ) return; // no force (stationary point): mixing is undefined, skip

      const double scale = std::sqrt( (*veloc_sqnorm) / (*force_sqnorm) );
      const double a = *alpha;
      const double one_minus_a = 1.0 - a;

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
          const auto* __restrict__ fx = cells[i][field::fx];
          const auto* __restrict__ fy = cells[i][field::fy];
          const auto* __restrict__ fz = cells[i][field::fz];
          auto* __restrict__ vx = cells[i][field::vx];
          auto* __restrict__ vy = cells[i][field::vy];
          auto* __restrict__ vz = cells[i][field::vz];
#         pragma omp simd
          for(size_t j=0;j<n;j++)
          {
            vx[j] = one_minus_a * vx[j] + a * fx[j] * scale;
            vy[j] = one_minus_a * vy[j] + a * fy[j] * scale;
            vz[j] = one_minus_a * vz[j] + a * fz[j] * scale;
          }
        }
        GRID_OMP_FOR_END
      }
    }
  };

  template<class GridT> using FIREVelocityMixTmpl = FIREVelocityMix<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(fire_velocity_mix)
  {
   OperatorNodeFactory::instance()->register_factory( "fire_velocity_mix", make_grid_variant_operator< FIREVelocityMixTmpl > );
  }

}
