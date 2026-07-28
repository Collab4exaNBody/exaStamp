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
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  // Moves particle positions along the conjugate-gradient direction h (vx,vy,vz):
  //   r += inv_xform * h * (alpha / dir_max)
  // alpha is the requested maximum atomic displacement for this trial (a length), dir_max is
  // max|h| as computed by cg_direction, so the largest single-atom displacement this step
  // produces is exactly alpha.
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz >
    >
  struct CGTrialStep : public OperatorNode
  {
    ADD_SLOT( GridT  , grid    , INPUT_OUTPUT );
    ADD_SLOT( Domain , domain  , INPUT , REQUIRED );
    ADD_SLOT( double , alpha   , INPUT , REQUIRED );
    ADD_SLOT( double , dir_max , INPUT , REQUIRED );

    inline void execute () override final
    {
      if( *dir_max <= 0.0 ) return; // no direction (e.g. already at a stationary point)

      GridT& grid = *(this->grid);
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();
      auto cells = grid.cells();

      const Mat3d inv_xform = domain->inv_xform();
      const double scale = (*alpha) / (*dir_max);

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          auto* __restrict__ rx = cells[i][field::rx];
          auto* __restrict__ ry = cells[i][field::ry];
          auto* __restrict__ rz = cells[i][field::rz];
          const auto* __restrict__ hx = cells[i][field::vx];
          const auto* __restrict__ hy = cells[i][field::vy];
          const auto* __restrict__ hz = cells[i][field::vz];
#         pragma omp simd
          for(size_t j=0;j<n;j++)
          {
            const Vec3d d = inv_xform * ( Vec3d{ hx[j], hy[j], hz[j] } * scale );
            rx[j] += d.x;
            ry[j] += d.y;
            rz[j] += d.z;
          }
        }
        GRID_OMP_FOR_END
      }
    }
  };

  template<class GridT> using CGTrialStepTmpl = CGTrialStep<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(cg_trial_step)
  {
   OperatorNodeFactory::instance()->register_factory( "cg_trial_step", make_grid_variant_operator< CGTrialStepTmpl > );
  }

}
