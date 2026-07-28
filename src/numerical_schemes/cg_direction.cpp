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

#include <mpi.h>
#include <memory>
#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  // Updates the conjugate-gradient search direction h (stored in the unused vx,vy,vz fields,
  // since this is a static minimization with no velocity/time integration) using
  // h_new = F + beta * h_old, and reduces two scalars needed by the caller:
  //   dir_max  = max |h_new| over all particles (used to normalize the trial step length)
  //   slope    = -sum(F . h_new) / dir_max       (directional derivative dE/dalpha at alpha=0,
  //                                                used by the Armijo line-search condition)
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx,field::_fy,field::_fz, field::_vx,field::_vy,field::_vz >
    >
  struct CGDirection : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi     , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT    , grid    , INPUT_OUTPUT );
    ADD_SLOT( double   , beta    , INPUT , REQUIRED );
    ADD_SLOT( double   , dir_max , OUTPUT );
    ADD_SLOT( double   , slope   , OUTPUT );

    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();
      auto cells = grid.cells();
      const double beta_val = *beta;

      double local_max = 0.0;
      double local_dot = 0.0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) reduction(max:local_max) reduction(+:local_dot) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          const auto* __restrict__ fx = cells[i][field::fx];
          const auto* __restrict__ fy = cells[i][field::fy];
          const auto* __restrict__ fz = cells[i][field::fz];
          auto* __restrict__ hx = cells[i][field::vx];
          auto* __restrict__ hy = cells[i][field::vy];
          auto* __restrict__ hz = cells[i][field::vz];
          for(size_t j=0;j<n;j++)
          {
            const double Fx = fx[j], Fy = fy[j], Fz = fz[j];
            const double Hx = Fx + beta_val * hx[j];
            const double Hy = Fy + beta_val * hy[j];
            const double Hz = Fz + beta_val * hz[j];
            hx[j] = Hx; hy[j] = Hy; hz[j] = Hz;
            local_dot += Fx*Hx + Fy*Hy + Fz*Hz;
            local_max = std::max( local_max , std::sqrt( Hx*Hx + Hy*Hy + Hz*Hz ) );
          }
        }
        GRID_OMP_FOR_END
      }

      double global_max = local_max;
      double global_dot = local_dot;
      MPI_Allreduce( MPI_IN_PLACE , &global_max , 1 , MPI_DOUBLE , MPI_MAX , *mpi );
      MPI_Allreduce( MPI_IN_PLACE , &global_dot , 1 , MPI_DOUBLE , MPI_SUM , *mpi );

      *dir_max = global_max;
      *slope = ( global_max > 0.0 ) ? ( -global_dot / global_max ) : 0.0;

      ldbg << "CGDirection: beta="<<beta_val<<", dir_max="<<global_max<<", slope="<<(*slope)<<std::endl;
    }
  };

  template<class GridT> using CGDirectionTmpl = CGDirection<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(cg_direction)
  {
   OperatorNodeFactory::instance()->register_factory( "cg_direction", make_grid_variant_operator< CGDirectionTmpl > );
  }

}
