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

  // Reduces per-particle force (fx,fy,fz, aliased to ax,ay,az) into two scalars used to drive
  // conjugate-gradient minimization: the max force magnitude (convergence test) and the sum of
  // squared force magnitudes (Fletcher-Reeves beta numerator/denominator).
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz >
    >
  struct CGForceStats : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi          , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT    , grid         , INPUT );
    ADD_SLOT( double   , force_max    , OUTPUT );
    ADD_SLOT( double   , force_sqnorm , OUTPUT );

    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();
      auto cells = grid.cells();

      double local_max = 0.0;
      double local_sqnorm = 0.0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) reduction(max:local_max) reduction(+:local_sqnorm) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          const auto* __restrict__ fx = cells[i][field::fx];
          const auto* __restrict__ fy = cells[i][field::fy];
          const auto* __restrict__ fz = cells[i][field::fz];
          for(size_t j=0;j<n;j++)
          {
            const double f2 = fx[j]*fx[j] + fy[j]*fy[j] + fz[j]*fz[j];
            local_sqnorm += f2;
            local_max = std::max( local_max , std::sqrt(f2) );
          }
        }
        GRID_OMP_FOR_END
      }

      double global_max = local_max;
      double global_sqnorm = local_sqnorm;
      MPI_Allreduce( MPI_IN_PLACE , &global_max    , 1 , MPI_DOUBLE , MPI_MAX , *mpi );
      MPI_Allreduce( MPI_IN_PLACE , &global_sqnorm , 1 , MPI_DOUBLE , MPI_SUM , *mpi );

      *force_max = global_max;
      *force_sqnorm = global_sqnorm;

      ldbg << "CGForceStats: force_max="<<global_max<<", force_sqnorm="<<global_sqnorm<<std::endl;
    }
  };

  template<class GridT> using CGForceStatsTmpl = CGForceStats<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(cg_force_stats)
  {
   OperatorNodeFactory::instance()->register_factory( "cg_force_stats", make_grid_variant_operator< CGForceStatsTmpl > );
  }

}
