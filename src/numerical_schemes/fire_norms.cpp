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

namespace exaStamp
{
  using namespace exanb;

  // Force and velocity squared-norms (summed over all atoms), computed mid-step per FIRE 2.0's
  // Velocity-Verlet algorithm (Guenole et al. 2020, Algorithm 6): called right after the first
  // half-kick, so force_sqnorm is |F(t)|^2 (pre-step force, unchanged since the half-kick only
  // touches velocity) and veloc_sqnorm is |v(t+dt/2)|^2 (post-half-kick velocity) - together they
  // give fire_velocity_mix's |v|/|F| mixing ratio. The end-of-step force_max used for the
  // convergence test and progress display is computed separately, by reusing cg_force_stats.
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_vx, field::_vy, field::_vz >
    >
  struct FIRENorms : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi          , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT    , grid         , INPUT );
    ADD_SLOT( double   , force_sqnorm , OUTPUT );
    ADD_SLOT( double   , veloc_sqnorm , OUTPUT );

    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();
      auto cells = grid.cells();

      double local_force_sqnorm = 0.0;
      double local_veloc_sqnorm = 0.0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) reduction(+:local_force_sqnorm) reduction(+:local_veloc_sqnorm) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          const auto* __restrict__ fx = cells[i][field::fx];
          const auto* __restrict__ fy = cells[i][field::fy];
          const auto* __restrict__ fz = cells[i][field::fz];
          const auto* __restrict__ vx = cells[i][field::vx];
          const auto* __restrict__ vy = cells[i][field::vy];
          const auto* __restrict__ vz = cells[i][field::vz];
          for(size_t j=0;j<n;j++)
          {
            local_force_sqnorm += fx[j]*fx[j] + fy[j]*fy[j] + fz[j]*fz[j];
            local_veloc_sqnorm += vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j];
          }
        }
        GRID_OMP_FOR_END
      }

      double global_force_sqnorm = local_force_sqnorm;
      double global_veloc_sqnorm = local_veloc_sqnorm;
      MPI_Allreduce( MPI_IN_PLACE , &global_force_sqnorm , 1 , MPI_DOUBLE , MPI_SUM , *mpi );
      MPI_Allreduce( MPI_IN_PLACE , &global_veloc_sqnorm , 1 , MPI_DOUBLE , MPI_SUM , *mpi );

      *force_sqnorm = global_force_sqnorm;
      *veloc_sqnorm = global_veloc_sqnorm;

      ldbg << "FIRENorms: force_sqnorm="<<global_force_sqnorm<<", veloc_sqnorm="<<global_veloc_sqnorm<<std::endl;
    }
  };

  template<class GridT> using FIRENormsTmpl = FIRENorms<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(fire_norms)
  {
   OperatorNodeFactory::instance()->register_factory( "fire_norms", make_grid_variant_operator< FIRENormsTmpl > );
  }

}
