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

  // FIRE's power P = F.v, summed over all atoms (fx,fy,fz aliased to ax,ay,az; vx,vy,vz hold
  // the real velocity here, unlike the CG scheme which repurposes them for a search direction).
  // Sign of P drives the whole FIRE adaptation state machine (fire_adapt): P>0 means the
  // trajectory is still going downhill and dt/alpha keep adapting, P<=0 means it just went
  // uphill and needs correcting.
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_vx, field::_vy, field::_vz >
    >
  struct FIREPower : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi   , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT    , grid  , INPUT );
    ADD_SLOT( double   , power , OUTPUT );

    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();
      auto cells = grid.cells();

      double local_power = 0.0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) reduction(+:local_power) )
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
            local_power += fx[j]*vx[j] + fy[j]*vy[j] + fz[j]*vz[j];
          }
        }
        GRID_OMP_FOR_END
      }

      double global_power = local_power;
      MPI_Allreduce( MPI_IN_PLACE , &global_power , 1 , MPI_DOUBLE , MPI_SUM , *mpi );

      *power = global_power;

      ldbg << "FIREPower: power="<<global_power<<std::endl;
    }
  };

  template<class GridT> using FIREPowerTmpl = FIREPower<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(fire_power)
  {
   OperatorNodeFactory::instance()->register_factory( "fire_power", make_grid_variant_operator< FIREPowerTmpl > );
  }

}
