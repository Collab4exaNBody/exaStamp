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

#include <memory>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <onika/parallel/random.h>
#include <exanb/grid_cell_particles/particle_region.h>
#include <exaStamp/unit_system.h>
#include <exaStamp/compute/thermodynamic_state.h>

#include <memory>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  class VirialFdotRNode : public OperatorNode
  {
    
    ADD_SLOT(MPI_Comm       , mpi          , INPUT, MPI_COMM_WORLD);
    ADD_SLOT(GridT          , grid         , INPUT_OUTPUT);
    ADD_SLOT(Mat3d          , stress_tensor, OUTPUT);
    ADD_SLOT(double         , volume      , INPUT);
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      if( grid->number_of_cells() == 0 ) return;
      GridT& grid              = *(this->grid);
      auto cells = grid.cells();
      IJK dims = grid.dimension();

      const ThermodynamicState& sim_info = *thermodynamic_state;
      Mat3d stress = make_zero_matrix();
      // ReduceStressTensorFunctor func;
      // reduce_cell_particles(*grid, false, func, stress, reduce_field_set, parallel_execution_context(), {}, rcpo); 

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:stress) )
        {
          
          size_t cell_i = grid_ijk_to_index( dims , cell_loc );

          const auto* __restrict__ rx = cells[cell_i][field::rx];
          const auto* __restrict__ ry = cells[cell_i][field::ry];
          const auto* __restrict__ rz = cells[cell_i][field::rz];

          const auto* __restrict__ fx = cells[cell_i][field::fx];
          const auto* __restrict__ fy = cells[cell_i][field::fy];
          const auto* __restrict__ fz = cells[cell_i][field::fz];
          
          const unsigned int n = cells[cell_i].size();

          for(unsigned int j=0;j<n;j++)
          {
            stress += tensor( Vec3d{fx[j],fy[j],fz[j]}, Vec3d{rx[j],ry[j],rz[j]} );
          }

        }
        GRID_OMP_FOR_END
	    }

      {
        double buff[9] = {stress.m11, stress.m12, stress.m13, stress.m21, stress.m22, stress.m23, stress.m31, stress.m32, stress.m33};
        MPI_Allreduce(MPI_IN_PLACE, buff, 9, MPI_DOUBLE, MPI_SUM, *mpi);
        stress = Mat3d{buff[0], buff[1], buff[2], buff[3], buff[4], buff[5], buff[6], buff[7], buff[8]} / *volume;
        *stress_tensor = stress;
        
        Mat3d virial_classic = sim_info.stress_tensor();
        std::cout << "FDOTR version = " << stress << std::endl;
        std::cout << "CLASS version = " << virial_classic << std::endl;        
      //   double tmp[6] = { sum_Te, sum_dTe, sum_Se, sum_Si };
      //   MPI_Allreduce(MPI_IN_PLACE,tmp,4,MPI_DOUBLE,MPI_SUM,*mpi);
      //   sum_Te = tmp[0];
      //   sum_dTe = tmp[1];
      //   sum_Se = tmp[2];
      //   sum_Si = tmp[3];
      }

    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
doc TODO
)EOF";
    }

  };

  template<class GridT> using VirialFdotRNodeTmpl = VirialFdotRNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(virial_fdotr)
  {
   OperatorNodeFactory::instance()->register_factory(
    "virial_fdotr",
    make_grid_variant_operator< VirialFdotRNodeTmpl >
    );
  }

}
