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

#include <exanb/core/grid.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <exanb/core/parallel_grid_algorithm.h>

#include <iostream>
#include <mpi.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_charge >
    >
  struct SumChargesPCOperator : public OperatorNode
  {      
    static constexpr size_t SIMD_VECTOR_SIZE = GridT::CellParticles::ChunkSize ;

    // ========= I/O slots =======================
    ADD_SLOT( MPI_Comm, mpi, INPUT);
    ADD_SLOT( GridT  , grid              , INPUT );
    ADD_SLOT( double , sum_charge        , OUTPUT );
    ADD_SLOT( double , sum_square_charge , OUTPUT );

    // Operator execution
    inline void execute () override final
    {
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();

      double sc = 0.0;
      double sc2 = 0.0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, reduction(+:sc) reduction(+:sc2) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          const auto* __restrict__ charges = cells[i][field::charge]; ONIKA_ASSUME_ALIGNED(charges);
          double lsc = 0.0;
          double lsc2 = 0.0;


#         pragma omp simd reduction(+:lsc) reduction(+:lsc2)
          for(size_t j=0;j<n;j++)
          {
            double c = charges[j];
            lsc += c;
            lsc2 += c*c;
          }
          sc += lsc;
          sc2 += lsc2;
        }
        GRID_OMP_FOR_END
      }
    
      {
        double tmp[2] = { sc , sc2 };
        MPI_Allreduce(MPI_IN_PLACE,tmp,2,MPI_DOUBLE,MPI_SUM,*mpi);

       *sum_charge = tmp[0];
       *sum_square_charge = tmp[1]; 
      }

  //    *sum_charge = sc;
 //     *sum_square_charge = sc2;


      ldbg<<" Sum charge : "<<*sum_charge<<std::endl<<std::flush;
      ldbg<<" Sum square charge : "<<*sum_square_charge<<std::endl<<std::flush;

      if(*sum_charge > 1.e-10){
        lout<<"######## ERROR in sum_charges_pc.cpp : the system is'nt neutral : sum charge :  ="<<*sum_charge<<std::endl;
        std::abort();
      }

    }

  };

  template<class GridT> using SumChargesPC = SumChargesPCOperator<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(sum_charges_pc)
  {  
    OperatorNodeFactory::instance()->register_factory( "sum_charges_pc" , make_grid_variant_operator< SumChargesPC > );
  }

}


