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

#include <cmath>
#include <algorithm>
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
    ADD_SLOT( uint64_t       , natoms           , OUTPUT , DocString{"global number of particles"} );

    // Operator execution
    inline void execute () override final
    {
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();

      double sc = 0.0;
      double sc2 = 0.0;
      size_t np = 0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, reduction(+:sc) reduction(+:sc2) reduction(+:np) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          np += n;
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
        unsigned long long n_all = np;
        MPI_Allreduce(MPI_IN_PLACE,&n_all,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,*mpi);
        *natoms = n_all;

       *sum_charge = tmp[0];
       *sum_square_charge = tmp[1]; 
      }



      ldbg<<" Sum charge : "<<*sum_charge<<std::endl<<std::flush;
      ldbg<<" Sum square charge : "<<*sum_square_charge<<std::endl<<std::flush;

      // non neutral system : Ewald and PPPM add the neutralizing background term, short range methods do not
      int rank = 0;
      MPI_Comm_rank( *mpi , &rank );
      if( rank == 0 && std::abs( *sum_charge ) > 1.e-10 * std::max( 1.0 , std::sqrt( *sum_square_charge ) ) )
      {
        lout << "Warning: the system is not neutral, total charge = "<< *sum_charge
             << " e. Ewald and PPPM add a neutralizing background term, short range coulombic methods do not." << std::endl;
      }

    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Total charge, sum of squared charges and number of particles of the system (all MPI ranks), computed from the per particle charge field.
Required by coulombic_ewald_init and coulombic_pppm_init (call it at the end of setup_system, before them).
Use sum_charges when charges come from the species definitions instead.
Prints a warning when the system is not neutral.
)EOF";
    }
  };

  template<class GridT> using SumChargesPC = SumChargesPCOperator<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(sum_charges_pc)
  {  
    OperatorNodeFactory::instance()->register_factory( "sum_charges_pc" , make_grid_variant_operator< SumChargesPC > );
  }

}


