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

#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <exaStamp/potential/coulombic/ewald.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  class EwaldInitOperator : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( double     , accuracy_relative , INPUT , REQUIRED );
    ADD_SLOT( double     , g_ewald       , INPUT , REQUIRED );
    ADD_SLOT( double     , radius      , INPUT , REQUIRED );
    ADD_SLOT( long       , kmax        , INPUT , REQUIRED );
    ADD_SLOT( Domain     , domain      , INPUT , OPTIONAL );
    ADD_SLOT( EwaldParameters , ewald_config, INPUT_OUTPUT );
    ADD_SLOT( double     , rcut        , OUTPUT );
    ADD_SLOT( double     , sum_square_charge, INPUT );
    ADD_SLOT( double     , sum_charge, INPUT );
    ADD_SLOT( uint64_t   , natoms , INPUT , OPTIONAL);
    ADD_SLOT( double     , rcut_max    , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT );
             
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute () override final
    {

      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);
      
      using ewald_constants::fpe0;
      using ewald_constants::epsilonZero;
      
      MPI_Barrier(*mpi);
      
      if( domain.has_value() )
      {
        auto & p = *ewald_config;
        
        // Check if initialization is needed based on input parameters (not current state)
        // This ensures all ranks make the same decision
        bool need_init = false;
        
        // Always initialize if volume is zero (first time)
        if( p.volume == 0.0 )
        {
          need_init = true;
        }
        // Or if any input parameter has changed
        else if( (*g_ewald > 0.0 && *g_ewald != p.g_ewald) || 
                 *radius != p.radius || 
                 *accuracy_relative != p.accuracy_relative || 
                 (*kmax > 0 && *kmax != p.kmax) )
        {
          need_init = true;
        }
        
        // Synchronize the decision across all ranks to ensure consistency
        int need_init_flag = need_init ? 1 : 0;
        int need_init_global = 0;
        MPI_Allreduce(&need_init_flag, &need_init_global, 1, MPI_INT, MPI_MAX, *mpi);
        need_init = (need_init_global != 0);
        
        if( need_init )
        {
          auto domainSize = domain->bounds_size();
          auto xform = domain->xform();
          if( ! is_diagonal( xform ) )
          {
            lerr << "Domain XForm is not diagonal, cannot compute domain box size" << std::endl;
            std::abort();
          }
          domainSize = xform * domainSize;
          if( ! ( domain->periodic_boundary_x() && domain->periodic_boundary_y() && domain->periodic_boundary_z() ) )
          {
            lerr << "Domain must be entierly periodic, cannot initialize ewald." << std::endl;
            std::abort();
          }

          // Initialize parameters on rank 0
          if( rank == 0 )
          {
            ewald_init_parameters( *g_ewald , *radius , *accuracy_relative , *kmax , domainSize, *natoms, *sum_square_charge, *sum_charge, p , ldbg<<"" );
          }
          
          // Broadcast all scalar parameters from rank 0 to all other ranks
          MPI_Bcast(&p.g_ewald, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.radius, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.accuracy_relative, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.kmax, 1, MPI_LONG_LONG_INT, 0, *mpi);
          MPI_Bcast(&p.kxmax, 1, MPI_LONG_LONG_INT, 0, *mpi);
          MPI_Bcast(&p.kymax, 1, MPI_LONG_LONG_INT, 0, *mpi);
          MPI_Bcast(&p.kzmax, 1, MPI_LONG_LONG_INT, 0, *mpi);
          MPI_Bcast(&p.nk, 1, MPI_LONG_LONG_INT, 0, *mpi);
          MPI_Bcast(&p.nknz, 1, MPI_LONG_LONG_INT, 0, *mpi);
          MPI_Bcast(&p.gm, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.gm_sr, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.qqr2e, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.bt_sr, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.volume, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.unitk.x, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.unitk.y, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.unitk.z, 1, MPI_DOUBLE, 0, *mpi);
          
          // Broadcast the approximation constants (though these should be the same)
          MPI_Bcast(&p.EWALD_F, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.EWALD_P, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.A1, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.A2, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.A3, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.A4, 1, MPI_DOUBLE, 0, *mpi);
          MPI_Bcast(&p.A5, 1, MPI_DOUBLE, 0, *mpi);
          
          // Resize and broadcast Gdata vector
          if( rank != 0 )
          {
            p.Gdata.resize(p.nknz);
          }
          
          if( p.nknz > 0 )
          {
            MPI_Bcast(p.Gdata.data(), p.nknz * sizeof(EwaldCoeffs), MPI_BYTE, 0, *mpi);
          }
          
          // Only print from rank 0 to avoid duplicate output
          if( rank == 0 && p.volume > 0.0 )
          {
            lout << "====== Ewald configuration ======" << std::endl;      
            lout << "size    = "<< domainSize << std::endl;
            lout << "period. = "<< std::boolalpha << domain->periodic_boundary_x() 
                                << " , " << std::boolalpha << domain->periodic_boundary_y() 
                                << " , " << std::boolalpha << domain->periodic_boundary_z() << std::endl;                    
            lout << "g_ewald = "<<p.g_ewald << std::endl;
            lout << "radius  = "<<p.radius << std::endl;
            lout << "accuracy_relative = "<<p.accuracy_relative << std::endl;
            lout << "kmax    = "<<p.kmax << std::endl;
            lout << "nk      = "<<p.nk << std::endl;
            lout << "nknz    = "<<p.nknz << std::endl;
            lout << "gm      = "<<p.gm << std::endl;
            lout << "lm      = "<< p.unitk.x <<','<<p.unitk.y<<','<<p.unitk.z << std::endl;
            lout << "volume  = "<< p.volume << std::endl;
            lout << "=================================" << std::endl;
          }
        }
      }
      // -----------------------------------------------
      // -----------------------------------------------

      *rcut = *radius;
      *rcut_max = std::max( *rcut_max , *radius );
      
      MPI_Barrier(*mpi);
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_ewald_init)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_ewald_init" , make_simple_operator<EwaldInitOperator> );
  }

}
