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
#include <mpi.h>

namespace exaStamp
{
inline namespace coulombic_ewald // distinct symbols from the legacy ewald plugin (plugins are loaded RTLD_GLOBAL)
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
      int rank=0;
      MPI_Comm_rank(*mpi, &rank);
      
      if( domain.has_value() )
      {
        auto & p = *ewald_config;

        auto domainSize = domain->bounds_size();
        const auto xform = domain->xform();
        if( ! is_diagonal( xform ) )
        {
          fatal_error() << "Domain XForm is not diagonal, cannot compute domain box size" << std::endl;
        }
        domainSize = xform * domainSize;

        // called before the system is set up (e.g. in init_parameters, to get rcut) : nothing to initialize yet
        const bool empty_domain = ( domainSize.x * domainSize.y * domainSize.z ) == 0.0;

        // all inputs are global (identical on all ranks), so every rank takes the same decision and
        // computes the same parameters : no broadcast needed.
        const bool need_init = ! empty_domain && (
             p.volume == 0.0
          || ( *g_ewald > 0.0 && *g_ewald != p.g_ewald )
          || *radius != p.radius
          || *accuracy_relative != p.accuracy_relative
          || ( *kmax > 0 && *kmax != p.kmax )
          || domainSize != p.box ); // box changed (NPT, deformation) : k vectors must be rebuilt
        
        if( need_init )
        {
          if( ! ( domain->periodic_boundary_x() && domain->periodic_boundary_y() && domain->periodic_boundary_z() ) )
          {
            fatal_error() << "Domain must be entierly periodic, cannot initialize ewald." << std::endl;
          }
          if( ! natoms.has_value() )
          {
            fatal_error() << "coulombic_ewald_init : natoms is required" << std::endl;
          }

          const bool first_init = ( p.volume == 0.0 );
          ewald_init_parameters( *g_ewald , *radius , *accuracy_relative , *kmax , domainSize, *natoms, *sum_square_charge, *sum_charge, p , ldbg<<"" );
          
          if( rank == 0 && p.volume > 0.0 && first_init )
          {
            lout << "====== Ewald configuration ======" << std::endl;      
            lout << "size    = "<< domainSize << std::endl;
            lout << "g_ewald = "<<p.g_ewald << std::endl;
            lout << "radius  = "<<p.radius << std::endl;
            lout << "accuracy_relative = "<<p.accuracy_relative << std::endl;
            lout << "kmax    = "<<p.kmax << " ("<<p.kxmax<<","<<p.kymax<<","<<p.kzmax<<")" << std::endl;
            lout << "nknz    = "<<p.nknz << std::endl;
            lout << "qsum    = "<<p.qsum << std::endl;
            lout << "volume  = "<< p.volume << std::endl;
            lout << "=================================" << std::endl;
          }
          else
          {
            ldbg << "Ewald re-initialized : size="<<domainSize<<" g_ewald="<<p.g_ewald<<" nknz="<<p.nknz<<std::endl;
          }
        }
      }

      *rcut = *radius;
      *rcut_max = std::max( *rcut_max , *radius );
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_ewald_init)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_ewald_init" , make_simple_operator<EwaldInitOperator> );
  }

}
}
