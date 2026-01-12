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

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute () override final
    {
      
      using ewald_constants::fpe0;
      using ewald_constants::epsilonZero;

      if( domain.has_value() )
      {
        auto & p = *ewald_config;
        
        if( ( *g_ewald > 0.0 && *g_ewald != p.g_ewald ) || p.volume==0.0 || *radius != p.radius || *accuracy_relative != p.accuracy_relative || ( *kmax > 0 && *kmax != p.kmax ) )
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

          ewald_init_parameters( *g_ewald , *radius , *accuracy_relative , *kmax , domainSize, *natoms, *sum_square_charge, *sum_charge, p , ldbg<<"" );
          if( p.volume > 0.0 )
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
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_ewald_init)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_ewald_init" , make_simple_operator<EwaldInitOperator> );
  }

}


