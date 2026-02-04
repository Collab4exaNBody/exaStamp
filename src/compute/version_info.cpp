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
#include <onika/log.h>

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  class ApplicationVersionInfo : public OperatorNode
  {  
    ADD_SLOT( bool , show_details, INPUT , true , DocString{"Print additional details"}  );

  public:
    inline void execute () override final
    {
      lout << "exaStamp v"<< EXASTAMP_VERSION << " ("
#     ifndef NDEBUG
      << "debug"
#     else
      << "release"
#     endif
      <<")" << std::endl;
      if( *show_details )
      {
        int mpi_maj=0, mpi_min=0;
        MPI_Get_version( &mpi_maj, &mpi_min );
        lout << " , MPI v"<<mpi_maj<<'.'<<mpi_min << std::endl;
      }
    }
  };
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(version_info)
  {
   OperatorNodeFactory::instance()->register_factory( "version_info", make_simple_operator< ApplicationVersionInfo > );
  }

}

