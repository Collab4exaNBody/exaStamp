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
#include <exanb/core/domain.h>

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  class DomainVolume : public OperatorNode
  {  
    ADD_SLOT( MPI_Comm, mpi, INPUT, MPI_COMM_WORLD, DocString{"MPI communicator"} );
    ADD_SLOT( Domain  , domain, INPUT, REQUIRED  , DocString{"Simulation domain"});
    ADD_SLOT( double  , volume, INPUT_OUTPUT      , DocString{"Computed volume"}  );

  public:
    inline bool is_sink() const override final { return true; }

    inline void execute () override final
    {
       *volume = domain->bounds().bmax.x - domain->bounds().bmin.x;
       lout << "Volume = " << *volume << std::endl;
    }

    inline std::string documentation() const override final
    {
      return "This is the documentation of domain_volume component. (c) 2023 T. Carrard";
    }
  };
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(domain_volume)
  {
   OperatorNodeFactory::instance()->register_factory( "domain_volume", make_simple_operator< DomainVolume > );
  }

}

