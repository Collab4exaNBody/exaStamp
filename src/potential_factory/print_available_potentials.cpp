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

#include <exaStamp/potential_factory/pair_potential_factory.h>

namespace exaStamp
{
  using namespace exanb;

  struct PrintAvailablePotentials : public OperatorNode
  {  
    ADD_SLOT( std::string , header , INPUT , std::string("")   , DocString{"Message header to print"} );

    // -----------------------------------------------
    // ----------- Operator documentation ------------
    inline std::string documentation() const override final
    {
      return R"EOF(
        Prints available pair potentials to the standard output.
        )EOF";
    }

    inline void execute () override final
    {
      lout << *header ;
      for( const auto & s : PairPotentialFactory::available_potentials() )
      {
        lout << s << std::endl;
      }
    }

  };
    
  // === register factories ===  
  ONIKA_AUTORUN_INIT(print_available_potentials)
  {
    OperatorNodeFactory::instance()->register_factory( "print_available_potentials", make_simple_operator<PrintAvailablePotentials> );
  }

}

