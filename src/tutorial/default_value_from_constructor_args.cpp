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
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>

#include <vector>

namespace exaStamp
{
  using namespace exanb;

  static const std::vector<double> myvec = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };

  class DefaultValueFromCTorArgs : public OperatorNode
  {      
    // ========= I/O slots =======================
    ADD_SLOT( std::vector<double> , input1, INPUT, std::make_tuple(size_t(10),5.0) );
    ADD_SLOT( std::vector<double> , input2, INPUT, std::make_tuple(myvec.begin()+1,myvec.begin()+6) );

  public:
    // Operator execution
    inline void execute () override final
    {      
      lout<<"input1 : length = " << input1->size()<<", content =";
      for(const auto& x:*input1) { lout << " " << x; }
      lout << std::endl;

      lout<<"input2 : length = " << input2->size()<<", content =";
      for(const auto& x:*input2) { lout << " " << x; }
      lout << std::endl;
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(default_value_from_constructor_args)
  {  
    OperatorNodeFactory::instance()->register_factory( "default_slot_value_from_ctor_args" , make_simple_operator< DefaultValueFromCTorArgs > );
  }

}


