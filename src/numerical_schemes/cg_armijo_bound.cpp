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

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  // Armijo sufficient-decrease bound: bound = ref_energy + c1 * alpha * slope
  // (slope < 0 for a descent direction, so bound < ref_energy). A trial step is accepted
  // when its energy is <= this bound.
  class CGArmijoBoundNode : public OperatorNode
  {
    ADD_SLOT( double , ref_energy , INPUT , REQUIRED );
    ADD_SLOT( double , alpha      , INPUT , REQUIRED );
    ADD_SLOT( double , slope      , INPUT , REQUIRED );
    ADD_SLOT( double , c1         , INPUT , 1.0e-4 );
    ADD_SLOT( double , bound      , OUTPUT );

  public:
    inline void execute () override final
    {
      *bound = *ref_energy + (*c1) * (*alpha) * (*slope);
      ldbg << "CGArmijoBound: bound="<<(*bound)<<std::endl;
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(cg_armijo_bound)
  {
    OperatorNodeFactory::instance()->register_factory( "cg_armijo_bound", make_compatible_operator< CGArmijoBoundNode > );
  }

}
