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

#include "pair_potential_template.h"

#include <exaStamp/potential_factory/pair_potential_factory.h>
#include <exaStamp/potential_factory/pair_potential.h>

#undef USTAMP_POTENTIAL_WITH_VIRIAL
#include "pair_potential_force_op.h"

#define USTAMP_POTENTIAL_WITH_VIRIAL 1
#include "pair_potential_force_op.h"

#define POTENTIAL_REGISTER_INIT() _POTENTIAL_REGISTER_INIT( CONSTRUCTOR_FUNC_NAME )
#define CONSTRUCTOR_FUNC_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_CLASS,_init)
#define _POTENTIAL_REGISTER_INIT(name) CONSTRUCTOR_ATTRIB void MAKE_UNIQUE_NAME(name,_,__LINE__,ONIKA_CURRENT_PACKAGE_NAME) ()

namespace exaStamp
{
  using namespace exanb;


  class USTAMP_POTENTIAL_CLASS final : public PairPotential
  {
  public:

    // constructor
    USTAMP_POTENTIAL_CLASS(const USTAMP_POTENTIAL_PARAMS& p) : m_p(p) {}
    
    // factory
    static std::shared_ptr<PairPotential> make_potential(YAML::Node node)
    {
      return std::make_shared<USTAMP_POTENTIAL_CLASS>( node.as<USTAMP_POTENTIAL_PARAMS>() );
    }

    // specialized for multimaterial, called from traversal method for each atom pair type
    std::shared_ptr<PairPotentialComputeOperator> force_op() override final
    {
      return std::make_shared< USTAMP_POTENTIAL_FORCE_OP_NAME >(m_p);
    }

    std::shared_ptr<PairPotentialComputeVirialOperator> force_virial_op() override final
    {
      return std::make_shared< USTAMP_POTENTIAL_FORCE_VIRIAL_OP_NAME >(m_p);
    }

  private:
    const USTAMP_POTENTIAL_PARAMS m_p;
  };

  // === register potential factory ===  
  // ONIKA_AUTORUN_INIT(pair_potential_class)
  POTENTIAL_REGISTER_INIT()
  {
    PairPotentialFactory::register_factory( USTAMP_POTENTIAL_STRING , USTAMP_POTENTIAL_CLASS::make_potential );
  }

}


