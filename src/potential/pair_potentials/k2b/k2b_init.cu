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

#include <algorithm>

#include "potential.h"

// Registers the k2b cutoff in init_parameters, so that rcut_max is known before setup_system builds
// ghosts and neighbor lists. parameters and rcut are passed through for compute_descriptor_k2b.
namespace exaStamp
{
  using namespace exanb;

  class K2bInit : public OperatorNode
  {
    ADD_SLOT( K2bPotentialParameters , parameters , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double                 , rcut       , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double                 , rcut_max   , INPUT_OUTPUT , 0.0 );

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Provides the k2b parameters and cutoff to compute_descriptor_k2b and raises rcut_max to rcut.
Place it in init_parameters, after species.

Usage example:

init_parameters:
  - species
  - k2b_init:
      rcut: 6.0 ang
      parameters: { n_rbf: 8, r_min: 0.5, r_cut: 6.0, sigma: 0.4, delta: 1.5, w: [ -1.057, -2.0949, 0.90561, -2.56538, 0.21529, -0.80587, -2.65201, 0.04461 ] }

)EOF";
    }

    inline void execute() override final
    {
      ldbg << "Initializing k2b potential, rcut=" << *rcut << std::endl;
      *rcut_max = std::max( *rcut_max, *rcut );
    }

  };

  ONIKA_AUTORUN_INIT(k2b_init)
  {
    OperatorNodeFactory::instance()->register_factory("k2b_init", make_simple_operator<K2bInit>);
  }

}
