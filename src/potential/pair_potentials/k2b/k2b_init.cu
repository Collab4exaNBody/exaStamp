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

// Early (init_parameters, before setup_system) k2b rcut_max registration -- mirrors pod_init.cu /
// mtp_init.cu / snap_init.cu, but needs no context object of its own: K2bPotentialParameters is
// already a plain, already-parsed struct with no derived/cached engine state to build (k2b's
// descriptor math is computed inline, no external library like POD/MTP's). `parameters`/`rcut` are
// plain INPUT_OUTPUT pass-through slots (same idiom as rebo_init.cu's own `parameters` slot) --
// compute_descriptor_k2b.cu keeps reading them by the same names, now just supplied earlier, so it
// needs zero source changes. Only job here: set rcut_max in time for ghost/neighbor setup, removing
// the manual `rcut_max: X ang` workaround every k2b .msp previously needed.
namespace exaStamp
{
  using namespace exanb;

  class K2bInit : public OperatorNode
  {
    ADD_SLOT( K2bPotentialParameters , parameters , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double                 , rcut       , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double                 , rcut_max   , INPUT_OUTPUT , 0.0 );

  public:

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
