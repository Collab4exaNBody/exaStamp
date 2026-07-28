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

  // Fletcher-Reeves conjugate-gradient coefficient: beta = ||F_new||^2 / ||F_old||^2 .
  // force_sqnorm_prev persists across outer CG iterations (same operator instance is re-run
  // each iteration by the enclosing loop) and doubles as the "no previous iteration yet" flag
  // via its 0.0 default, giving beta=0 (plain steepest descent) on the first iteration.
  //
  // Periodic restart (beta forced to 0, i.e. a fresh steepest-descent direction) every
  // restart_period iterations: plain Fletcher-Reeves search directions degrade after many
  // iterations without one, a well-known FR-CG failure mode that shows up as force/energy
  // plateauing well above the requested tolerance instead of continuing to converge.
  class CGBetaNode : public OperatorNode
  {
    ADD_SLOT( double , force_sqnorm       , INPUT , REQUIRED );
    ADD_SLOT( double , force_sqnorm_prev  , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( long   , restart_period     , INPUT , 20 );
    ADD_SLOT( long   , iters_since_restart, INPUT_OUTPUT , 0 );
    ADD_SLOT( double , beta               , OUTPUT );

  public:
    inline void execute () override final
    {
      bool restart = ( *force_sqnorm_prev <= 0.0 ) || ( *iters_since_restart >= *restart_period );
      *beta = restart ? 0.0 : ( *force_sqnorm / *force_sqnorm_prev );
      *force_sqnorm_prev = *force_sqnorm;
      *iters_since_restart = restart ? 1 : ( *iters_since_restart + 1 );
      ldbg << "CGBeta: beta="<<(*beta)<<", restart="<<restart<<std::endl;
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(cg_beta)
  {
    OperatorNodeFactory::instance()->register_factory( "cg_beta", make_compatible_operator< CGBetaNode > );
  }

}
