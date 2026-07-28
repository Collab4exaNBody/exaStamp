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

#include <algorithm>
#include <memory>

namespace exaStamp
{
  using namespace exanb;

  // FIRE 2.0 adaptive dt/alpha state machine (Guenole et al. 2020, Algorithm 2, lines 8-30 -
  // everything except the MD integration itself, which the enclosing .msp loop does with the
  // ordinary push_f_v_r/push_f_v Verlet operators rebound to this operator's own dt).
  //
  // On a downhill step (power>0): after >delaystep consecutive downhill steps, grow dt (capped
  // at dt_start*dt_max_factor) and shrink alpha (multiplicatively towards 0, biasing v towards
  // pure inertial motion as the trajectory settles into a smooth descent).
  // On an uphill/stationary step (power<=0): reset the downhill counter, and - unless still
  // within the initial_delay grace period - shrink dt (floored at dt_start*dt_min_factor) and
  // reset alpha to alpha_start. Either way, "uphill" is reported so the caller can apply the
  // half-step-back position correction and zero the velocity.
  // fire_not_stalled goes false once too many consecutive uphill steps (pneg_max) have piled up
  // without a single downhill step in between - a sign further relaxation isn't possible - and
  // feeds into the trajectory loop's own continue condition, same role as CG's line_search_ok.
  class FIREAdaptNode : public OperatorNode
  {
    ADD_SLOT( double , power              , INPUT , REQUIRED );

    ADD_SLOT( double , dt_start           , INPUT , REQUIRED );
    ADD_SLOT( double , dt_max_factor      , INPUT , 10.0 );
    ADD_SLOT( double , dt_min_factor      , INPUT , 0.02 );
    ADD_SLOT( long   , delaystep          , INPUT , 20 );
    ADD_SLOT( double , dt_grow            , INPUT , 1.1 );
    ADD_SLOT( double , dt_shrink          , INPUT , 0.5 );
    ADD_SLOT( double , alpha_start        , INPUT , 0.25 );
    ADD_SLOT( double , alpha_shrink       , INPUT , 0.99 );
    ADD_SLOT( long   , pneg_max           , INPUT , 2000 );
    ADD_SLOT( bool   , initial_delay      , INPUT , true );

    ADD_SLOT( double , dt                 , INPUT_OUTPUT );
    ADD_SLOT( double , alpha              , INPUT_OUTPUT );
    ADD_SLOT( long   , n_ppos             , INPUT_OUTPUT , 0 );
    ADD_SLOT( long   , n_pneg             , INPUT_OUTPUT , 0 );
    ADD_SLOT( long   , iter               , INPUT_OUTPUT , 0 );

    ADD_SLOT( bool   , uphill             , OUTPUT );
    ADD_SLOT( bool   , fire_not_stalled   , OUTPUT );

  public:
    inline void execute () override final
    {
      *iter = *iter + 1;
      const double dt_max = (*dt_start) * (*dt_max_factor);
      const double dt_min = (*dt_start) * (*dt_min_factor);

      const bool up = ( *power <= 0.0 );
      *uphill = up;

      if( ! up )
      {
        *n_ppos = *n_ppos + 1;
        *n_pneg = 0;
        if( *n_ppos > *delaystep )
        {
          *dt = std::min( (*dt) * (*dt_grow) , dt_max );
          *alpha = (*alpha) * (*alpha_shrink);
        }
        *fire_not_stalled = true;
      }
      else
      {
        *n_ppos = 0;
        *n_pneg = *n_pneg + 1;
        *fire_not_stalled = ( *n_pneg <= *pneg_max );

        const bool skip_shrink = (*initial_delay) && ( *iter < *delaystep );
        if( ! skip_shrink )
        {
          if( (*dt) * (*dt_shrink) >= dt_min )
          {
            *dt = (*dt) * (*dt_shrink);
          }
          *alpha = *alpha_start;
        }
      }

      ldbg << "FIREAdapt: power="<<(*power)<<", uphill="<<up<<", dt="<<(*dt)<<", alpha="<<(*alpha)
           <<", n_ppos="<<(*n_ppos)<<", n_pneg="<<(*n_pneg)<<", not_stalled="<<(*fire_not_stalled)<<std::endl;
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(fire_adapt)
  {
    OperatorNodeFactory::instance()->register_factory( "fire_adapt", make_compatible_operator< FIREAdaptNode > );
  }

}
