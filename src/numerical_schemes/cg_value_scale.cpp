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
#include <algorithm>
#include <limits>

namespace exaStamp
{
  using namespace exanb;

  // Generic scalar rescale with clamping: value = clamp( value * factor , min , max ), in place.
  // Used to grow/shrink the conjugate-gradient line-search step size (two separate instances,
  // rebound to different factor/min/max values) without a dedicated operator for each direction.
  // Named cg_value_scale (not value_scale) to avoid colliding with exaNBody's own
  // value_scale_operator.cpp, which registers "value_scale" with a different (non-in-place,
  // non-clamping) in_value/out_value signature.
  class CGValueScaleNode : public OperatorNode
  {
    ADD_SLOT( double , value  , INPUT_OUTPUT );
    ADD_SLOT( double , factor , INPUT , REQUIRED );
    ADD_SLOT( double , min    , INPUT , std::numeric_limits<double>::lowest() );
    ADD_SLOT( double , max    , INPUT , std::numeric_limits<double>::max() );

  public:
    inline void execute () override final
    {
      double v = (*value) * (*factor);
      v = std::max( *min , std::min( *max , v ) );
      *value = v;
      ldbg << "CGValueScale: value="<<v<<std::endl;
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(cg_value_scale)
  {
    OperatorNodeFactory::instance()->register_factory( "cg_value_scale", make_compatible_operator< CGValueScaleNode > );
  }

}
