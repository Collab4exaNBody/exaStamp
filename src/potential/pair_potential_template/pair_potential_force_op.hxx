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

#include <onika/log.h>
#include <onika/math/basic_types_operators.h>

#include "pair_potential_force_op.h"

 // do we need pair Weighting version ?
#if ! defined(USTAMP_POTENTIAL_WITH_WEIGHTS) || EXASTAMP_ENABLE_PAIR_WEIGHTING

#ifdef USTAMP_POTENTIAL_WITH_WEIGHTS
#define UseWeights true
#else
#define UseWeights false
#endif

namespace exaStamp
{
  using namespace exanb;

  void _CLASS_NAME::operator() (
    ComputePairBuffer2<UseWeights,false>& tab,
    double& ep,
    double& ax,
    double& ay,
    double& az
#   ifdef USTAMP_POTENTIAL_WITH_VIRIAL
    , Mat3d& virial
#   endif
    ) const noexcept
  {
    const auto & p = m_potential_params;
    const auto & pair_params = m_pair_params;
    const unsigned int n = tab.count;
    const double ecut = m_ecut;

#   ifndef USTAMP_POTENTIAL_WITH_VIRIAL
    Mat3d virial;
#   endif

#   include "force_op_impl2.hxx"
  }
  
}

#undef _CLASS_NAME
#undef _CLASS_BASE
#undef UseWeights

#endif // do we need pair Weighting version ?


