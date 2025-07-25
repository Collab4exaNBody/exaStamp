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

#pragma once

#include <cmath>
#include <utility>

#include <onika/cuda/cuda.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>

namespace exaStamp
{
  using namespace exanb;

  ONIKA_HOST_DEVICE_FUNC inline void ljexp6rf_energy(const LJExp6RFParms& p_rc, const PairPotentialMinimalParameters& p_pair, double r, double& _e, double& _de)
  {
    const auto [ de , e ] = p_rc.compute_force_energy ( r , p_pair.m_atom_a.m_charge * p_pair.m_atom_b.m_charge );
    _e = e; _de = de;
  }
}

#define USTAMP_POTENTIAL_NAME     ljexp6rf
#define USTAMP_POTENTIAL_PARAMS   LJExp6RFParms
#define USTAMP_POTENTIAL_COMPUTE  ljexp6rf_energy

// only atom charges are meaningful for LJExp6RF
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) std::make_pair(p.m_atom_a.m_charge,p.m_atom_b.m_charge)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1
#define USTAMP_POTENTIAL_ENABLE_RIGIDMOL 1  // define to 1 to generate the rigid molecule variant

