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
#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential/ewald/ewald.h>
#include <onika/flat_tuple.h>

namespace exaStamp
{
  using namespace exanb;

  template<class EwaldParmsT>
  ONIKA_HOST_DEVICE_FUNC static inline void ewald_compute_energy(const EwaldParmsT& p, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    ewald_compute_energy( p, p_pair.m_atom_a.m_charge, p_pair.m_atom_b.m_charge, r, e, de );
  }
}

#define USTAMP_POTENTIAL_NAME     ewald_short_range
#define USTAMP_POTENTIAL_PARAMS   EwaldParms
#define USTAMP_POTENTIAL_COMPUTE  ewald_compute_energy
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) onika::make_flat_tuple(p.m_atom_a.m_charge,p.m_atom_b.m_charge)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1

