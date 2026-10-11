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

// Per atom charge front-end of the reaction field potential. The kernel is defined once, in the coul_rf pair potential
// (also used by the ljrf, exp6rf and ljexp6rf pair potentials and the ljexp6rf_pc molecular force field).
#include <exaStamp/potential/pair_potentials/coul_rf/coul_rf.h>

namespace exaStamp
{
  using ReactionFieldParameters = ReactionFieldParms;

  struct ReactionFieldKernel
  {
    ReactionFieldParameters m_params;
    ReactionFieldKernel() = default;
    inline ReactionFieldKernel(const ReactionFieldParameters& p) : m_params(p) {}
    static inline const char* documentation() { return R"EOF(
Reaction field coulomb potential with per particle charges (dielectric continuum beyond rc) :
E = qi.qj/(4.pi.eps0) [1/r + k_rf.r^2 - c_rf] for r < rc, k_rf = (eps-1)/(2.eps+1)/rc^3, c_rf = 1/rc + k_rf.rc^2.
parameters : { epsilon: <dielectric constant> , rc: <distance> }. Pair style with species charges : coul_rf.
)EOF"; }
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double c, double r, double& e, double& de) const { reaction_field_compute_energy( m_params, c, r, e, de ); }
  };
}
