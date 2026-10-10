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

// Per atom charge front-end of the DSF potential. The kernel is defined once, in the coul_dsf pair potential.
#include <exaStamp/potential/pair_potentials/coul_dsf/coul_dsf.h>

namespace exaStamp
{
  using DsfParameters = CoulDsfParms;

  struct DsfKernel
  {
    DsfParameters m_params;
    DsfKernel() = default;
    inline DsfKernel(const DsfParameters& p) : m_params(p) {}
    static inline const char* documentation() { return R"EOF(
Damped shifted force coulomb potential with per particle charges (Fennell & Gezelter) :
energy and force shifted to 0 at rc. parameters : { alpha: <1/distance> , rc: <distance> }.
The self energy -(e_shift/2 + alpha/sqrt(pi)).q^2/(4.pi.eps0) is included (self_energy: true).
Pair style with species charges : coul_dsf.
)EOF"; }
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double c, double r, double& e, double& de) const { coul_dsf_kernel( m_params, c, r, e, de ); }
    ONIKA_HOST_DEVICE_FUNC inline double self_energy(double q) const { return coul_dsf_self_energy( m_params, q ); }
  };
}
