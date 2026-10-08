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

// Per atom charge front-end of the Wolf potential. The kernel is defined once, in the coul_wolf pair potential.
#include <exaStamp/potential/pair_potentials/coul_wolf/coul_wolf.h>

namespace exaStamp
{
  using WolfParameters = CoulWolfParms;

  struct WolfKernel
  {
    WolfParameters m_params;
    WolfKernel() = default;
    inline WolfKernel(const WolfParameters& p) : m_params(p) {}
    static inline const char* documentation() { return R"EOF(
Wolf damped shifted coulomb potential with per particle charges (LAMMPS pair_style coul/wolf) :
E = qi.qj/(4.pi.eps0) [erfc(alpha.r)/r - erfc(alpha.rc)/rc] for r < rc, force shifted to 0 at rc.
parameters : { alpha: <1/distance> , rc: <distance> }. The self energy -(e_shift/2 + alpha/sqrt(pi)).q^2/(4.pi.eps0) is
included (self_energy: true). Pair style with species charges : coul_wolf.
)EOF"; }
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double c, double r, double& e, double& de) const { coul_wolf_kernel( m_params, c, r, e, de ); }
    ONIKA_HOST_DEVICE_FUNC inline double self_energy(double q) const { return coul_wolf_self_energy( m_params, q ); }
  };
}
