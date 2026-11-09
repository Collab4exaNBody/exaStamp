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

#include <exaStamp/unit_system.h>

namespace exaStamp
{
  // Coulomb constant 1/(4.pi.epsilon0), same value as LAMMPS metal units (force->qqr2e), so that coulombic
  // potentials (ewald, wolf, dsf, coul_cut) give the same results as LAMMPS.
  // Note : the reaction field potentials use onika::physics::epsilonZero instead (14.3996454784 eV.ang/e-^2).
  static constexpr double COULOMB_CONSTANT_EV_ANG = 14.399645; // eV.ang/e-^2
  static constexpr double COULOMB_CONSTANT = EXASTAMP_CONST_QUANTITY( COULOMB_CONSTANT_EV_ANG * eV * ang / (ec^2) ); // internal units
}
