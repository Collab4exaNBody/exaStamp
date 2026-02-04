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

#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>

#define USTAMP_POTENTIAL_NAME         lj
#define USTAMP_POTENTIAL_PARAMS       LennardJonesParms
#define USTAMP_POTENTIAL_COMPUTE      lj_compute_energy

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

#define USTAMP_POTENTIAL_ENABLE_CUDA 1      // define to 1 to enable Cuda/HIP kernel compilation
#define USTAMP_POTENTIAL_ENABLE_RIGIDMOL 1  // define to 1 to generate the rigid molecule variant

