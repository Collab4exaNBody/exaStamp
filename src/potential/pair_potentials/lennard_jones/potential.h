#pragma once

#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>

#define USTAMP_POTENTIAL_NAME         lj
#define USTAMP_POTENTIAL_PARAMS       LennardJonesParms
#define USTAMP_POTENTIAL_COMPUTE      lj_compute_energy

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

#define USTAMP_POTENTIAL_ENABLE_CUDA 1      // define to 1 to enable Cuda/HIP kernel compilation
#define USTAMP_POTENTIAL_ENABLE_RIGIDMOL 1  // define to 1 to generate the rigid molecule variant

