#pragma once

#include <exaStamp/potential/pair_potentials/exp6/exp6.h>

#define USTAMP_POTENTIAL_NAME     exp6
#define USTAMP_POTENTIAL_PARAMS   Exp6Parms
#define USTAMP_POTENTIAL_COMPUTE  exp6_compute_energy
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

#define USTAMP_POTENTIAL_ENABLE_CUDA 1

