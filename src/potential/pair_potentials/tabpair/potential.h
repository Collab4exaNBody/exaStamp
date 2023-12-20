#pragma once

#include "tabulated_pair.h"

#define USTAMP_POTENTIAL_NAME     tabpair
#define USTAMP_POTENTIAL_PARAMS   TabPairPotentialParms
#define USTAMP_POTENTIAL_COMPUTE  tab_pair_potential_energy
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

//#define USTAMP_POTENTIAL_ENABLE_CUDA 1

