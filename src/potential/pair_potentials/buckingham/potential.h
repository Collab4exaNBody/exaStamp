#pragma once

#include "buckingham.h"

#define USTAMP_POTENTIAL_NAME     buckingham
#define USTAMP_POTENTIAL_PARAMS   BuckinghamParms
#define USTAMP_POTENTIAL_COMPUTE  buckingham_energy

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

//#define USTAMP_POTENTIAL_ENABLE_CUDA 1

