#pragma once

#include "coul_cut.h"

#define USTAMP_POTENTIAL_NAME     coul_cut
#define USTAMP_POTENTIAL_PARAMS   CoulCutParms
#define USTAMP_POTENTIAL_COMPUTE  coul_cut_energy

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

//#define USTAMP_POTENTIAL_ENABLE_CUDA 1

