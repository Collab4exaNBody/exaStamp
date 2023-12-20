#pragma once

#include <exaStamp/potential/pair_potentials/coul_wolf_pair/coul_wolf_pair.h>

#define USTAMP_POTENTIAL_NAME     coul_wolf_pair
#define USTAMP_POTENTIAL_PARAMS   CoulWolfParms
#define USTAMP_POTENTIAL_COMPUTE  coul_wolf_pair_energy

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

//#define USTAMP_POTENTIAL_ENABLE_CUDA 1

