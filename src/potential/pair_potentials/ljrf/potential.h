#pragma once

#include "ljrf.h"

#define USTAMP_POTENTIAL_NAME         ljrf
#define USTAMP_POTENTIAL_PARAMS       LJRFParms
#define USTAMP_POTENTIAL_COMPUTE      ljrf_energy

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) std::make_pair(p.m_atom_a.m_charge,p.m_atom_b.m_charge)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1
#define USTAMP_POTENTIAL_ENABLE_RIGIDMOL 1  // define to 1 to generate the rigid molecule variant

