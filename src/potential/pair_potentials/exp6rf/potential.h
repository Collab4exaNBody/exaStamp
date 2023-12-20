#pragma once

#include "exp6rf.h"

#define USTAMP_POTENTIAL_NAME         exp6rf
#define USTAMP_POTENTIAL_PARAMS       Exp6RFParms
#define USTAMP_POTENTIAL_COMPUTE      exp6rf_energy

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) std::make_pair(p.m_atom_a.m_charge,p.m_atom_b.m_charge)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1

