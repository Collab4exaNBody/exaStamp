#pragma once

#include "ljwolf.h"

#define USTAMP_POTENTIAL_NAME         ljwolf
#define USTAMP_POTENTIAL_PARAMS       LJWOLFParms
#define USTAMP_POTENTIAL_COMPUTE      ljwolf_energy

//#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) std::make_pair(p.m_atom_a.m_charge,p.m_atom_b.m_charge)
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1

