#pragma once

#include <cmath>
#include <utility>

#include <onika/cuda/cuda.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>

namespace exaStamp
{
  using namespace exanb;

  ONIKA_HOST_DEVICE_FUNC inline void ljexp6rf_energy(const LJExp6RFParms& p_rc, const PairPotentialMinimalParameters& p_pair, double r, double& _e, double& _de)
  {
    const auto [ de , e ] = p_rc.compute_force_energy ( r , p_pair.m_atom_a.m_charge * p_pair.m_atom_b.m_charge );
    _e = e; _de = de;
  }
}

#define USTAMP_POTENTIAL_NAME     ljexp6rf
#define USTAMP_POTENTIAL_PARAMS   LJExp6RFParms
#define USTAMP_POTENTIAL_COMPUTE  ljexp6rf_energy

// only atom charges are meaningful for LJExp6RF
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) std::make_pair(p.m_atom_a.m_charge,p.m_atom_b.m_charge)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1
#define USTAMP_POTENTIAL_ENABLE_RIGIDMOL 1  // define to 1 to generate the rigid molecule variant

