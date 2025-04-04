#pragma once

#include <cmath>
#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential/ewald/ewald.h>
#include <onika/flat_tuple.h>

namespace exaStamp
{
  using namespace exanb;

  template<class EwaldParmsT>
  ONIKA_HOST_DEVICE_FUNC static inline void ewald_compute_energy(const EwaldParmsT& p, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    ewald_compute_energy( p, p_pair.m_atom_a.m_charge, p_pair.m_atom_b.m_charge, r, e, de );
  }
}

#define USTAMP_POTENTIAL_NAME     ewald_short_range
#define USTAMP_POTENTIAL_PARAMS   EwaldParms
#define USTAMP_POTENTIAL_COMPUTE  ewald_compute_energy
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) onika::make_flat_tuple(p.m_atom_a.m_charge,p.m_atom_b.m_charge)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1

