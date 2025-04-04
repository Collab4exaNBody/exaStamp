#pragma once

#include <cmath>
#include <utility>

#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential/reaction_field/reaction_field.h>

#include <onika/cuda/cuda.h>
#include <onika/flat_tuple.h>

namespace exaStamp
{
  using namespace exanb;

  ONIKA_HOST_DEVICE_FUNC inline void reaction_field_compute_energy(const ReactionFieldParms& p_rc, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    reaction_field_compute_energy( p_rc, p_pair.m_atom_a.m_charge * p_pair.m_atom_b.m_charge, r, e, de );
  }
}

#define USTAMP_POTENTIAL_NAME     reaction_field
#define USTAMP_POTENTIAL_PARAMS   ReactionFieldParms
#define USTAMP_POTENTIAL_COMPUTE  reaction_field_compute_energy

// only atom charges are meaningful for reaction field
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) onika::make_flat_tuple(p.m_atom_a.m_charge,p.m_atom_b.m_charge)

#define USTAMP_POTENTIAL_ENABLE_CUDA 1

