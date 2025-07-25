/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

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

