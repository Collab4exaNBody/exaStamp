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
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>
#include <exaStamp/coulomb_constant.h>
#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // Damped shifted force coulomb potential (LAMMPS pair coul/dsf).
  // Single definition, used by the pair potential template (coul_dsf) and by coulombic_dsf (per atom charges).
  // Note : the pair potential template subtracts e(rcut) from every pair. e(rc) is ~1.8e-7 eV per unit charge product here
  // (the A&S erfc is not exact), so coul_dsf energies differ from LAMMPS by that constant per pair (forces are identical).
  // coulombic_dsf does not shift and matches LAMMPS energies.
  struct CoulDsfParms
  {
    double alpha = 0.0;
    double rc = 0.0;
    double qqrd2e = COULOMB_CONSTANT_EV_ANG;
    double e_shift = 0.0;
    double f_shift = 0.0;
  };

  // pair term for charge product c = qi.qj : e (energy) and de = de/dr, internal units.
  // As in LAMMPS, erfc(alpha.r) uses the Abramowitz & Stegun approximation, the shifts are computed once with exact erfc
  ONIKA_HOST_DEVICE_FUNC inline void coul_dsf_kernel(const CoulDsfParms& p, double c, double r, double& e, double& de)
  {
    assert( r > 0. );
    constexpr double EWALD_P = 0.3275911;
    constexpr double A1 = 0.254829592;
    constexpr double A2 = -0.284496736;
    constexpr double A3 = 1.421413741;
    constexpr double A4 = -1.453152027;
    constexpr double A5 = 1.061405429;
    const double MY_PIS = sqrt(M_PI);
    const double rsq = r*r;

    const double prefactor = p.qqrd2e * c / r;
    const double erfcd = exp(-p.alpha * p.alpha * rsq);
    const double t = 1.0 / (1.0 + EWALD_P * p.alpha * r);
    const double erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;

    const double forcecoul = prefactor * (erfcc / r + 2.0 * p.alpha / MY_PIS * erfcd + r * p.f_shift) * r;
    const double fpair = -forcecoul / r;
    const double ecoul = prefactor * (erfcc - r * p.e_shift - rsq * p.f_shift);
    
    e = EXASTAMP_QUANTITY( ecoul * eV );
    de = EXASTAMP_QUANTITY( fpair * eV / ang );
  }

  // self energy of a particle with charge q (LAMMPS pair coul/dsf e_self), internal units.
  // Not included in the pair term : use coulombic_dsf_self.
  ONIKA_HOST_DEVICE_FUNC inline double coul_dsf_self_energy(const CoulDsfParms& p, double q)
  {
    const double e_self = -( p.e_shift / 2.0 + p.alpha / sqrt(M_PI) ) * q * q * p.qqrd2e;
    return EXASTAMP_QUANTITY( e_self * eV );
  }

  // pair potential template adapter
  ONIKA_HOST_DEVICE_FUNC inline void coul_dsf_pair_energy(const CoulDsfParms& p, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    coul_dsf_kernel( p, p_pair.m_atom_a.m_charge * p_pair.m_atom_b.m_charge, r, e, de );
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::CoulDsfParms>
  {
    static bool decode(const Node& node, exaStamp::CoulDsfParms& v)
    {
      using onika::physics::Quantity;
      if( !node.IsMap() ) { return false; }
      v.alpha = node["alpha"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      const double erfcc = erfc(v.alpha * v.rc);
      const double erfcd = exp(-v.alpha * v.alpha * v.rc * v.rc);
      v.f_shift = -(erfcc / (v.rc * v.rc) + 2.0 / sqrt(M_PI) * v.alpha * erfcd / v.rc);
      v.e_shift = erfcc / v.rc - v.f_shift * v.rc;
      return true;
    }
  };
}
