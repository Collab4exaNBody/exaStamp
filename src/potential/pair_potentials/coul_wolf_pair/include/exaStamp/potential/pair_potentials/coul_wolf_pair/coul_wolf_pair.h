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

  // Wolf damped shifted coulomb potential (LAMMPS pair coul/wolf).
  // Single definition, used by the pair potential template (coul_wolf_pair, ljwolf) and by coulombic_wolf (per atom charges).
  struct CoulWolfParms
  {
    double alpha = 0.0;
    double rc = 0.0;
    double qqrd2e = COULOMB_CONSTANT_EV_ANG;
    double e_shift = 0.0;
    double f_shift = 0.0;    
  };

  // pair term for charge product c = qi.qj : e (energy) and de = de/dr, internal units
  ONIKA_HOST_DEVICE_FUNC inline void coul_wolf_kernel(const CoulWolfParms& p, double c, double r, double& e, double& de)
  {
    assert( r > 0. );
    
    const double prefactor = p.qqrd2e * c / r;
    const double r2 = r * r;
    const double alpha2 = p.alpha * p.alpha;
    
    const double erfcc = erfc(p.alpha * r);
    const double erfcd = exp(-alpha2 * r2);
    const double v_sh = (erfcc - p.e_shift * r) * prefactor;

    e = EXASTAMP_QUANTITY( v_sh * eV );
    
    const double dvdrr = (erfcc / r2 + 2.0 * p.alpha / sqrt(M_PI) * erfcd / r) + p.f_shift;
    const double forcecoul = dvdrr * r2 * prefactor;
    const double fpair = -forcecoul / r;
    
    de = EXASTAMP_QUANTITY( fpair * eV / ang );
  }

  // self energy of a particle with charge q (LAMMPS pair coul/wolf e_self), internal units.
  // Not included in the pair term : use coulombic_wolf_self.
  ONIKA_HOST_DEVICE_FUNC inline double coul_wolf_self_energy(const CoulWolfParms& p, double q)
  {
    const double e_self = -( p.e_shift / 2.0 + p.alpha / sqrt(M_PI) ) * q * q * p.qqrd2e;
    return EXASTAMP_QUANTITY( e_self * eV );
  }

  // pair potential template adapter
  ONIKA_HOST_DEVICE_FUNC inline void coul_wolf_pair_energy(const CoulWolfParms& p, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    coul_wolf_kernel( p, p_pair.m_atom_a.m_charge * p_pair.m_atom_b.m_charge, r, e, de );
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::CoulWolfParms>
  {
    static bool decode(const Node& node, exaStamp::CoulWolfParms& v)
    {
      using onika::physics::Quantity;
      if( !node.IsMap() ) { return false; }
      v.alpha = node["alpha"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      v.e_shift = erfc(v.alpha * v.rc) / v.rc;
      v.f_shift = -(v.e_shift + 2.0 * v.alpha / sqrt(M_PI) * exp(-v.alpha * v.alpha * v.rc * v.rc)) / v.rc;      
      return true;
    }
  };
}
