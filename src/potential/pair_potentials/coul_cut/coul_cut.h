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
#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // Coul_cut Parameters
  struct CoulCutParms
  {
    double dielectric = 1.0;
    bool shift = true;
  };

  ONIKA_HOST_DEVICE_FUNC inline void coul_cut_energy(const CoulCutParms& p, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    static constexpr double qqr2e = EXASTAMP_CONST_QUANTITY( 14.399645 * eV * ang / (ec^2) ) ; // units = eV.ang/e-^2 from LAMMPS

    assert( r > 0. );

    const double C = qqr2e / p.dielectric;

    e = C * p_pair.m_atom_a.m_charge * p_pair.m_atom_b.m_charge / r;
    de = -e/r;
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::CoulCutParms>
  {
    static bool decode(const Node& node, exaStamp::CoulCutParms& v)
    {
      using onika::physics::Quantity;

      if( !node.IsMap() ) { return false; }
      v.dielectric = node["dielectric"].as<Quantity>().convert();
      v.shift = node["shift"].as<bool>();
      std::cout << "shift val = " << v.shift << std::endl;
      return true;
    }
  };
}

