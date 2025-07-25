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
#include <onika/physics/constants.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/unit_system.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // CoulWolf Parameters
  struct CoulWolfParms
  {
    double alpha = 0.0;
    double rc = 0.0;
    double qqrd2e = 14.399645;
    double e_shift = 0.0;
    double f_shift = 0.0;
  };
  
  // core computation kernel for coul wolf potential
  ONIKA_HOST_DEVICE_FUNC inline void coul_wolf_compute_energy(const CoulWolfParms& p, double c1, double c2, double r, double& e, double& de)
  {
    assert( r > 0. );
    
    // LAMMPS
    double prefactor = p.qqrd2e * c1 * c2 / r;
    double erfcc = erfc(p.alpha * r);
    double erfcd = exp(-p.alpha * p.alpha * r * r);
    double v_sh = (erfcc - p.e_shift * r) * prefactor;
    double dvdrr = (erfcc / ( r * r ) + 2.0 * p.alpha / sqrt(M_PI) * erfcd / r) + p.f_shift;
    double forcecoul = dvdrr * r * r * prefactor;
    double fpair = -forcecoul / r;
    
    e = EXASTAMP_QUANTITY( v_sh * eV );
    de = EXASTAMP_QUANTITY( fpair * eV / ang );
    
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::CoulWolfParms;
  
  using onika::physics::Quantity;

  template<> struct convert<CoulWolfParms>
  {
    static bool decode(const Node& node, CoulWolfParms& v)
    {
      if( !node.IsMap() ) { return false; }
      v.alpha = node["alpha"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      v.e_shift = erfc(v.alpha * v.rc) / v.rc;
      v.f_shift = -(v.e_shift + 2.0 * v.alpha / sqrt(M_PI) * exp(-v.alpha * v.alpha * v.rc * v.rc)) / v.rc;
      
      return true;
    }
  };
}

