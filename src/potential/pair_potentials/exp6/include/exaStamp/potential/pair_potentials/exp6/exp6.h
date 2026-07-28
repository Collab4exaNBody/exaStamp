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

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;


  // Exp6 Parameters
  struct Exp6Parms
  {
    double A = 0.0;
    double B = 1.0;
    double C = 0.0;
    double D = 0.0;
  };

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::Exp6Parms;
  using onika::physics::Quantity;

  template<> struct convert<Exp6Parms>
  {
    static inline bool decode(const Node& node, Exp6Parms& v)
    {
      if( !node.IsMap() ) { return false; }

      v.A = node["A"].as<Quantity>().convert();
      v.B = node["B"].as<Quantity>().convert();
      v.C = node["C"].as<Quantity>().convert();
      v.D = node["D"].as<Quantity>().convert();

      return true;
    }
  };
}

namespace exaStamp
{
  using namespace exanb;

  ONIKA_HOST_DEVICE_FUNC inline void exp6_compute_energy(const Exp6Parms& p_exp6, const PairPotentialMinimalParameters&, double r, double& e, double& de)
  {
    assert(r>0);
    double one_rB = 1/(r*p_exp6.B);
    // pow(x,6)/pow(x,12) with an integer exponent still lowers to a full
    // exp(y*log(x)) sequence on the SFU on CUDA - replaced with exact
    // repeated squaring, same as lj_compute_energy.
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double Cr6 = p_exp6.C / r6;
    double twelve_rB = 12 * one_rB;
    double twelve_rB2 = twelve_rB * twelve_rB;
    double twelve_rB4 = twelve_rB2 * twelve_rB2;
    double twelve_rB12 = twelve_rB4 * twelve_rB4 * twelve_rB4;
    double Dtwelve_rB12 = p_exp6.D * twelve_rB12;
    double AexpmBr = p_exp6.A * exp(-p_exp6.B * r);

    e  =   AexpmBr       - Cr6   + Dtwelve_rB12;
    de = - p_exp6.B * AexpmBr + (6 * Cr6 - 12 * Dtwelve_rB12)/r;
  }

}


