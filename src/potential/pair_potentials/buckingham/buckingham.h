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
#include <onika/cuda/cuda.h>
#include <exaStamp/potential_factory/pair_potential.h>

namespace exaStamp
{
  using namespace exanb;

  // Buckingham Parameters
  struct BuckinghamParms
  {
    double A = 0.0;
    double Rho = 0.0;
    double C = 0.0;
  };

  ONIKA_HOST_DEVICE_FUNC inline void buckingham_energy(const BuckinghamParms& p, const PairPotentialMinimalParameters&, double x, double& e, double& de)
  {
    // derivative from Wolfram : derivative of A*exp(-x/R)-(C/x^6)
    // A => p.A , R => p.Rho , C => p.C
    // https://www.wolframalpha.com/input/?i=derivative+of+A*exp%28-x%2FR%29-%28C%2Fx%5E6%29
    assert( x > 0. );
    double x6 = pow(x,6);
    double x7 = pow(x,7);
    e = p.A * exp( -x / p.Rho ) - ( p.C / x6 );
    de = ( 6 * p.C / x7 ) - ( p.A * std::exp( -x / p.Rho ) / p.Rho );
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::BuckinghamParms;
  using onika::physics::Quantity;

  template<> struct convert<BuckinghamParms>
  {
    static bool decode(const Node& node, BuckinghamParms& v)
    {
      if( !node.IsMap() ) { return false; }
      v.A = node["A"].as<Quantity>().convert();
      v.Rho = node["Rho"].as<Quantity>().convert();
      v.C = node["C"].as<Quantity>().convert();
      return true;
    }
  };
}

