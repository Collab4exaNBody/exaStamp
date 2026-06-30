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
#include <onika/cuda/cuda.h>
#include <exaStamp/unit_system.h>
#include <exaStamp/potential_factory/pair_potential.h>

namespace exaStamp
{
  using namespace exanb;

  // Yukawa Parameters
  struct YukawaParms
  {
    double A = 0.0;
    double kappa = 0.0;
  };

  ONIKA_HOST_DEVICE_FUNC
  inline void yukawa_compute_energy(const YukawaParms& p, const PairPotentialMinimalParameters&, double r, double& e, double& de)
  {
    assert( r > 0. );
    double ratio = p.A / r;
    double rinv = 1. / r;
    e  = ratio * exp( - p.kappa * r );
    de = e * ( rinv - p.kappa );
  }
  
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::YukawaParms;
  using onika::physics::Quantity;

  template<> struct convert<YukawaParms>
  {
    static bool decode(const Node& node, YukawaParms& v)
    {
      v = YukawaParms{};
      if( !node.IsMap() ) { return false; }
      v.A     = node["A"].as<Quantity>().convert();
      v.kappa = node["kappa"].as<Quantity>().convert();
      return true;
    }
  };
}

