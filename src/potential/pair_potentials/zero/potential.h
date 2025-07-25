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

struct ZeroPotentialParameters
{
};

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::ZeroPotentialParameters>
  {
    static bool decode(const Node& node, exaStamp::ZeroPotentialParameters& v)
    {
      return true;
    }
  };
}

namespace exaStamp
{

ONIKA_HOST_DEVICE_FUNC inline void zero_potential_compute_force(const ZeroPotentialParameters&, const PairPotentialMinimalParameters&, double r, double& e, double& de)
{
  assert( r > 0. );
  e = 0.0 ;
  de = 0.0 ;
}

}

#define USTAMP_POTENTIAL_NAME     zero
#define USTAMP_POTENTIAL_PARAMS   ZeroPotentialParameters
#define USTAMP_POTENTIAL_COMPUTE  zero_potential_compute_force

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

#define USTAMP_POTENTIAL_ENABLE_CUDA 1

