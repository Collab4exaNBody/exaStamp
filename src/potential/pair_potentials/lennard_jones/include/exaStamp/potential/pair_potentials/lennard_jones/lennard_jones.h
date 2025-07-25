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

  // LennardJones Parameters
  struct LennardJonesParms
  {
    double epsilon = 0.0;
    double sigma = 0.0;
  };

//# pragma omp declare simd uniform(p,ppp) notinbranch
  ONIKA_HOST_DEVICE_FUNC
  inline void lj_compute_energy(const LennardJonesParms& p, const PairPotentialMinimalParameters&, double r, double& e, double& de)
  {
    assert( r > 0. );
    double ratio = p.sigma / r;
    double ratio6 = pow(ratio,6);    // attractive
    double ratio12 = pow(ratio,12);  // repulsive
    e = 4. * p.epsilon * (ratio12-ratio6) ;
    de = ( -24. * p.epsilon * (2.*ratio12-ratio6) ) / r;
  }
  
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::LennardJonesParms;
  using onika::physics::Quantity;

  template<> struct convert<LennardJonesParms>
  {
    static bool decode(const Node& node, LennardJonesParms& v)
    {
      v = LennardJonesParms{};
      if( !node.IsMap() ) { return false; }
      v.epsilon = node["epsilon"].as<Quantity>().convert();
      v.sigma   = node["sigma"]  .as<Quantity>().convert();
      //std::cout<<"epsilon="<<v.epsilon<<", sigma="<<v.sigma<<std::endl;
      return true;
    }
  };
}

