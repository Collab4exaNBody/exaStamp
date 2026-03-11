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

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

namespace exaStamp
{
  using namespace exanb;

  struct DsfParameters
  {
    double alpha = 0.0;
    double rc = 0.0;
    double qqrd2e = 14.399645;
    double e_shift = 0.0;
    double f_shift = 0.0;

    ONIKA_HOST_DEVICE_FUNC
    inline bool is_null() const { return alpha==0.0 && e_shift==0.0 && f_shift==0.0; }    
    
  };
  
  ONIKA_HOST_DEVICE_FUNC
  inline void dsf_compute_energy(const DsfParameters& p, double c, double r, double& e, double& de)
  {
    assert( r > 0. );
    double MY_PIS = sqrt(M_PI);
    double rsq = r*r;
    double cut_coul = p.rc;
    double cut_coulsq = cut_coul * cut_coul;    

    double erfcc = erfc(p.alpha * cut_coul);
    double erfcd = exp(-p.alpha * p.alpha * cut_coul * cut_coul);
    double f_shift = -(erfcc / cut_coulsq + 2.0 / MY_PIS * p.alpha * erfcd / cut_coul);
    double e_shift = erfcc / cut_coul - f_shift * cut_coul;

    double prefactor = p.qqrd2e * c / r;
    erfcd = exp(-p.alpha * p.alpha * rsq);
    double t = 1.0 / (1.0 + EWALD_P * p.alpha * r);
    erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;

    double forcecoul = prefactor * (erfcc / r + 2.0 * p.alpha / MY_PIS * erfcd + r * f_shift) * r;
    double fpair = -forcecoul / r;
    double ecoul = prefactor * (erfcc - r * e_shift - rsq * f_shift);
    
    e = EXASTAMP_QUANTITY( ecoul * eV );
    de = EXASTAMP_QUANTITY( fpair * eV / ang );
    
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::DsfParameters;
  
  using onika::physics::Quantity;

  template<> struct convert<DsfParameters>
  {
    static bool decode(const Node& node, DsfParameters& v)
    {
      if( !node.IsMap() ) { return false; }
      v.alpha = node["alpha"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();

      double MY_PIS = sqrt(M_PI);
      double cut_coul = v.rc;
      double cut_coulsq = cut_coul * cut_coul;
      double erfcc = erfc(v.alpha * cut_coul);
      double erfcd = exp(-v.alpha * v.alpha * cut_coul * cut_coul);
      v.f_shift = -(erfcc / cut_coulsq + 2.0 / MY_PIS * v.alpha * erfcd / cut_coul);
      v.e_shift = erfcc / cut_coul - v.f_shift * cut_coul;
        
      return true;
    }
  };
}

