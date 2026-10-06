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
  
  // LAMMPS pair coul/dsf : erfc(alpha.r) from Abramowitz & Stegun approximation, shifts computed once with exact erfc
  ONIKA_HOST_DEVICE_FUNC
  inline void dsf_compute_energy(const DsfParameters& p, double c, double r, double& e, double& de)
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

  // self energy of a particle with charge q (LAMMPS pair coul/dsf e_self)
  ONIKA_HOST_DEVICE_FUNC
  inline double dsf_self_energy(const DsfParameters& p, double q)
  {
    const double e_self = -( p.e_shift / 2.0 + p.alpha / sqrt(M_PI) ) * q * q * p.qqrd2e;
    return EXASTAMP_QUANTITY( e_self * eV );
  }

  struct DsfKernel
  {
    DsfParameters m_params;
    DsfKernel() = default;
    inline DsfKernel(const DsfParameters& p) : m_params(p) {}
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double c, double r, double& e, double& de) const { dsf_compute_energy( m_params, c, r, e, de ); }
    ONIKA_HOST_DEVICE_FUNC inline double self_energy(double q) const { return dsf_self_energy( m_params, q ); }
  };
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

