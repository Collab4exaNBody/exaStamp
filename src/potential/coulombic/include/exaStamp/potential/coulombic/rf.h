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

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/cuda/cuda.h>
#include <exaStamp/unit_system.h>

namespace exaStamp
{
  using namespace exanb;

  struct ReactionFieldParameters
  {
    double rc = 12.5; 
    double RF0 = 0.0;
    double RF1 = 0.0;
    double RF2 = 0.0;
    double ecut = 0.0;
    
    ONIKA_HOST_DEVICE_FUNC
    inline bool is_null() const { return RF0==0.0 && RF1==0.0 && RF2==0.0 && ecut==0.0; }    
  };

  ONIKA_HOST_DEVICE_FUNC
  inline void rf_compute_energy(const ReactionFieldParameters& p_rc, double c, double r, double& e, double& de)
  {
    assert( r > 0. );
    
    const double r2 = r * r;
    
    e  = ( ( p_rc.RF0/r - p_rc.RF2 + p_rc.RF1*r2 ) - p_rc.ecut ) * c;
    de = ( -p_rc.RF0/r2 + 2.*p_rc.RF1 * r ) * c;
  }

  inline void init_rf(ReactionFieldParameters& v, double rc, double epsilon)
  {
    static constexpr double Epsilon0 = EXASTAMP_CONST_QUANTITY( onika::physics::epsilonZero * ( C^2 ) * ( s^2 ) / ( m^3 ) / ( kg^1 ) );
    static constexpr double one_FourPiEpsilon0 = 1. / 4. / M_PI / Epsilon0;
    v.rc = rc;
    if( rc==0.0 || epsilon==0.0 )
    {
      v.rc = 1e-18;
      v.RF0 = 0.0;
      v.RF1 = 0.0;
      v.RF2 = 0.0;
      v.ecut = 0.0;
    }
    else
    {
      v.RF0     = one_FourPiEpsilon0;
      v.RF1     = one_FourPiEpsilon0 / pow(v.rc,3) * (epsilon - 1.)/(2.*epsilon + 1.);
      v.RF2     = one_FourPiEpsilon0 / v.rc * (3. * epsilon)/(2.*epsilon + 1.);
      v.ecut    = 0.0;
      double e=0.0, de=0.0;
      rf_compute_energy( v , 1.0 , v.rc , e , de );
      v.ecut = e;
    }      
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::ReactionFieldParameters;
  
  using onika::physics::Quantity;

  template<> struct convert<ReactionFieldParameters>
  {
    static bool decode(const Node& node, ReactionFieldParameters& v)
    {
      if( !node.IsMap() ) { return false; }
      init_rf( v, node["rc"].as<Quantity>().convert() , node["epsilon"].as<double>() );
      return true;
    }
  };
}

