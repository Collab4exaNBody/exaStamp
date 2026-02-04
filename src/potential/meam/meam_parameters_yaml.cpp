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

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include "meam_parameters_yaml.h"

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  bool convert<MeamParameters>::decode(const Node& node, MeamParameters& v)
  {
    using onika::physics::Quantity;
    if( !node.IsMap() ) { return false; }

#   define MEAM_PARAMETER(x) do{ if(!node[#x]){return false;} v.x = node[#x].as<Quantity>().convert(); }while(0)
    MEAM_PARAMETER(rmax);
    MEAM_PARAMETER(rmin);
    MEAM_PARAMETER(Ecoh);
    MEAM_PARAMETER(E0);
    MEAM_PARAMETER(A);
    MEAM_PARAMETER(r0);
    MEAM_PARAMETER(alpha);
    MEAM_PARAMETER(delta);
    MEAM_PARAMETER(beta0);
    MEAM_PARAMETER(beta1);
    MEAM_PARAMETER(beta2);
    MEAM_PARAMETER(beta3);
    MEAM_PARAMETER(t0);
    MEAM_PARAMETER(t1);
    MEAM_PARAMETER(t2);
    MEAM_PARAMETER(t3);
    MEAM_PARAMETER(s0);
    MEAM_PARAMETER(s1);
    MEAM_PARAMETER(s2);
    MEAM_PARAMETER(s3);
    MEAM_PARAMETER(Cmin);
    MEAM_PARAMETER(Cmax);
    MEAM_PARAMETER(Z);
    MEAM_PARAMETER(rc);
    MEAM_PARAMETER(rp);
#   undef MEAM_PARAMETER

    // Why ??? because exaStamp v1 does the same, so we put it there for backward compatibility
    if( v.r0 != 0.0 ) {  v.r0 = 1.0 / v.r0; }
    
    return true;
  }

}

