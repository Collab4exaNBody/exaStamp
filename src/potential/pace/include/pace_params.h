/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>

#include <string>

using namespace exanb;
// PACE Parameters
struct PaceParams
{
  std::string pace_coef;
  bool recursive = true;
  int nt = 1;
  int chunksize = 2048;
};

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<PaceParams>
  {
    static bool decode(const Node& node, PaceParams& v)
    {
      if( !node.IsMap() ) { return false; }
      if( ! node["coef"] ) { return false; }
      if( node["coef"] ) { v.pace_coef = node["coef"].as<std::string>(); }
      if( node["algorithm"] ) {
        v.recursive = ( node["algorithm"].as<std::string>() == "recursive" );
      }
      if( node["nspecies"] ) { v.nt = node["nspecies"].as<int>(); }
      if( node["chunksize"] ) { v.nt = node["chunksize"].as<int>(); }
      return true;
    }
  };
}
