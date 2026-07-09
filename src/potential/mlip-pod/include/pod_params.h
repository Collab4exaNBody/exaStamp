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
#include <string>

struct PodParams
{
  std::string pod_file;    // POD descriptor file (.pod or .txt)
  std::string coeff_file;  // trained coefficient file
};

namespace YAML
{
  template<> struct convert<PodParams>
  {
    static bool decode(const Node& node, PodParams& v)
    {
      if( !node.IsMap() )       { return false; }
      if( !node["pod_file"] )   { return false; }
      if( !node["coeff_file"] ) { return false; }
      v.pod_file   = node["pod_file"].as<std::string>();
      v.coeff_file = node["coeff_file"].as<std::string>();
      return true;
    }
  };
}
