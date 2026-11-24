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
#include <onika/log.h>
#include <initializer_list>
#include <string>

namespace exaStamp
{

  // Fatal error on a key of a YAML map that the decoder does not read: a misspelled or
  // misplaced parameter (e.g. coul_wolf instead of cw) would otherwise be ignored silently.
  inline void check_yaml_keys(const YAML::Node& node, std::initializer_list<const char*> keys, const char* context)
  {
    if( ! node.IsMap() ) return;
    for(const auto& kv : node)
    {
      const std::string key = kv.first.as<std::string>();
      bool known = false;
      for(const char* k : keys) known = known || key == k;
      if( ! known )
      {
        std::string expected;
        for(const char* k : keys) expected += std::string( expected.empty() ? "" : ", " ) + k;
        onika::fatal_error() << context << ": unknown key '" << key << "' (expected " << expected << ")" << std::endl;
      }
    }
  }

}
