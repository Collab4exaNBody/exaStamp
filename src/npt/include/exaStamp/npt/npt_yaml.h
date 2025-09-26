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

#include <exaStamp/npt/npt.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/physics/constants.h>
#include <yaml-cpp/yaml.h>

// bool hasValue(const YAML::Node& node) {
//   return node && node.IsDefined() && !node.IsNull();
// };

namespace YAML
{
  using exaStamp::NPTConfig;
  
  using onika::physics::Quantity;
  
  template<> struct convert<NPTConfig>
  {
    static inline bool decode(const Node& node, NPTConfig& v)
    {
      if(!node.IsMap() ) { return false; }
      v = NPTConfig{};
      if( node["Tstart"] )
        {
          v.m_Tstart = node["Tstart"].as<Quantity>().convert();
        }
      if( node["Tend"] )
        {
          v.m_Tend = node["Tend"].as<Quantity>().convert();
        }      
      if( node["Tdamp"] )
        {
          v.m_Tdamp = node["Tdamp"].as<Quantity>().convert();
        }
      if( node["Pstart"] )
        {
          v.m_Pstart = node["Pstart"].as<double>();
        }
      if( node["Pend"] )
        {
          v.m_Pend = node["Pend"].as<double>();
        }
      if( node["Pdamp"] )
        {
          v.m_Pdamp = node["Pdamp"].as<Quantity>().convert();
        }
      bool modedefined = node["mode"].IsDefined() && !node["mode"].IsNull();
      if( !modedefined )
        {
          std::cout << "A Mode is required either NVT or NPT" << std::endl;
          std::cout << "Please include the following :" << std::endl;
          std::cout << "mode: NVT" << std::endl;
          std::cout << "or" << std::endl;
          std::cout << "mode: NPT" << std::endl;
          std::cout << "to the init_npt YAML block." << std::endl;
          std::abort();
        }
      else
        {
          v.m_mode = node["mode"].as<std::string>();
        }
      if ( node["Pmode"] )
        {
          v.m_Pmode = node["Pmode"].as<std::string>();
        }
      return true;
    }
  };
}

