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

#include <exaStamp/parrinellorahman/parrinellorahman.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/physics/constants.h>
#include <yaml-cpp/yaml.h>

namespace YAML
{
  using exaStamp::ParrinelloRahmanConfig;
  
  using onika::physics::Quantity;
    
  template<> struct convert<ParrinelloRahmanConfig>
  {
    static inline bool decode(const Node& node, ParrinelloRahmanConfig& v)
    {
      if(!node.IsMap() ) { return false; }
      v = ParrinelloRahmanConfig{};
      if( node["Text"] )
      {
        v.m_Text = node["Text"].as<Quantity>().convert();
      }
      if( node["masseNVT"] )
      {
        v.m_masseNVT = node["masseNVT"].as<Quantity>().convert();
      }
      if( node["Pext"] )
      {
        v.m_Pext = node["Pext"].as<Quantity>().convert();
      }
      if( node["masseB"] )
      {
        v.m_masseB = node["masseB"].as<Quantity>().convert();
      }
      if( node["hmask"] )
      {
        v.m_hmask = node["hmask"].as<Mat3d>();
      }
      if( node["hblend"] )
      {
        v.m_hblend = node["hblend"].as<Mat3d>();
      }
      return true;
    }
  };
}

