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

#include <map>
#include <string>

//#include <onika/memory/allocator.h>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>

namespace exaStamp
{
  using DistanceMap = std::map<std::string,double>;
}

namespace YAML
{

  template<> struct convert< exaStamp::DistanceMap >
  {
    static inline Node encode(const exaStamp::DistanceMap& dmap)
    {
      using namespace exanb;
      Node node;
      for(const auto& p:dmap)
      {
        node[ p.first ] = p.second;
      }
      return node;
    }
    
    static inline bool decode(const Node& node, exaStamp::DistanceMap& dmap)
    {
      using onika::physics::Quantity;
      dmap.clear();
      if( ! node.IsMap() ) return false;
      //std::cout << "convert DistanceMap" << std::endl;
      for(auto p: node)
      {
        //std::cout << p.first.as<std::string>() << " -> " << p.second.as<std::string>() << std::endl;
        dmap[ p.first.as<std::string>() ] = p.second.as<Quantity>().convert();
      }
      return true;
    }
    
  };
  
}
