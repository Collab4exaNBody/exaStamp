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
#include "KIM_SimulatorHeaders.hpp"
#include "KIM_SupportedExtensions.hpp"

#include <string>

using namespace exanb;

// NNP Parameters
struct KIMParams
{
  std::string model;
  double rcut;  
};

struct KIMThreadContext
{
  KIM::Model * kim_model = nullptr;
};

struct KIMContext
{
  KIM::Model * kim_model = nullptr;
  std::vector<KIMThreadContext> m_thread_ctx;
  //  std::vector<std::shared_ptr<KIM::Model*>> m_test;
  //  std::vector<KIM::Model> m_test;  
};

// Yaml conversion operators, allows to read NNP parameters from config file for n2p2
namespace YAML
{
  template<> struct convert<KIMParams>
  {
    static inline bool decode(const Node& node, KIMParams& v)
    {
      if( !node.IsMap()   ) { return false; }
      if( ! node["model"] ) { return false; }
      if(   node["model"] ) { v.model = node["model"].as<std::string>(); }
      if( ! node["rcut"]  ) { return false; }
      if(   node["rcut"]  ) { v.rcut = node["rcut"].as<double>(); }
      return true;
    }
  };
}


