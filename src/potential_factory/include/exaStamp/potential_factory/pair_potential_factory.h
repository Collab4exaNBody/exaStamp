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

#include <exaStamp/potential_factory/pair_potential.h>

#include <string>
#include <memory>
#include <functional>
#include <map>
#include <list>
#include <yaml-cpp/yaml.h>

namespace exaStamp
{
  using namespace exanb;

  struct PairPotentialFactory
  {
    using Instance = std::shared_ptr<PairPotential>;
    static bool register_factory( const std::string& name, std::function<Instance(YAML::Node)> pot_factory );
    static Instance make_instance( YAML::Node );
    static std::list<std::string> available_potentials();
    private:
    static std::map< std::string, std::function<Instance(YAML::Node)> > s_factories;
  };

}

