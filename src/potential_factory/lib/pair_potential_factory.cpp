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

#include <exaStamp/potential_factory/pair_potential_factory.h>

#include <exanb/core/grid.h>
#include <onika/log.h>
#include <onika/plugin.h>
#include <onika/yaml/yaml_utils.h>

#include <cassert>
#include <iostream>
#include <memory>

namespace exaStamp
{
  using namespace exanb;

  using PairPotentialCreator = std::function<std::shared_ptr<PairPotential>(YAML::Node)>;

  bool PairPotentialFactory::register_factory( const std::string& name, PairPotentialCreator creator )
  {
    if( ! onika::quiet_plugin_register() ) { lout<<"  potential   "<<name<<std::endl; }
    onika::plugin_db_register( "potential" , name );
    s_factories[ name ] = creator;
    return true;
  }

  std::shared_ptr<PairPotential> PairPotentialFactory::make_instance( YAML::Node node )
  {
    if( ! node["potential"] )
    {
      lerr << "Fatal: Yaml node does not have a 'potential' entry." << std::endl;
      std::abort();
    }

    std::string pot_name = node["potential"].as<std::string>();
    onika::check_load_plugins_for( "potential" , pot_name );

    auto it = s_factories.find( pot_name );
    if( it == s_factories.end() )
    {
      fatal_error()<<"Could not find a potential factory for name '"<<pot_name<<"'"<<std::endl;
    }

    YAML::Node tmp;
    if( node["parameters"] )
    {
      tmp = node["parameters"];
    }
    else
    {
      tmp = "null";
    }
    return it->second ( tmp );
  }

  std::list<std::string> PairPotentialFactory::available_potentials()
  {
    std::list<std::string> result;
    for(const auto& f : s_factories)
    {
	    result.push_back( f.first );
    }
    return result;
  }

  std::map< std::string, PairPotentialCreator > PairPotentialFactory::s_factories;
}

