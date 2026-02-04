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

#include <exaStamp/particle_species/particle_specie_yaml.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <exanb/mpi/cost_weight_map.h>
#include <onika/physics/units.h>

namespace exaStamp
{
    using UserTypeCostMap = std::map<std::string,double>;
}

namespace YAML
{

  template<> struct convert< exaStamp::UserTypeCostMap >
  {
    static inline bool decode(const Node& node, exaStamp::UserTypeCostMap & cmap)
    {
      if( ! node.IsMap() )
      {
        exanb::fatal_error() << "type cost must be a map" << std::endl;
        return false;
      }
      cmap.clear();
      for(auto p: node)
      {
        cmap[ p.first.as<std::string>() ] = p.second.as<onika::physics::Quantity>().convert();
      }
      return true;
    }
  };
}

namespace exaStamp
{
  using namespace exanb;

  class TypeCostMap : public OperatorNode
  {  
    ADD_SLOT( ParticleSpecies , species , INPUT , REQUIRED );
    ADD_SLOT( UserTypeCostMap , type_costs , INPUT , REQUIRED );
    ADD_SLOT( CostWeightMap , cost_weight_map , INPUT_OUTPUT );

  public:
    inline void execute () override final
    {
      std::unordered_map<std::string,unsigned int> type_map;
      for(unsigned int i=0;i<species->size();i++)
      {
        type_map[species->at(i).name()] = i;
        (*cost_weight_map) [ i ] = 1.0;
      }
      for(const auto& it: *type_costs)
      {
        if( type_map.find(it.first) == type_map.end() )
        {
          fatal_error() << "type '"<<it.first<<"' doesn't exists"<<std::endl;
        }
        (*cost_weight_map) [ type_map[it.first] ] = it.second;
      }
    }
    
    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsMap() && !node["type_costs"] )
      {
        tmp["type_costs"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize(tmp);
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(type_cost_map)
  {
    OperatorNodeFactory::instance()->register_factory( "type_costs", make_simple_operator<TypeCostMap> );
  }

}

