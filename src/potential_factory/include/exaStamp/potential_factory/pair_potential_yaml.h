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
#include <exaStamp/potential_factory/pair_potential_factory.h>

#include <onika/physics/units.h>

#include <yaml-cpp/yaml.h>

namespace YAML
{    

  template<> struct convert< exaStamp::PairPotentialUserParameters >
  {
    static inline bool decode(const Node& node, exaStamp::PairPotentialUserParameters& v)
    {
      if( ! node.IsMap() ) { return false; }
      v.m_type_a = node["type_a"].as<std::string>();
      v.m_type_b = node["type_b"].as<std::string>();
      if( node["rcut"] )
      {
        v.m_rcut = node["rcut"].as<onika::physics::Quantity>().convert();
      }
      v.m_pair_potential = exaStamp::PairPotentialFactory::make_instance( node );
      return true;
    }
  };

  template<> struct convert< exaStamp::UserMultiPairPotentials >
  {
    static inline bool decode(const Node& node, exaStamp::UserMultiPairPotentials& v)
    {
      v.m_compiled_potentials.m_rcuts.clear();
      v.m_compiled_potentials.m_force_ops.clear();
      v.m_compiled_potentials.m_pair_id_map.clear();
      
      v.m_compiled_potentials_virial.m_rcuts.clear();
      v.m_compiled_potentials_virial.m_force_ops.clear();
      v.m_compiled_potentials_virial.m_pair_id_map.clear();
     
      return convert< std::vector<exaStamp::PairPotentialUserParameters> >::decode( node , v.m_user_potentials );
      // v.m_user_potentials = node.as< std::vector<PairPotentialUserParameters> >()
      // return true;
    }
  };

}

