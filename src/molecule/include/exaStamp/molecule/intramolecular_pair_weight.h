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
#include <vector>
#include <algorithm>
#include <string>

#include <onika/physics/units.h>
#include <onika/log.h>


namespace exaStamp
{

  struct MolecularPairWeight
  {
    double m_bond_weight = 1.0;
    double m_bend_weight = 1.0;
    double m_torsion_weight = 1.0;
    double m_rf_bond_weight = 1.0;
    double m_rf_bend_weight = 1.0;
    double m_rf_torsion_weight = 1.0;
  };

  struct IntramolecularPairWeighting
  {
    std::map< std::string , MolecularPairWeight > m_molecule_weight;
  };

}

// Yaml conversion operators, allows to read bonds potentials parameters from config file
namespace YAML
{
  using exanb::lerr;
  using exaStamp::MolecularPairWeight;
  using exaStamp::IntramolecularPairWeighting;

  template<> struct convert<IntramolecularPairWeighting>
  {
    static bool decode(const Node& node, IntramolecularPairWeighting& ipw)
    {
      if( !node.IsMap() )
      {
        lerr << "IntramolecularPairWeighting type is not a map as expected" << std::endl;
        return false;
      }

      ipw.m_molecule_weight.clear();

      for(auto p : node)
      {
        std::string molname = p.first.as<std::string>();
        ipw.m_molecule_weight[molname] = MolecularPairWeight{};
        if( p.second["bond"      ] ) ipw.m_molecule_weight[molname].m_bond_weight       = p.second["bond"      ].as<double>();
        if( p.second["bond_rf"   ] ) ipw.m_molecule_weight[molname].m_rf_bond_weight    = p.second["bond_rf"   ].as<double>();

        if( p.second["bend"      ] ) ipw.m_molecule_weight[molname].m_bend_weight       = p.second["bend"      ].as<double>();
        if( p.second["bend_rf"   ] ) ipw.m_molecule_weight[molname].m_rf_bend_weight    = p.second["bend_rf"   ].as<double>();

        if( p.second["torsion"   ] ) ipw.m_molecule_weight[molname].m_torsion_weight    = p.second["torsion"   ].as<double>();
        if( p.second["torsion_rf"] ) ipw.m_molecule_weight[molname].m_rf_torsion_weight = p.second["torsion_rf"].as<double>();
      }

      return true;
    }
  };


}
