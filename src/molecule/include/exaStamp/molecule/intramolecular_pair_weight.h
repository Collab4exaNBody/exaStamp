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
