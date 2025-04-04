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

