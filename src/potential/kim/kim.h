#pragma once

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>

#include <string>

using namespace exanb;

// NNP Parameters
struct KIMParams
{
  std::string model;
  double rcut;  
};

// Yaml conversion operators, allows to read NNP parameters from config file for n2p2
namespace YAML
{
  template<> struct convert<KIMParams>
  {
    static inline bool decode(const Node& node, KIMParams& v)
    {
      if( !node.IsMap()       ) { return false; }
      if( ! node["model"      ] ) { return false; }
      v.model = node["model"      ].as<std::string>();
      if( ! node["rcut"      ] ) { return false; }
      v.rcut = node["rcut"      ].as<double>();
      return true;
    }
  };
}


