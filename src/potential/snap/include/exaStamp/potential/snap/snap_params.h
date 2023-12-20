#pragma once

#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>

#include <string>

namespace exaStamp
{
  using namespace exanb;


  // Snap Parameters
  struct SnapParms
  {
    std::string lammps_param;
    std::string lammps_coef;
    int nt = 2;
  };

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::SnapParms>
  {
    static bool decode(const Node& node, exaStamp::SnapParms& v)
    {
      if( !node.IsMap() ) { return false; }
      if( ! node["param"] ) { return false; }
      if( ! node["coef"] ) { return false; }
      v.lammps_param = node["param"].as<std::string>();
      v.lammps_coef  = node["coef"].as<std::string>();
      if( node["nt"] ) { v.nt = node["nt"].as<int>(); }
      return true;
    }
  };
}


