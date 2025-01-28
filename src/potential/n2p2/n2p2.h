#pragma once

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>

#include <string>

using namespace exanb;

// NNP Parameters
struct NNPParams
{
  std::string dir;
  bool        showew;
  bool        resetew;
  int         showewsum;
  int         maxew;
  double      cflength;
  double      cfenergy;
  double      cutoff;
  bool		dumpout; // To enable/disable i/o on stdout
};

#define NNP_INTERFACE_VERSION     "v0.0.1"
#define NNP_INTERFACE_NAME        "Interface ExaStamp-n2p2 : a Neural Network Potential Package"
#define NNP_FLAG_DEBUG            "### [DEBUG - ExaStamp/n2p2]"
#define NNP_FLAG_INFO             "[ExaStamp/n2p2]"

// Yaml conversion operators, allows to read NNP parameters from config file for n2p2
namespace YAML
{
  template<> struct convert<NNPParams>
  {
    static inline bool decode(const Node& node, NNPParams& v)
    {
      if( !node.IsMap()       ) { return false; }
      
      if( ! node["dir"      ] ) { return false; }
      if( ! node["showew"   ] ) { return false; }
      if( ! node["resetew"  ] ) { return false; }
      if( ! node["showewsum"] ) { return false; }
      if( ! node["maxew"    ] ) { return false; }
      if( ! node["cflength" ] ) { return false; }
      if( ! node["cfenergy" ] ) { return false; }
      if( ! node["cutoff"   ] ) { return false; }
      if( ! node["dump_out" ] ) { return false; }

      static const double conv_energy_inv =  1e-4 * onika::physics::elementaryCharge / onika::physics::atomicMass;
      
      v.dir       = node["dir"      ].as<std::string>();
      v.showew    = node["showew"   ].as<bool>();
      v.resetew   = node["resetew"  ].as<bool>();
      v.showewsum = node["showewsum"].as<int>();
      v.maxew     = node["maxew"    ].as<int>();
      v.cflength  = node["cflength" ].as<double>();
      v.cfenergy  = node["cfenergy" ].as<double>() / conv_energy_inv; 
      v.cutoff    = node["cutoff"   ].as<double>();
      v.dumpout   = node["dump_out" ].as<bool>();

      return true;
    }
  };
}


