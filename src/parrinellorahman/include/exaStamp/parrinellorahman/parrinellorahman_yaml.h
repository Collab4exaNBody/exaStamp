#pragma once

#include <exaStamp/parrinellorahman/parrinellorahman.h>
#include <onika/math/basic_types_yaml.h>
#include <exanb/core/physics_constants.h>
#include <yaml-cpp/yaml.h>

namespace YAML
{
  using exaStamp::ParrinelloRahmanConfig;
  using exanb::UnityConverterHelper;
  using exanb::Quantity;
    
  template<> struct convert<ParrinelloRahmanConfig>
  {
    static inline bool decode(const Node& node, ParrinelloRahmanConfig& v)
    {
      if(!node.IsMap() ) { return false; }
      v = ParrinelloRahmanConfig{};
      if( node["Text"] )
      {
        v.m_Text = node["Text"].as<Quantity>().convert();
      }
      if( node["masseNVT"] )
      {
        v.m_masseNVT = node["masseNVT"].as<Quantity>().convert();
      }
      if( node["Pext"] )
      {
        v.m_Pext = node["Pext"].as<Quantity>().convert();
      }
      if( node["masseB"] )
      {
        v.m_masseB = node["masseB"].as<Quantity>().convert();
      }
      if( node["hmask"] )
      {
        v.m_hmask = node["hmask"].as<Mat3d>();
      }
      if( node["hblend"] )
      {
        v.m_hblend = node["hblend"].as<Mat3d>();
      }
      return true;
    }
  };
}

