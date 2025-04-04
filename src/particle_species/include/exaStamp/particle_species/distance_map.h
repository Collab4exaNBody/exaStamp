#pragma once

#include <map>
#include <string>

//#include <onika/memory/allocator.h>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>

namespace exaStamp
{
  using DistanceMap = std::map<std::string,double>;
}

namespace YAML
{

  template<> struct convert< exaStamp::DistanceMap >
  {
    static inline Node encode(const exaStamp::DistanceMap& dmap)
    {
      using namespace exanb;
      Node node;
      for(const auto& p:dmap)
      {
        node[ p.first ] = p.second;
      }
      return node;
    }
    
    static inline bool decode(const Node& node, exaStamp::DistanceMap& dmap)
    {
      using onika::physics::Quantity;
      dmap.clear();
      if( ! node.IsMap() ) return false;
      //std::cout << "convert DistanceMap" << std::endl;
      for(auto p: node)
      {
        //std::cout << p.first.as<std::string>() << " -> " << p.second.as<std::string>() << std::endl;
        dmap[ p.first.as<std::string>() ] = p.second.as<Quantity>().convert();
      }
      return true;
    }
    
  };
  
}
