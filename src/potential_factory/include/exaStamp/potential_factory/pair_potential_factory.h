#pragma once

#include <exaStamp/potential_factory/pair_potential.h>

#include <string>
#include <memory>
#include <functional>
#include <map>
#include <list>
#include <yaml-cpp/yaml.h>

namespace exaStamp
{
  using namespace exanb;

  struct PairPotentialFactory
  {
    using Instance = std::shared_ptr<PairPotential>;
    static bool register_factory( const std::string& name, std::function<Instance(YAML::Node)> pot_factory );
    static Instance make_instance( YAML::Node );
    static std::list<std::string> available_potentials();
    private:
    static std::map< std::string, std::function<Instance(YAML::Node)> > s_factories;
  };

}

