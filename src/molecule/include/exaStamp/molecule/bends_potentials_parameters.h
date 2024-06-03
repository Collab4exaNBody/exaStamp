#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <vector>
#include <algorithm>
#include <string>

#include <exanb/core/quantity_yaml.h>
#include <exaStamp/particle_species/particle_specie.h>
//#include "exanb/container_utils.h"
#include <exanb/core/log.h>

#include <exaStamp/molecule/potential_functional.h>

namespace exaStamp
{
  struct BendPotential
  {
    std::string type;
    std::array<std::string,3> species;
    std::shared_ptr<IntraMolecularPotentialFunctional> m_potential_function;
  };

  struct BendsPotentialParameters
  {
    std::vector<BendPotential> m_potentials;
    std::unordered_map< uint64_t, std::shared_ptr<IntraMolecularPotentialFunctional> > m_type_to_potential;
  };

  //using BendsPotentialParameters = std::vector<BendPotential>;
}

// Yaml conversion operators, allows to read bends potentials parameters from config file
namespace YAML
{
  using exanb::lerr;
  using exaStamp::BendsPotentialParameters;
  using exaStamp::BendPotential;
  using exaStamp::IntraMolecularPotentialFunctional;
  using exaStamp::IntraMolecularHarmFunctional;
  using exaStamp::IntraMolecularBondOPLSFunctional;
  using exaStamp::IntraMolecularQuarFunctional;

  template<> struct convert<BendsPotentialParameters>
  {
    static bool decode(const Node& node, BendsPotentialParameters& bs)
    {
      if( !node.IsSequence() )
        {
          lerr << "BendPotential type is not readable from config file." << std::endl;
          abort();
        }

      for(size_t n=0;n<node.size();++n)
        {
          bs.m_potentials.push_back(node[n].as<BendPotential>());
        }

      return true;
    }
  };

  template<> struct convert<BendPotential>
  {
    static bool decode(const Node& node, BendPotential& b)
    {
      // potential gestion
      if(!node["potential"])
        {
          lerr << "No potential type given in the bond interaction. Choice are harm_bend, bend and no_potential." << std::endl;
          abort();
        }
      b.type = node["potential"].as<std::string>();

      // Types gestion
      if( !node["types"] )
        {
          lerr << "No types elements given in the bend interaction (example types: [C, O, C])." << std::endl;
          abort();
        }
      b.species = {node["types"][0].as<std::string>(), node["types"][1].as<std::string>(), node["types"][2].as<std::string>()};

      if( !node["parameters"] && b.type!="no_potential" )
        {
          lerr << "No parameters elements given in the bend interaction." << std::endl;
          abort();
        }

      // Energy and forces gestion
      if(b.type=="harm_bend")
        {
          double k  = node["parameters"]["k" ].as<Quantity>().convert();
          double theta0 = node["parameters"]["theta0"].as<Quantity>().convert();
          b.m_potential_function = std::make_shared<IntraMolecularHarmFunctional>(k,theta0);
        }
      if(b.type=="opls_bend")
        {
          double k  = node["parameters"]["k" ].as<Quantity>().convert();
          double theta0 = node["parameters"]["theta0"].as<Quantity>().convert();
          b.m_potential_function = std::make_shared<IntraMolecularBondOPLSFunctional>(k,theta0);
        }
      else if(b.type=="quar_bend")
        {
          double k2 = node["parameters"]["k2"].as<Quantity>().convert();
          double k3 = node["parameters"]["k3"].as<Quantity>().convert();
          double k4 = node["parameters"]["k4"].as<Quantity>().convert();
          double theta0 = node["parameters"]["theta0"].as<Quantity>().convert();
          b.m_potential_function = std::make_shared<IntraMolecularQuarFunctional>(k2,k3,k4,theta0);
        }
      else if(b.type=="no_potential")
        {
          b.m_potential_function = std::make_shared<IntraMolecularPotentialFunctional>();
        }
      else
        {
          lerr << "Type of bend potential " << b.type << " unknown. Type can be harm_bend, quar_bend, no_potential" << std::endl;
          abort();
        }
      return true;
    }
  };
}
