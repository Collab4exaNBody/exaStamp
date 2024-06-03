#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <functional>
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
  struct BondPotential
  {
    std::string type;
    std::array<std::string,2> species;
    std::shared_ptr<IntraMolecularPotentialFunctional> potential;
  };


  struct BondsPotentialParameters
  {
    std::vector<BondPotential> m_bond_desc;
    std::unordered_map< uint64_t, std::shared_ptr<IntraMolecularPotentialFunctional> > m_type_pair_to_potential;
  };

}

// Yaml conversion operators, allows to read bonds potentials parameters from config file
namespace YAML
{
  using exanb::lerr;
  using exaStamp::BondsPotentialParameters;
  using exaStamp::BondPotential;
  using exaStamp::IntraMolecularPotentialFunctional;
  using exaStamp::IntraMolecularHarmFunctional;
  using exaStamp::IntraMolecularBondOPLSFunctional;
  using exaStamp::IntraMolecularQuarFunctional;

  template<> struct convert<BondsPotentialParameters>
  {
    static bool decode(const Node& node, BondsPotentialParameters& bs)
    {
      if( !node.IsSequence() )
        {
          lerr << "BondPotentialParameters type is not readable from config file." << std::endl;
          abort();
        }

      for(size_t n=0;n<node.size();++n)
        {
          bs.m_bond_desc.push_back(node[n].as<BondPotential>());
        }

      return true;
    }
  };

  template<> struct convert<BondPotential>
  {
    static bool decode(const Node& node, BondPotential& b)
    {
      if(!node["potential"])
        {
          lerr << "No potential type given in the bond interaction. Choice are harmBond and bond." << std::endl;
          abort();
        }

      // potential gestion
      b.type = node["potential"].as<std::string>();


      // Types gestion
      if( !node["types"] )
        {
          lerr << "No types elements given in the bond interaction (example types: [C, O])." << std::endl;
          abort();
        }
      b.species = {node["types"][0].as<std::string>(), node["types"][1].as<std::string>()};

      if( !node["parameters"] && b.type!="no_potential" )
        {
          lerr << "No parameters elements given in the bond interaction." << std::endl;
          abort();
        }

      // Energy and forces gestion
      if(b.type=="harm_bond")
        {
          double k  = node["parameters"]["k" ].as<Quantity>().convert();
          double r0 = node["parameters"]["r0"].as<Quantity>().convert();
          // e = 1/2 k * (r-r0)^2
          // F = k * (r-r0)
          b.potential = std::make_shared<IntraMolecularHarmFunctional>(k,r0);
        }
      if(b.type=="opls_bond") // same as harm_bond without 1/2 factor for energy
        {
          double k  = node["parameters"]["k" ].as<Quantity>().convert();
          double r0 = node["parameters"]["r0"].as<Quantity>().convert();
          // e = k * (r-r0)^2
          // F = k * (r-r0)
          b.potential = std::make_shared<IntraMolecularBondOPLSFunctional>(k,r0);
        }
      else if(b.type=="quar_bond")
        {
          double k2 = node["parameters"]["k2"].as<Quantity>().convert();
          double k3 = node["parameters"]["k3"].as<Quantity>().convert();
          double k4 = node["parameters"]["k4"].as<Quantity>().convert();
          double r0 = node["parameters"]["r0"].as<Quantity>().convert();

          // e = k2 * (r-r0)^2 + k3 * (r-r0)^3 + k4 * (r-r0)^4
          // F = 2*k2*(r-r0) + 3*k3*(r-r0)^2 + 4*k4*(r-r0)^3
          b.potential = std::make_shared<IntraMolecularQuarFunctional>(k2,k3,k4,r0);
        }
      else if(b.type=="no_potential")
        {
          // e = 0
          // F = 0
          b.potential = std::make_shared<IntraMolecularPotentialFunctional>();
        }
      else
        {
          lerr << "type of potential " << b.type << " unknown. Type can be harm_bond, quar_bond, no_potential" << std::endl;
          abort();
        }

      return true;
    }
  };
}
