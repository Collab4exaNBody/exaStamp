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
  struct ImproperPotential
  {
    std::string type;
    std::array<std::string,4> species;
  
    //std::vector<double> params;
    //std::function<double(double)> f_energy;
    //std::function<double(double)> f_forces;
    std::shared_ptr<IntraMolecularPotentialFunctional> m_potential_function;
  };
  
  struct ImpropersPotentialParameters
  {
    std::vector<ImproperPotential> m_potentials;
    std::unordered_map< uint64_t, std::shared_ptr<IntraMolecularPotentialFunctional> > m_type_to_potential;
  };
}

// Yaml conversion operators, allows to read torsions potentials parameters from config file
namespace YAML
{
  using exanb::lerr;
  using exaStamp::ImproperPotential;
  using exaStamp::ImpropersPotentialParameters;
  using exaStamp::IntraMolecularPotentialFunctional;
  using exaStamp::IntraMolecularHarmFunctional;
  using exaStamp::IntraMolecularCosOPLSFunctional;

  template<> struct convert<ImpropersPotentialParameters>
  {
    static bool decode(const Node& node, ImpropersPotentialParameters& bs)
    {
      if( !node.IsSequence() )
      {
        lerr << "ImproperPotential type is not readable from config file." << std::endl;
        abort();
      }
      bs.m_potentials.clear();
      for(size_t n=0;n<node.size();++n)
      {
        bs.m_potentials.push_back(node[n].as<ImproperPotential>());
      }
      return true;
    }
  };

  template<> struct convert<ImproperPotential>
  {
    static bool decode(const Node& node, ImproperPotential& b)
    {
      // potential gestion
      if(!node["potential"])
        {
          lerr << "No potential type given in the improper interaction. Choice are harm_improper, no_potential." << std::endl;
          abort();
        }
      b.type = node["potential"].as<std::string>();

      // Types gestion
      if( !node["types"] )
        {
          lerr << "No types elements given in the improper interaction (example types: [C, C, C, H])." << std::endl;
          abort();
        }
      b.species = {node["types"][0].as<std::string>(), node["types"][1].as<std::string>(), node["types"][2].as<std::string>(), node["types"][3].as<std::string>()};

      if( !node["parameters"]  && b.type!="no_potential")
        {
          lerr << "No parameters elements given in the improper interaction." << std::endl;
          abort();
        }

      // Energy and forces gestion
      if(b.type=="harm_improper")
        {

          //b.params.push_back(node["parameters"]["k"   ].as<Quantity>().convert());
          //b.params.push_back(node["parameters"]["chi0"].as<Quantity>().convert());
          double k    = node["parameters"]["k"   ].as<Quantity>().convert(); //b.params.at(0);
          double chi0 = node["parameters"]["chi0"].as<Quantity>().convert(); //b.params.at(1);


          // e = 1/2 k * (chi-chi0)^2
          // F = k * (chi-chi0)
          /*
          b.f_energy = [k,chi0](const double chi){return 0.5 * k * pow((chi - chi0),2);};
          b.f_forces = [k,chi0](const double chi){return       k *     (chi - chi0);};
          */
          b.m_potential_function = std::make_shared<IntraMolecularHarmFunctional>(k,chi0);
        }
      else if(b.type=="opls_improper")
        {

          //b.params.push_back(node["parameters"]["k1"   ].as<Quantity>().convert());
          //b.params.push_back(node["parameters"]["k2"   ].as<Quantity>().convert());
          //b.params.push_back(node["parameters"]["k3"   ].as<Quantity>().convert());
          double k1    = node["parameters"]["k1"].as<Quantity>().convert(); //b.params.at(0);
          double k2    = node["parameters"]["k2"].as<Quantity>().convert(); //b.params.at(1);
          double k3    = node["parameters"]["k3"].as<Quantity>().convert(); //b.params.at(2);


          // e = k1/2(1 + cos(chi)) + k2/2(1-cos(2 chi)) + k3/2(1+cos(3 chi))
          // F = -k1/2 sin(chi) + k2 sin(2 chi) - 3/2 k3 sin(3 chi)
          //b.f_energy = [k1,k2,k3](const double chi){return   0.5 * k1 * (1 + cos(  chi))
          //                                          +      0.5 * k2 * (1 - cos(2*chi))
          //                                          +      0.5 * k3 * (1 + cos(3*chi));};
          //b.f_forces = [k1,k2,k3](const double chi){return  -0.5 * k1 * sin(  chi)
          //                                          +            k2 * sin(2*chi)  
          //                                          -      1.5 * k3 * sin(3*chi);};

          b.m_potential_function = std::make_shared<IntraMolecularCosOPLSFunctional>(k1,k2,k3);
        }
      else if(b.type=="no_potential")
        {
          // e = 0
          // F = 0
          b.m_potential_function = std::make_shared<IntraMolecularPotentialFunctional>();
          //b.f_energy = [](const double){return 0;};
          //b.f_forces = [](const double){return 0;};
        }
      else
        {
          lerr << "Type of improper potential " << b.type << " unknown. Type can be harm_improper, harm_improper1.0, opls_improper, no_potential" << std::endl;
          abort();
        }
      return true;
    }
  };
}
