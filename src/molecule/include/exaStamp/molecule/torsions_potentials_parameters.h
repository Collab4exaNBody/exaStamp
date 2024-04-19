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
  struct TorsionPotential
  {
    std::string type;
    std::array<std::string,4> species;

    //std::vector<double> params;
  //  std::function<double(double)> f_energy;
  //  std::function<double(double)> f_forces;
    std::shared_ptr<IntraMolecularPotentialFunctional> m_potential_function;
  };

  struct TorsionsPotentialParameters
  {
    std::vector<TorsionPotential> m_potentials;
    std::unordered_map< uint64_t, std::shared_ptr<IntraMolecularPotentialFunctional> > m_type_to_potential;
  };

}

// Yaml conversion operators, allows to read torsions potentials parameters from config file
namespace YAML
{
  using exanb::lerr;
  using exanb::fatal_error;
  using exaStamp::TorsionPotential;
  using exaStamp::TorsionsPotentialParameters;
  using exaStamp::IntraMolecularPotentialFunctional;
  using exaStamp::IntraMolecularHarmFunctional;
  using exaStamp::IntraMolecularCompassFunctional;
  using exaStamp::IntraMolecularHalfCompassFunctional;
  using exaStamp::IntraMolecularOPLSFunctional;
  using exaStamp::IntraMolecularCosTwoFunctional;

  template<> struct convert<TorsionPotential>
  {
    static bool decode(const Node& node, TorsionPotential& b)
    {
      // potential gestion
      if(!node["potential"])
      {
        fatal_error() << "No potential type given in the torsion interaction. Choice are harm_torsion, torsion, cos_two, half_torsion, opls_torsion, torsion_gaff." << std::endl;
      }
      b.type = node["potential"].as<std::string>();

      // Types gestion
      if( !node["types"] )
      {
        fatal_error() << "No types elements given in the torsion interaction (example types: [H, C, C, H])." << std::endl;
      }
      b.species = {node["types"][0].as<std::string>(), node["types"][1].as<std::string>(), node["types"][2].as<std::string>(), node["types"][3].as<std::string>()};

      if( !node["parameters"] && b.type!="no_potential" )
      {
        fatal_error() << "No parameters elements given in the torsion interaction." << std::endl;
      }

      // Energy and forces gestion
      if(b.type=="harm_torsion")
      {
        //b.params.push_back(node["parameters"]["k"   ].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["phi0"].as<Quantity>().convert());
        double k    = node["parameters"]["k"   ].as<Quantity>().convert(); //b.params.at(0);
        double phi0 = node["parameters"]["phi0"].as<Quantity>().convert(); //b.params.at(1);

        // e = 1/2 k * (phi-phi0)^2
        // F = k * (phi-phi0)^2
        /*
        b.f_energy = [k,phi0](const double phi){return 0.5 * k * pow((phi - phi0),2);};
        b.f_forces = [k,phi0](const double phi){return       k *     (phi - phi0);};
        */
        b.m_potential_function = std::make_shared<IntraMolecularHarmFunctional>(k,phi0);
      }
      else if(b.type=="compass_torsion")
      {
        //b.params.push_back(node["parameters"]["k1"].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["k2"].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["k3"].as<Quantity>().convert());
        double k1 = node["parameters"]["k1"].as<Quantity>().convert(); //b.params.at(0);
        double k2 = node["parameters"]["k2"].as<Quantity>().convert(); //b.params.at(1);
        double k3 = node["parameters"]["k3"].as<Quantity>().convert(); //b.params.at(2);

        // e = k1 * (1-cos(phi)) +     k2 * (1-cos(2(phi)))  +     k3 * (1-cos(3(phi)))
        // F = k1 * sin(phi)     + 2 * k2 * sin(2*(phi))     + 3 * k3 * sin(3(phi))
/*
        b.f_energy = [k1,k2,k3](const double phi){
                       return   k1 * ( 1 - cos(  phi) )
                         +      k2 * ( 1 - cos(2*phi) )
                         +      k3 * ( 1 - cos(3*phi) );};
        b.f_forces = [k1,k2,k3](const double phi){
                       return       k1 * sin(  phi)
                         +      2 * k2 * sin(2*phi)
                         +      3 * k3 * sin(3*phi);};
*/
        b.m_potential_function = std::make_shared<IntraMolecularCompassFunctional>(k1,k2,k3);
      }
      else if(b.type=="0.5compass_torsion")
      {
        //b.params.push_back(node["parameters"]["k1"].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["k2"].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["k3"].as<Quantity>().convert());
        double k1 = node["parameters"]["k1"].as<Quantity>().convert(); //b.params.at(0);
        double k2 = node["parameters"]["k2"].as<Quantity>().convert(); //b.params.at(1);
        double k3 = node["parameters"]["k3"].as<Quantity>().convert(); //b.params.at(2);

        // e = k1/2 * (1-cos(phi)) +     k2/2 * (1-cos(2(phi)))  +     k3/2 * (1-cos(3(phi)))
        // F = k1/2 * sin(phi)     +     k2 * sin(2*(phi))       + 3/2 * k3 * sin(3(phi))
/*
        b.f_energy = [k1,k2,k3](const double phi){
                       return   k1/2 * ( 1 - cos(  phi) )
                         +      k2/2 * ( 1 - cos(2*phi) )
                         +      k3/2 * ( 1 - cos(3*phi) );};
        b.f_forces = [k1,k2,k3](const double phi){
                       return 0.5 * k1 * sin(  phi)
                         +          k2 * sin(2*phi)
                         +    1.5 * k3 * sin(3*phi);};
*/
        b.m_potential_function = std::make_shared<IntraMolecularHalfCompassFunctional>(k1,k2,k3);
      }
      else if(b.type=="opls_torsion")
      {
        //b.params.push_back(node["parameters"]["k1"].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["k2"].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["k3"].as<Quantity>().convert());
        double k1 = node["parameters"]["k1"].as<Quantity>().convert(); //b.params.at(0);
        double k2 = node["parameters"]["k2"].as<Quantity>().convert(); //b.params.at(1);
        double k3 = node["parameters"]["k3"].as<Quantity>().convert(); //b.params.at(2);

        // e =   0.5 * k1 * (1+cos(phi)) + 0.5* k2 * (1-cos(2(phi))) + 0.5 * k3 * (1+cos(3(phi))) 
        // F = - 0.5 * k1 * sin(phi)     +  k2 * sin(2*phi)          - 3/2 * k3 * sin(3*phi)     
/*
        b.f_energy = [k1,k2,k3](const double phi){
                       return   0.5 * k1 * ( 1 + cos(  phi) )
                         +      0.5 * k2 * ( 1 - cos(2*phi) )
                         +      0.5 * k3 * ( 1 + cos(3*phi) );};
        b.f_forces = [k1,k2,k3](const double phi){
                       return   - 0.5  * k1 * sin(  phi)
                         +               k2 * sin(2*phi)
                         -        1.5  * k3 * sin(3*phi);};
*/
        b.m_potential_function = std::make_shared<IntraMolecularOPLSFunctional>(k1,k2,k3);
      }
    else if(b.type=="cos_two")
      {
        //b.params.push_back(node["parameters"]["k"].as<Quantity>().convert());
        //b.params.push_back(node["parameters"]["phi0"].as<Quantity>().convert());
        double k = node["parameters"]["k"].as<Quantity>().convert(); //b.params.at(0);
        double phi0 = node["parameters"]["phi0"].as<Quantity>().convert(); //b.params.at(1);

        // e = k/2 * (1-cos(2(phi-phi0)))
        // F = k * sin(2*(phi-phi0))
/*
        b.f_energy = [k,phi0](const double phi){
                       double tmp = phi-phi0;
                       return   0.5 * k * ( 1 - cos(2*tmp) );};
        b.f_forces = [k,phi0](const double phi){
                       double tmp = phi-phi0;
                       return         k * sin(2*tmp);};
*/
        b.m_potential_function = std::make_shared<IntraMolecularCosTwoFunctional>(k,phi0);
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
        fatal_error() << "Type of torsion potential " << b.type << " unknown. Type can be harm_torsion, compass_torsion, 0.5compass_torsion, opls_torsion, cos_two, no_potential" << std::endl;
      }
      return true;
    }
  };
  
  template<> struct convert<TorsionsPotentialParameters>
  {
    static bool decode(const Node& node, TorsionsPotentialParameters& bs)
    {
      if( !node.IsSequence() )
      {
        fatal_error() << "TorsionPotential type is not readable from config file." << std::endl;
      }
      bs.m_potentials.clear();
      for(size_t n=0;n<node.size();++n)
      {
        bs.m_potentials.push_back(node[n].as<TorsionPotential>());
      }
      return true;
    }
  };
  
  
}
