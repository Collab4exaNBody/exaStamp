/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <yaml-cpp/yaml.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/log.h>

#include <cassert>

namespace YAML
{

  template<> struct convert< exaStamp::ParticleSpecie >
  {
    static inline Node encode(const exaStamp::ParticleSpecie& atom)
    {
      using namespace exanb;

      Node node;
      node[atom.name()]["mass"] = atom.m_mass;
      node[atom.name()]["z"] = atom.m_z;
      node[atom.name()]["charge"] = atom.m_charge;
      if( ! atom.molecule_name().empty() )
      {
        node[atom.name()]["molecule"] = atom.molecule_name();
      }
      return node;
    }
    
    static inline bool decode(const Node& node, exaStamp::ParticleSpecie& atom)
    {
      using namespace exanb;

      if( ! node.IsMap() ) { return false; }
      if( node.size() != 1 ) { return false; }
      atom.set_name( node.begin()->first.as<std::string>() );

      atom.m_rigid_atom_count = 0;
      //atom.m_rigid_atom_names.clear();
      
      YAML::Node params = node.begin()->second;

      atom.m_mass = params["mass"].as<Quantity>().convert(); 
      if(params["charge"]  ) { atom.m_charge = params["charge"].as<Quantity>().convert(); }
      if(params["z"]       ) { atom.m_z = params["z"].as<int>(); }
      if(params["molecule"]) { atom.set_molecule_name( params["molecule"].as<std::string>() ); }
      if(params["rigid_molecule"])
      {
        if( ! params["rigid_molecule"].IsSequence() )
        {
          lerr << "rigid_molecule description in particle '"<<atom.name()<<"' must be a list"<<std::endl;
          return false;
        }
        if( params["rigid_molecule"].size() < 2 )
        {
          lerr << "rigid_molecule must be described with at least 2 atoms"<< std::endl;
          return false;
        }
        for(auto p: params["rigid_molecule"])
        {
          if( !p.IsMap() || p.size()!=1 )
          {
            lerr << "rigid molecule atom should be described with a map of the form 'name: [ rx, ry, rz ]'" << std::endl;
            return false;
          }
          atom.m_rigid_atoms[ atom.m_rigid_atom_count ] = { p.begin()->second.as<exanb::Vec3d>() , -1 };
          atom.set_rigid_atom_name( atom.m_rigid_atom_count , p.begin()->first.as<std::string>() );
          ++ atom.m_rigid_atom_count ;
          //assert( atom.m_rigid_atom_names.size() == atom.m_rigid_atom_count );
        }
      }
      else // single atom case
      {
        // refers itself as the only rigid molecule atom
        atom.m_rigid_atoms[ atom.m_rigid_atom_count ] = { Vec3d{0.,0.,0.} , -1 };
        atom.set_rigid_atom_name( atom.m_rigid_atom_count , atom.name() );
        atom.m_rigid_atom_count = 1;
        //assert( atom.m_rigid_atom_names.size() == atom.m_rigid_atom_count );
      }
      return true;
    }

  };


  template<> struct convert< exaStamp::ParticleSpecies >
  {
    static inline Node encode(const exaStamp::ParticleSpecies& v)
    {
      Node node;
      for(const auto& x:v) node.push_back(x);
      return node;
    }

    static inline bool decode(const Node& node, exaStamp::ParticleSpecies& v)
    {
      if( ! node.IsSequence() ) return false;
      v.clear();
      for(auto i: node) v.push_back( i.as<exaStamp::ParticleSpecie>() );
      return true;
    }
  };

}

