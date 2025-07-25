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

#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>
#include <iomanip>

#include <onika/math/basic_types_def.h>
#include <onika/math/basic_types_operators.h>
//#include <exaStamp/atom_bond_connectivity.h>

#include <onika/memory/allocator.h>
#include <exanb/core/particle_type_id.h>

#ifndef XSTAMP_MAX_RIGID_MOLECULE_ATOMS
#define XSTAMP_MAX_RIGID_MOLECULE_ATOMS 4
#endif

namespace exaStamp
{
  static constexpr size_t MAX_PARTICLE_SPECIES = exanb::MAX_PARTICLE_TYPES;
  static constexpr size_t MAX_ATOM_NAME_LEN = 16;
  static constexpr size_t LEGACY_MAX_RIGID_MOLECULE_ATOMS = 16; // WARNING, changing this constant will make the dump reader fail
  static constexpr size_t MAX_RIGID_MOLECULE_ATOMS = XSTAMP_MAX_RIGID_MOLECULE_ATOMS;

  struct alignas(8) RigidMoleculeAtom
  {
    exanb::Vec3d m_pos {0.,0.,0.}; // relative position
    int m_atom_type = -1;
  };

  template<size_t _MaxRigidMolAtoms = MAX_RIGID_MOLECULE_ATOMS>
  struct ParticleSpecieTmpl
  {
    static inline constexpr size_t MaxRigidMolAtoms = _MaxRigidMolAtoms;
    static inline constexpr size_t MAX_STR_LEN = MAX_ATOM_NAME_LEN;
  
    double m_mass = 0.;
    double m_charge = 0.;
    unsigned int m_z = 0;    
    unsigned int m_rigid_atom_count = 0;
/*
    unsigned int m_molecule_type = 0;     // type of molecule this Atom belongs to
    unsigned int m_molecule_place = 0;    // place in molecule of this atom
*/
    exanb::Vec3d m_minert {0.,0.,0.};
    RigidMoleculeAtom m_rigid_atoms[MaxRigidMolAtoms];

    char m_rigid_atom_names[MaxRigidMolAtoms][MAX_STR_LEN];    
    char m_name[MAX_STR_LEN] = {'\0'};
    char m_molecule_name[MAX_STR_LEN] = {'\0'};
    
    inline void set_name(const std::string& s)
    {
      if( s.length() >= MAX_STR_LEN ) { std::cerr<<"Atom name too long : length="<<s.length()<<", max="<<(MAX_STR_LEN-1)<<"\n"; std::abort(); }
      std::strncpy(m_name,s.c_str(),MAX_STR_LEN); m_name[MAX_STR_LEN-1]='\0';
    }
    inline void set_molecule_name(const std::string& s)
    {
      if( s.length() >= MAX_STR_LEN ) { std::cerr<<"Molecule name too long : length="<<s.length()<<", max="<<(MAX_STR_LEN-1)<<"\n"; std::abort(); }
      std::strncpy(m_molecule_name,s.c_str(),MAX_STR_LEN); m_molecule_name[MAX_STR_LEN-1]='\0';
    }
    inline void set_rigid_atom_name(size_t i, const std::string& s)
    {
      if( s.length() >= MAX_STR_LEN ) { std::cerr<<"Rigid atom name too long : length="<<s.length()<<", max="<<(MAX_STR_LEN-1)<<"\n"; std::abort(); }
      std::strncpy(m_rigid_atom_names[i],s.c_str(),MAX_STR_LEN); m_rigid_atom_names[i][MAX_STR_LEN-1]='\0';
    }
    
    inline std::string name() const { return m_name; }
    inline std::string molecule_name() const { return m_molecule_name; }
    inline std::string rigid_atom_name(size_t i) const { return m_rigid_atom_names[i]; }
    
    template< size_t MaxRA >
    inline void copy_from( const ParticleSpecieTmpl<MaxRA>& from )
    {
      m_mass = from.m_mass;
      m_charge = from.m_charge;
      m_z = from.m_z;
      m_rigid_atom_count = from.m_rigid_atom_count;
      m_minert = from.m_minert;
      set_name( from.name() );
      set_molecule_name( from.molecule_name() );
      for(unsigned int i=0;i<m_rigid_atom_count;i++)
      {
        m_rigid_atoms[i] = from.m_rigid_atoms[i];
        set_rigid_atom_name(i, from.rigid_atom_name(i) );
      }
    }
  };

  using ParticleSpecie = ParticleSpecieTmpl<>;
  using ParticleSpecies = onika::memory::CudaMMVector<ParticleSpecie>;

  struct MinimalParticleSpecies
  {
    double m_mass = 0.;
    double m_charge = 0.;
    unsigned int m_z = 0;
    MinimalParticleSpecies() = default;
    MinimalParticleSpecies(const MinimalParticleSpecies&) = default;
    MinimalParticleSpecies(MinimalParticleSpecies&&) = default;
    inline MinimalParticleSpecies(const ParticleSpecie& sp) : m_mass(sp.m_mass) , m_charge(sp.m_charge) , m_z(sp.m_z) {}
    MinimalParticleSpecies& operator = (const MinimalParticleSpecies& sp) = default;
    inline MinimalParticleSpecies& operator = (const ParticleSpecie& sp) { m_mass=sp.m_mass; m_charge=sp.m_charge; m_z=sp.m_z; return *this; }
  };


  template<class StreamT>
  inline StreamT& print_user_input(const ParticleSpecies& species, StreamT& out, int indent=0)
  {
    using onika::math::norm2;
    auto space = [indent](unsigned int n) -> std::string { return std::string((indent+n)*2,' '); } ;
    for(const auto & atom : species)
    {
      out << std::setprecision(15)
          << space(0) << "- "<<atom.name()<<":" << std::endl
          << space(2) << "mass: "<<atom.m_mass<<std::endl
          << space(2) << "z: "<<atom.m_z<<std::endl
          << space(2) << "charge: "<<atom.m_charge<<std::endl;
      if(!atom.molecule_name().empty()) { out<< space(2) << "molecule: "<<atom.molecule_name()<<std::endl; }
      if( atom.m_rigid_atom_count > 1 )
      {
        double max_radius = 0.0;
        for(unsigned int i=0;i<atom.m_rigid_atom_count;i++) max_radius = std::max( max_radius , norm2( atom.m_rigid_atoms[i].m_pos ) );
        out << space(2) << "radius: "<<sqrt(max_radius) << std::endl
            << space(2) << "rigid_molecule:" << std::endl;
        for(unsigned int i=0;i<atom.m_rigid_atom_count;i++)
        {
          out << space(3) << "- "<<atom.rigid_atom_name(i)<<": ["<<atom.m_rigid_atoms[i].m_pos<<"]"<<std::endl;
        }
      }
    }
    return out;
  }

}

