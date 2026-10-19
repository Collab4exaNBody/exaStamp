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

#include <string>
#include <cstdint>
#include <vector>

#include <exaStamp/particle_species/particle_specie.h>
#include <onika/log.h>
#include <onika/math/basic_types.h>

#include <onika/memory/allocator.h>
#include <onika/cuda/ro_shallow_copy.h>

#include <yaml-cpp/yaml.h>

#ifndef XSTAMP_MAX_MOLECULE_ATOMS
#define XSTAMP_MAX_MOLECULE_ATOMS 16
#endif

#ifndef XSTAMP_MAX_MOLECULE_BONDS
#define XSTAMP_MAX_MOLECULE_BONDS 16
#endif

#ifndef XSTAMP_MAX_MOLECULE_BENDS
#define XSTAMP_MAX_MOLECULE_BENDS 16
#endif

#ifndef XSTAMP_MAX_MOLECULE_TORSIONS
#define XSTAMP_MAX_MOLECULE_TORSIONS 16
#endif

#ifndef XSTAMP_MAX_MOLECULE_IMPROPERS
#define XSTAMP_MAX_MOLECULE_IMPROPERS 16
#endif

#ifndef XSTAMP_MAX_MOLECULE_PAIRS
#define XSTAMP_MAX_MOLECULE_PAIRS 32
#endif

namespace exaStamp
{
  using namespace exanb;

  static constexpr size_t MAX_MOLECULE_ATOMS = XSTAMP_MAX_MOLECULE_ATOMS;
  static constexpr size_t MAX_MOLECULE_BONDS = XSTAMP_MAX_MOLECULE_BONDS;
  static constexpr size_t MAX_MOLECULE_BENDS = XSTAMP_MAX_MOLECULE_BENDS;
  static constexpr size_t MAX_MOLECULE_TORSIONS = XSTAMP_MAX_MOLECULE_TORSIONS;
  static constexpr size_t MAX_MOLECULE_IMPROPERS = XSTAMP_MAX_MOLECULE_IMPROPERS;
  static constexpr size_t MAX_MOLECULE_PAIRS = 256; //XSTAMP_MAX_MOLECULE_PAIRS;

  struct alignas(8) MoleculeSpecies
  {
    unsigned int m_nb_atoms = 0;
    unsigned int m_nb_bonds = 0;
    unsigned int m_nb_bends = 0;
    unsigned int m_nb_torsions = 0;
    unsigned int m_nb_impropers = 0;

    int8_t m_atom_connectivity[MAX_MOLECULE_ATOMS][4]; // -1 for empty slots
    uint8_t m_torsions[MAX_MOLECULE_TORSIONS][4];
    uint8_t m_impropers[MAX_MOLECULE_IMPROPERS][4];
    uint8_t m_bonds[MAX_MOLECULE_BONDS][2];
    uint8_t m_bends[MAX_MOLECULE_BENDS][3];
    uint8_t m_atom_type[MAX_MOLECULE_ATOMS];
    
    char m_name[MAX_ATOM_NAME_LEN] = {'\0'};

    inline void set_name(const std::string& s)
    {
      if( s.length() >= MAX_ATOM_NAME_LEN ) { fatal_error()<<"Molecule name too long : length="<<s.length()<<", max="<<(MAX_ATOM_NAME_LEN-1)<<"\n"; std::abort(); }
      std::strncpy(m_name,s.c_str(),MAX_ATOM_NAME_LEN); m_name[MAX_ATOM_NAME_LEN-1]='\0';
    }
    inline std::string name() const { return m_name; }
    
    inline bool has_connectivity() const { return m_nb_atoms==0 || m_nb_bonds>0 || m_nb_bends>0 || m_nb_torsions>0 || m_nb_impropers>0; }
    void update_connectivity();
    
    template<class StreamT> inline void print(StreamT& out, const ParticleSpecies& species) const
    {
      out << "Molecule "<<name()<<" : "<<m_nb_atoms<<" atoms" <<std::endl;
      for(unsigned int i=0;i<m_nb_atoms;i++)
      {
        out << "  atom "<<i<<" ("<<species[m_atom_type[i]].name()<<") : {";
        for(int j=0;j<4;j++) if(m_atom_connectivity[i][j]!=-1) out << " " << int(m_atom_connectivity[i][j]);
        out << " }"<<std::endl;        
      }
      out << "  "<< m_nb_bonds<<" bonds :";
      for(unsigned int i=0;i<m_nb_bonds;i++)
      {
        out<<" {"<<int(m_bonds[i][0])<<","<<int(m_bonds[i][1])<<"}";
      }
      out << std::endl << "  "<<m_nb_bends <<" angles :";
      for(unsigned int i=0;i<m_nb_bends;i++)
      {
        out<<" {"<<int(m_bends[i][0])<<","<<int(m_bends[i][1])<<","<<int(m_bends[i][2])<<"}";
      }
      out << std::endl << "  "<< m_nb_torsions<<" torsions :";
      for(unsigned int i=0;i<m_nb_torsions;i++)
      {
        out<<" {"<<int(m_torsions[i][0])<<","<<int(m_torsions[i][1])<<","<<int(m_torsions[i][2])<<","<<int(m_torsions[i][3])<<"}";
      }
      out << std::endl << "  "<< m_nb_impropers<<" impropers :";
      for(unsigned int i=0;i<m_nb_impropers;i++)
      {
        out<<" {"<<int(m_impropers[i][0])<<","<<int(m_impropers[i][1])<<","<<int(m_impropers[i][2])<<","<<int(m_impropers[i][3])<<"}";
      }
      out<<std::endl;
    }
  };

  struct MoleculeSpeciesVector
  {
    onika::memory::CudaMMVector<MoleculeSpecies> m_molecules;
    onika::memory::CudaMMVector<MoleculeSpecies> m_bridge_molecules;  

    // molecule instances : each particle's idmol field is a dense molecule index, 0 <= idmol < number of molecules.
    // m_molecule_type[idmol] is the type of the molecule (index in m_molecules). An empty table means type 0 for all.
    std::vector<uint8_t> m_molecule_type;

    // continuous unfolding (molecule_unfold) : position of the anchor atom (smallest atom id) of each molecule at the
    // last unfolding, in grid space, indexed by idmol. m_unfold_anchor_id[idmol] is the anchor atom id + 1, 0 if unknown.
    std::vector<uint64_t> m_unfold_anchor_id;
    std::vector<onika::math::Vec3d> m_unfold_anchor_pos;

    inline size_t number_of_molecule_instances() const { return m_molecule_type.size(); }
    inline unsigned int molecule_type( uint64_t idmol ) const
    {
      if( m_molecule_type.empty() ) return 0;
      if( idmol >= m_molecule_type.size() )
      {
        fatal_error() << "molecule id "<<idmol<<" out of range (number of molecules = "<<m_molecule_type.size()<<")" << std::endl;
      }
      return m_molecule_type[idmol];
    }
    inline int molecule_type_index( const std::string& name ) const
    {
      for(size_t i=0;i<m_molecules.size();i++) if( m_molecules[i].name() == name ) return i;
      return -1;
    }
  };

  // folded positions saved by molecule_unfold, restored by molecule_restore_positions (owned particles, cell order)
  struct MoleculeFoldedPositions
  {
    std::vector<double> m_rx, m_ry, m_rz;
  };

  struct MoleculeSpeciesVectorRO
  {
    const MoleculeSpecies * __restrict__ m_molecules = nullptr;
    const MoleculeSpecies * __restrict__ m_bridge_molecules = nullptr;  

    MoleculeSpeciesVectorRO() = default;
    MoleculeSpeciesVectorRO(const MoleculeSpeciesVectorRO&) = default;
    MoleculeSpeciesVectorRO(MoleculeSpeciesVectorRO&&) = default;
    MoleculeSpeciesVectorRO& operator = (const MoleculeSpeciesVectorRO&) = default;
    MoleculeSpeciesVectorRO& operator = (MoleculeSpeciesVectorRO&&) = default;

    inline MoleculeSpeciesVectorRO( const MoleculeSpeciesVector& ml )
      : m_molecules( ml.m_molecules.data() )
      , m_bridge_molecules( ml.m_bridge_molecules.data() )
      {}
  };
  
  //using MoleculeSpeciesVector = std::vector<MoleculeSpecies>; // onika::memory::CudaMMVector<MoleculeSpecies>;
  

}

// specialize ReadOnlyShallowCopyType so MoleculeSpeciesVectorRO is the read only type for MoleculeSpeciesVector
namespace onika { namespace cuda { template<> struct ReadOnlyShallowCopyType< exaStamp::MoleculeSpeciesVector > { using type = exaStamp::MoleculeSpeciesVectorRO; }; } }


namespace YAML
{
  using exanb::fatal_error;
  using exaStamp::MoleculeSpecies;
  using exaStamp::MoleculeSpeciesVector;

  template<> struct convert<MoleculeSpecies>
  {
    static bool decode(const Node& node, MoleculeSpecies& p)
    {
      if( ! node.IsMap() ) return false;
      p = MoleculeSpecies{};
      if( node["name"] )
      {
        p.set_name( node["name"].as<std::string>() );
      }
      return true;
    }
  };

  template<> struct convert<MoleculeSpeciesVector>
  {
    static bool decode(const Node& node, MoleculeSpeciesVector& p)
    {
      p.m_molecules.clear();
      p.m_bridge_molecules.clear();
      if( ! node.IsSequence() ) return false;
      const auto mspvec = node.as< std::vector<MoleculeSpecies> >();
      p.m_molecules.reserve( mspvec.size() );
      for(const auto& msp : mspvec) p.m_molecules.push_back( msp );
      return true;
    }
  };
  
}
