#pragma once

#include <string>
#include <cstdint>

#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/log.h>

#include <onika/memory/allocator.h>
#include <onika/cuda/ro_shallow_copy.h>

#include <yaml-cpp/yaml.h>

namespace exaStamp
{
  using namespace exanb;

  static constexpr size_t MAX_MOLECULE_ATOMS = 16;
  static constexpr size_t MAX_MOLECULE_BONDS = 16;
  static constexpr size_t MAX_MOLECULE_BENDS = 16;
  static constexpr size_t MAX_MOLECULE_TORSIONS = 16;
  static constexpr size_t MAX_MOLECULE_IMPROPERS = 16;
  static constexpr size_t MAX_MOLECULE_PAIRS = 32;

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
    
    template<class StreamT> inline void print(StreamT& out, const ParticleSpecies& species)
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
  
  static inline uint8_t molecule_type_from_id( uint64_t idmol )
  {
    return idmol & 0xFF ;
  }

  static inline uint8_t molecule_place_from_id( uint64_t idmol )
  {
    return (idmol>>8) & 0xFF ;
  }

  static inline uint64_t molecule_instance_from_id( uint64_t idmol )
  {
    return idmol >> 16 ;
  }

  static inline uint64_t make_molecule_id( uint64_t mol_instance, unsigned int mol_place, unsigned int mol_type )
  {
    return ( mol_instance << 16 ) | ( (mol_place&0xFF) << 8 ) | ( mol_type&0xFF ) ;
  }

  // helper class to read and write molecule species to dump files
  struct MoleculeOptionalHeaderIO
  {
    double& m_bond_max_dist;
    double& m_bond_max_stretch;
    MoleculeSpeciesVector* m_molecules = nullptr;
    template<class WriteFuncT> inline size_t write_optional_header( WriteFuncT write_func )
    {
      static const uint64_t header_version = 100;
      const int vermaj = header_version/100;
      const int vermin = header_version%100;
      size_t n_bytes = 0;
      size_t n_molecules = 0;
      if( m_molecules != nullptr ) n_molecules = m_molecules->m_molecules.size();
      lout << "molecule ext. v"<<vermaj<<'.'<<vermin <<std::endl;
      lout << "mol. max dist = "<<m_bond_max_dist<<std::endl;
      lout << "mol. stretch  = "<<m_bond_max_stretch<<std::endl;
      lout << "mol. species  = "<<n_molecules<<std::endl;
      n_bytes += write_func( header_version );
      n_bytes += write_func( m_bond_max_dist );
      n_bytes += write_func( m_bond_max_stretch );
      n_bytes += write_func( n_molecules );
      for(size_t i=0;i<n_molecules;i++)
      {
        n_bytes += write_func( m_molecules->m_molecules.at(i) );
      }
      return n_bytes;
    }
    template<class ReadFuncT> inline size_t read_optional_header( ReadFuncT read_func )
    {
      uint64_t header_version = 100;
      size_t n_bytes = 0;
      size_t n_molecules = 0;
      m_bond_max_dist = m_bond_max_stretch = 0.0;
      n_bytes += read_func( header_version );
      const int vermaj = header_version/100;
      const int vermin = header_version%100;
      if( header_version > 100 ) { fatal_error() << "Molecule header version "<<vermaj<<'.'<<vermin<<" not supported by this software version"<<std::endl; }
      n_bytes += read_func( m_bond_max_dist );
      n_bytes += read_func( m_bond_max_stretch );
      n_bytes += read_func( n_molecules );
      lout << "molecule ext. v"<<vermaj<<'.'<<vermin <<std::endl;
      lout << "mol. max dist = "<<m_bond_max_dist<<std::endl;
      lout << "mol. stretch  = "<<m_bond_max_stretch<<std::endl;
      lout << "mol. species  = "<<n_molecules<<std::endl;
      if( n_molecules>0 && m_molecules==nullptr )
      {
        fatal_error() << "Missing molecules container to read molecule species" << std::endl;
      }
      if( m_molecules!=nullptr ) m_molecules->m_molecules.resize( n_molecules );
      for(size_t i=0;i<n_molecules;i++)
      {
        n_bytes += read_func( m_molecules->m_molecules.at(i) );
      }
      ldbg<<"Molecule dump header : bytes read ="<<n_bytes<<std::endl;
      return n_bytes;
    }
  };

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
