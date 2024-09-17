#pragma once

#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/impropers_potentials_parameters.h>
#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>
#include <exaStamp/molecule/intramolecular_pair_weight.h>

namespace exaStamp
{
  using namespace exanb;

  // helper class to read and write molecule species to dump files
  template<class LDBGT>
  struct MoleculeOptionalHeaderIO
  {
    static inline constexpr uint64_t header_version = 101; // Version 1.1 adds potential and weighting parameters in optional header
  
    // since version 1.0
    double& m_bond_max_dist;
    double& m_bond_max_stretch;
    MoleculeSpeciesVector * m_molecules = nullptr;
    
    LDBGT& ldbg;
    
    // since version 1.1
    BondsPotentialParameters * m_potentials_for_bonds = nullptr;
    BendsPotentialParameters * m_potentials_for_angles = nullptr;
    TorsionsPotentialParameters * m_potentials_for_torsions = nullptr;
    ImpropersPotentialParameters * m_potentials_for_impropers = nullptr;
    LJExp6RFMultiParms * m_potentials_for_pairs = nullptr;
    IntramolecularPairWeighting * m_mol_pair_weights = nullptr;
    
    struct IntramolecularFunctionalIO
    {
      char m_type[4][32];
      MoleculeGenericFuncParam m_params = {};
    };

    struct IntramolecularPairIO
    {
      char m_type_a[32];
      char m_type_b[32];
      LJExp6RFParms m_params;
    };
    
    struct IntramolecularPairWeightIO
    {
      char m_name[32];
      MolecularPairWeight m_weights;
    };

    template<class StreamT>
    inline void print(StreamT& out, const ParticleSpecies& species) const
    {
      const int vermaj = header_version/100;
      const int vermin = header_version%100;
      out << std::defaultfloat << std::endl;
      out << "molecule ext. v"<<vermaj<<'.'<<vermin <<std::endl;
      out << "mol. max dist = "<<m_bond_max_dist<<std::endl;
      out << "mol. stretch  = "<<m_bond_max_stretch<<std::endl;
      out << std::setprecision(6) << std::endl;
      if( m_molecules != nullptr )
      {
        out << "mol. species  = "<<m_molecules->m_molecules.size()<<std::endl;
        for(const auto& mol : m_molecules->m_molecules)
        {
          mol.print( out , species );
        }
      }
      if( header_version >= 101 )
      {
        if( m_potentials_for_bonds != nullptr )
        {
          out << "bond params = "<<m_potentials_for_bonds->m_bond_desc.size()<<std::endl;
          for(const auto& bond : m_potentials_for_bonds->m_bond_desc)
          {
            out<<"\t"<<bond.species[0]<<"/"<<bond.species[1]<<" : "<<bond.type<<" = "<< bond.potential->generic_parameters()<<std::endl;
          }
        }
        if( m_potentials_for_angles != nullptr )
        {
          out << "angle params  = "<<m_potentials_for_angles->m_potentials.size()<<std::endl;
          for(const auto& angle : m_potentials_for_angles->m_potentials)
          {
            out<<"\t"<<angle.species[0]<<"/"<<angle.species[1]<<"/"<<angle.species[2]<<" : "<<angle.type<<" = "<< angle.m_potential_function->generic_parameters() <<std::endl;
          }
        }
        if( m_potentials_for_torsions != nullptr )
        {
          out << "torsion params = "<<m_potentials_for_torsions->m_potentials.size()<<std::endl;
          for(const auto& torsion : m_potentials_for_torsions->m_potentials)
          {
            out<<"\t"<<torsion.species[0]<<"/"<<torsion.species[1]<<"/"<<torsion.species[2]<<"/"<<torsion.species[3]<<" : "<<torsion.type<<" = "<< torsion.m_potential_function->generic_parameters() <<std::endl;
          }
        }
        if( m_potentials_for_impropers != nullptr )
        {
          out << "improper params = "<<m_potentials_for_impropers->m_potentials.size()<<std::endl;
          for(const auto& improper : m_potentials_for_impropers->m_potentials)
          {
            out<<"\t"<<improper.species[0]<<"/"<<improper.species[1]<<"/"<<improper.species[2]<<"/"<<improper.species[3]<<" : "<<improper.type<<" = "<< improper.m_potential_function->generic_parameters() <<std::endl;
          }
        }
        if( m_potentials_for_pairs != nullptr )
        {
          out << "pair params = "<<m_potentials_for_pairs->m_potentials.size()<<std::endl;
          for(const auto& pairpot : m_potentials_for_pairs->m_potentials)
          {
            out<<"\t"<<pairpot.m_type_a[0]<<"/"<<pairpot.m_type_b[1]<<" = "; pairpot.m_params.to_stream(out); out<<std::endl;
          }
        }
        if( m_mol_pair_weights != nullptr )
        {
          out << "pair weights  = "<<m_mol_pair_weights->m_molecule_weight.size()<<std::endl;
          for(const auto& elem : m_mol_pair_weights->m_molecule_weight)
          {
            out <<"\t"<<elem.first<<" : PAIR : bond="<<elem.second.m_bond_weight<<", angle="<<elem.second.m_bend_weight<<", torsion="<<elem.second.m_torsion_weight
                                <<" : RF : bond="<<elem.second.m_rf_bond_weight<<", rf angle="<<elem.second.m_rf_bend_weight<<", rf torsion="<<elem.second.m_rf_torsion_weight<<std::endl;
          }
        }
      }    
    }

    template<class WriteFuncT> inline size_t write_optional_header( WriteFuncT write_func , bool verbose = true )
    {
      const int vermaj = header_version/100;
      const int vermin = header_version%100;
      size_t n_bytes = 0;
      size_t n_molecules = 0;
      if( m_molecules != nullptr ) n_molecules = m_molecules->m_molecules.size();
      if(verbose) lout << "molecule ext. v"<<vermaj<<'.'<<vermin <<std::endl;
      if(verbose) lout << "mol. max dist = "<<m_bond_max_dist<<std::endl;
      if(verbose) lout << "mol. stretch  = "<<m_bond_max_stretch<<std::endl;
      if(verbose) lout << "mol. species  = "<<n_molecules<<std::endl;
      n_bytes += write_func( header_version );
      n_bytes += write_func( m_bond_max_dist );
      n_bytes += write_func( m_bond_max_stretch );
      n_bytes += write_func( n_molecules );
      for(size_t i=0;i<n_molecules;i++)
      {
        n_bytes += write_func( m_molecules->m_molecules.at(i) );
      }

      if( header_version >= 101 )
      {
        const size_t n_bond_potentials = ( m_potentials_for_bonds != nullptr ) ? m_potentials_for_bonds->m_bond_desc.size() : 0 ;
        n_bytes += write_func( n_bond_potentials );
        if( m_potentials_for_bonds != nullptr )
        {
          ldbg << "bond params   = "<<n_bond_potentials<<std::endl;
          for(const auto& bond : m_potentials_for_bonds->m_bond_desc)
          {
            IntramolecularFunctionalIO descio;
            for(int i=0;i<2;i++) std::strncpy( descio.m_type[i], bond.species[i].c_str(), 32 );
            descio.m_params = bond.potential->generic_parameters();
            n_bytes += write_func( descio );
          }
        }

        const size_t n_angle_potentials = ( m_potentials_for_angles != nullptr ) ? m_potentials_for_angles->m_potentials.size() : 0 ;
        n_bytes += write_func( n_angle_potentials );
        if( m_potentials_for_angles != nullptr )
        {
          ldbg << "angle params  = "<<n_angle_potentials<<std::endl;
          for(const auto& angle : m_potentials_for_angles->m_potentials)
          {
            IntramolecularFunctionalIO descio;
            for(int i=0;i<3;i++) std::strncpy( descio.m_type[i], angle.species[i].c_str(), 32 );
            descio.m_params = angle.m_potential_function->generic_parameters();
            n_bytes += write_func( descio );
          }
        }

        const size_t n_torsion_potentials = ( m_potentials_for_torsions != nullptr ) ? m_potentials_for_torsions->m_potentials.size() : 0 ;
        n_bytes += write_func( n_torsion_potentials );
        if( m_potentials_for_torsions != nullptr )
        {
          ldbg << "torsion param = "<<n_torsion_potentials<<std::endl;
          for(const auto& torsion : m_potentials_for_torsions->m_potentials)
          {
            IntramolecularFunctionalIO descio;
            for(int i=0;i<4;i++) std::strncpy( descio.m_type[i], torsion.species[i].c_str(), 32 );
            descio.m_params = torsion.m_potential_function->generic_parameters();
            n_bytes += write_func( descio );
          }
        }

        const size_t n_improper_potentials = ( m_potentials_for_impropers != nullptr ) ? m_potentials_for_impropers->m_potentials.size() : 0 ;
        n_bytes += write_func( n_improper_potentials );
        if( m_potentials_for_impropers != nullptr )
        {
          ldbg << "improper parm = "<<n_improper_potentials<<std::endl;
          for(const auto& improper : m_potentials_for_impropers->m_potentials)
          {
            IntramolecularFunctionalIO descio;
            for(int i=0;i<4;i++) std::strncpy( descio.m_type[i], improper.species[i].c_str(), 32 );
            descio.m_params = improper.m_potential_function->generic_parameters();
            n_bytes += write_func( descio );
          }
        }
        
        const size_t n_pair_potentials = ( m_potentials_for_pairs != nullptr ) ? m_potentials_for_pairs->m_potentials.size() : 0 ;
        n_bytes += write_func( n_pair_potentials );
        if( m_potentials_for_pairs != nullptr )
        {
          ldbg << "pair params   = "<<n_pair_potentials<<std::endl;
          for(const auto& pairpot : m_potentials_for_pairs->m_potentials)
          {
            IntramolecularPairIO descio;
            std::strncpy( descio.m_type_a, pairpot.m_type_a.c_str(), 32 );
            std::strncpy( descio.m_type_b, pairpot.m_type_b.c_str(), 32 );
            descio.m_params = pairpot.m_params;
            n_bytes += write_func( descio );
          }
        }

        const size_t n_pair_weights = ( m_mol_pair_weights != nullptr ) ? m_mol_pair_weights->m_molecule_weight.size() : 0;
        n_bytes += write_func( n_pair_weights );
        if( m_mol_pair_weights != nullptr )
        {
          ldbg << "pair weights  = "<<n_pair_weights<<std::endl;
          for(const auto& elem : m_mol_pair_weights->m_molecule_weight)
          {
            IntramolecularPairWeightIO descio;
            std::strncpy( descio.m_name, elem.first.c_str(), 32 );
            descio.m_weights = elem.second;
            n_bytes += write_func( descio );
          }
        }
      }

      return n_bytes;
    }

    inline size_t serialize_molecule_header( std::vector<uint8_t>& buffer )
    {
      buffer.clear();
      auto bufwrite = [&buffer]( auto x ) -> size_t
      {
        const uint8_t * obj_ptr = reinterpret_cast<const uint8_t *>( & x );
        buffer.insert( buffer.end() , obj_ptr , obj_ptr + sizeof(x) );
        return sizeof(x);
      };
      this->write_optional_header( bufwrite , false );
      return buffer.size();
    }

    template<class ReadFuncT> inline size_t read_optional_header( ReadFuncT read_func , bool verbose = true )
    {
      uint64_t header_version = 100;
      size_t n_bytes = 0;
      size_t n_molecules = 0;
      m_bond_max_dist = m_bond_max_stretch = 0.0;
      n_bytes += read_func( header_version );
      const int vermaj = header_version/100;
      const int vermin = header_version%100;
      n_bytes += read_func( m_bond_max_dist );
      n_bytes += read_func( m_bond_max_stretch );
      n_bytes += read_func( n_molecules );
      if(verbose) lout << "mol. header   = v"<<vermaj<<'.'<<vermin <<std::endl;
      if(verbose) lout << "mol. max dist = "<<m_bond_max_dist<<std::endl;
      if(verbose) lout << "mol. stretch  = "<<m_bond_max_stretch<<std::endl;
      if(verbose) lout << "mol. species  = "<<n_molecules<<std::endl;
      if( n_molecules>0 && m_molecules==nullptr )
      {
        fatal_error() << "Missing molecules container to read molecule species" << std::endl;
      }
      if( m_molecules!=nullptr ) m_molecules->m_molecules.resize( n_molecules );
      for(size_t i=0;i<n_molecules;i++)
      {
        n_bytes += read_func( m_molecules->m_molecules.at(i) );
      }

      if( header_version >= 101 )
      {
        size_t n_bond_potentials = 0;
        n_bytes += read_func( n_bond_potentials );
        if( n_bond_potentials > 0 )
        {
          ldbg << "bond params   = "<<n_bond_potentials<<std::endl;
          assert( m_potentials_for_bonds != nullptr );
          m_potentials_for_bonds->m_bond_desc.clear();
          for(size_t i=0;i<n_bond_potentials;i++)
          {
            IntramolecularFunctionalIO descio;
            n_bytes += read_func( descio );
            BondPotential pot = { "generic_bond" , {descio.m_type[0],descio.m_type[1]} , std::make_shared<IntraMolecularDumpFunctional>(descio.m_params) };
            m_potentials_for_bonds->m_bond_desc.push_back( pot );
          }
        }

        size_t n_angle_potentials = 0;
        n_bytes += read_func( n_angle_potentials );
        if( n_angle_potentials > 0 )
        {
          ldbg << "angle params  = "<<n_angle_potentials<<std::endl;
          assert( m_potentials_for_angles != nullptr );
          m_potentials_for_angles->m_potentials.clear();
          for(size_t i=0;i<n_angle_potentials;i++)
          {
            IntramolecularFunctionalIO descio;
            n_bytes += read_func( descio );
            BendPotential pot = { "generic_angle" , {descio.m_type[0],descio.m_type[1],descio.m_type[2]} , std::make_shared<IntraMolecularDumpFunctional>(descio.m_params) };
            m_potentials_for_angles->m_potentials.push_back( pot );
          }
        }

        size_t n_torsion_potentials = 0;
        n_bytes += read_func( n_torsion_potentials );
        if( n_torsion_potentials > 0 )
        {
          ldbg << "torsion param = "<<n_torsion_potentials<<std::endl;
          assert( m_potentials_for_torsions != nullptr );
          m_potentials_for_torsions->m_potentials.clear();
          for(size_t i=0;i<n_torsion_potentials;i++)
          {
            IntramolecularFunctionalIO descio;
            n_bytes += read_func( descio );
            TorsionPotential pot = { "generic_torsion" , {descio.m_type[0],descio.m_type[1],descio.m_type[2],descio.m_type[3]} , std::make_shared<IntraMolecularDumpFunctional>(descio.m_params) };
            m_potentials_for_torsions->m_potentials.push_back( pot );
          }
        }

        size_t n_improper_potentials = 0;
        n_bytes += read_func( n_improper_potentials );
        if( n_improper_potentials > 0 )
        {
          ldbg << "improper parm = "<<n_improper_potentials<<std::endl;
          assert( m_potentials_for_impropers != nullptr );
          m_potentials_for_impropers->m_potentials.clear();
          for(size_t i=0;i<n_improper_potentials;i++)
          {
            IntramolecularFunctionalIO descio;
            n_bytes += read_func( descio );
            ImproperPotential pot = { "generic_improper" , {descio.m_type[0],descio.m_type[1],descio.m_type[2],descio.m_type[3]} , std::make_shared<IntraMolecularDumpFunctional>(descio.m_params) };
            m_potentials_for_impropers->m_potentials.push_back( pot );
          }
        }

        size_t n_pair_potentials = 0;
        n_bytes += read_func( n_pair_potentials );
        if( n_pair_potentials > 0 )
        {
          ldbg << "pair params   = "<<n_pair_potentials<<std::endl;
          assert( m_potentials_for_pairs != nullptr );
          m_potentials_for_pairs->m_potentials.clear();
          for(size_t i=0;i<n_pair_potentials;i++)
          {
            IntramolecularPairIO descio;
            n_bytes += read_func( descio );
            LJExp6RFMultiParmsPair pot = { descio.m_type_a , descio.m_type_b , descio.m_params };
            m_potentials_for_pairs->m_potentials.push_back( pot );
          }
        }

        size_t n_pair_weights = 0;
        n_bytes += read_func( n_pair_weights );
        if( n_pair_weights > 0 )
        {
          ldbg << "pair weights  = "<<n_pair_weights<<std::endl;
          assert( m_mol_pair_weights != nullptr );
          m_mol_pair_weights->m_molecule_weight.clear();
          for(size_t i=0;i<n_pair_weights;i++)
          {
            IntramolecularPairWeightIO descio;
            n_bytes += read_func( descio );
            m_mol_pair_weights->m_molecule_weight.insert( { descio.m_name , descio.m_weights } );
          }
        }
      }

      ldbg<<"Molecule dump header : bytes read ="<<n_bytes<<std::endl;
      return n_bytes;
    }
    
    inline size_t deserialize_molecule_header( const std::vector<uint8_t>& buffer )
    {
      uint8_t * bufptr = buffer.data();
      size_t bufsz = buffer.size();
      auto bufread = [&bufptr,&bufsz]( auto& x ) -> size_t
      {
        assert( bufsz >= sizeof(x) );
        std::memcpy( reinterpret_cast<uint8_t*>(&x) , bufptr , sizeof(x) );
        bufsz -= sizeof(x);
        bufptr += sizeof(x);
        return sizeof(x);
      };
      size_t n = this->read_optional_header( bufread , false );
      assert( n == buffer.size() );
      return n;
    }

    
  };

}

