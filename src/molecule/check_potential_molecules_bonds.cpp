#include <yaml-cpp/yaml.h>
#include <memory>
#include <utility>// std::pair
#include <iomanip>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>

#include <exanb/core/basic_types_yaml.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/log.h>

#include <exaStamp/molecule/bonds_potentials_parameters.h>

namespace exaStamp
{

  class CheckBondPotentials : public OperatorNode
  {

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( ChemicalBonds           , chemical_bonds        , INPUT, OPTIONAL );
    ADD_SLOT( BondsPotentialParameters, potentials_for_bonds  , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies         , species               , INPUT, REQUIRED );

  public:
    inline void execute() override final
    {

      // Initial checks
      if( ! chemical_bonds.has_value() )
      {
        lerr << "chemical_bonds input missing" << std::endl;
        std::abort();
      }

      // Re-build the map containing the bond potential parameters from potentials_for_bonds->m_bond_desc
      // (read automatically by a YAML routine) in case the block ReadBondPotentials was not called before
      if( potentials_for_bonds->m_type_pair_to_potential.empty() )
      {
        // Convert atom type name (string) to type (unsigned int)
        std::unordered_map<std::string, unsigned int> map_species_name_id;
        for(size_t i=0;i<species->size(); ++i)
        {
          map_species_name_id[species->at(i).m_name] = i;
        }

        // pre-build map of bond type pair to bond potential function
        for( const auto& b_type : potentials_for_bonds->m_bond_desc )
        {
          uint64_t a = map_species_name_id.at(b_type.species.at(0));
          uint64_t b = map_species_name_id.at(b_type.species.at(1));
          if( a > b ) std::swap( a , b );
          assert( a < 65536 );
          assert( b < 65536 );
          potentials_for_bonds->m_type_pair_to_potential[ (a<<16) | b ] = b_type.potential;
        }
      }

      // Charge the bond potential parameters
      const auto& map_bonds_potential = potentials_for_bonds->m_type_pair_to_potential;

      // Charge the bonds 
      ChemicalBonds&            bonds_list    = *chemical_bonds;
      size_t n_bonds = bonds_list.size();
      
      // Loop to check whether all bonds detected in the system have their bond potential parameters
#     pragma omp parallel //for shared(map_bonds_potential) private(cell,pos,type)
      {
        size_t       cell[2];
        size_t       pos[2];
        unsigned int type[2];

#       pragma omp for
        for(size_t i=0; i<n_bonds; ++i)
        {
          // ---- bonds_list contains the local ids of the bonded atoms -----
          //lerr << "In src/molecule/check_potential_molecules_bonds.cpp : bonds_list[i][0] =  " << (uint64_t) bonds_list[i].at(0) << " bonds_list[i][1] = " << (uint64_t) bonds_list[i].at(1) << std::endl << std::flush ;
          const uint64_t atom_to_decode_a = bonds_list[i][0]; //atom_from_idmap( bonds_list[i][0] ,id_map, id_map_ghosts ); 
          const uint64_t atom_to_decode_b = bonds_list[i][1]; //atom_from_idmap( bonds_list[i][1] ,id_map, id_map_ghosts ); 
          

          decode_cell_particle(atom_to_decode_a, cell[0], pos[0], type[0]);
          decode_cell_particle(atom_to_decode_b, cell[1], pos[1], type[1]);


          // Get potential type from species of atoms pair
          uint64_t type_a = type[0];
          uint64_t type_b = type[1];
          if( type_a > type_b ) std::swap( type_a , type_b );          
          auto it = map_bonds_potential.find( ( type_a << 16 ) | type_b );
          if (it == map_bonds_potential.end())
          {
            lerr << "ERROR: Parameters for BOND ["<< species->at(type_a).m_name << " ; " << species->at(type_b).m_name << "] are missing!" << std::endl;
            std::abort();
          }

        }
      }

    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "check_bond_potentials", make_simple_operator< CheckBondPotentials > );
  }

}

