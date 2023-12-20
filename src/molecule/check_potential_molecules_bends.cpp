#include <yaml-cpp/yaml.h>
#include <memory>
#include <utility>// std::pair

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>

#include <exanb/core/basic_types_yaml.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/log.h>

#include <exaStamp/molecule/bends_potentials_parameters.h>

namespace exaStamp
{

  class CheckBendPotentials : public OperatorNode
  {    
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( ChemicalAngles           , chemical_angles       , INPUT, OPTIONAL );
    ADD_SLOT( BendsPotentialParameters , potentials_for_angles , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies          , species               , INPUT, REQUIRED );

  public:
    inline void execute ()  override final
    {
      // Initial checks
      if( ! chemical_angles.has_value() )
      {
        lerr << "chemical_angles input missing" << std::endl;
        std::abort();
      }
      
      // Re-build the map containing the bending potential parameters from potentials_for_angles->m_potentials
      //(read automatically by a YAML routine) in case the block ReadBendPotentials was not called before
      if( potentials_for_angles->m_type_to_potential.empty() )
      {
        std::unordered_map<std::string, unsigned int> map_species_name_id;
        for(size_t i=0;i<species->size(); ++i)
        {
          map_species_name_id[species->at(i).m_name] = i;
        }        
        for(const auto& b_type : potentials_for_angles->m_potentials)
        {          
          uint64_t type_a = map_species_name_id.at(b_type.species.at(0));
          uint64_t type_b = map_species_name_id.at(b_type.species.at(1));
          uint64_t type_c = map_species_name_id.at(b_type.species.at(2));
          if( type_a > type_c ) std::swap( type_a, type_c );
          assert( type_a < 65536 );
          assert( type_b < 65536 );
          assert( type_c < 65536 );
          uint64_t key = (type_a<<32) | (type_b<<16) | type_c ;
          assert( potentials_for_angles->m_type_to_potential.find(key) == potentials_for_angles->m_type_to_potential.end() );
          potentials_for_angles->m_type_to_potential[key] = b_type.m_potential_function;
        }
      }

      // Charge the angle potential parameters
      const auto& map_bends_potential = potentials_for_angles->m_type_to_potential;

      // Charge the bending angles
      ChemicalAngles&           bends_list  = *chemical_angles;
      const size_t n_bends = bends_list.size();

      // Loop to check whther all the angles detected in the system have their corresponding angle potential parameters
#     pragma omp parallel
      {
        size_t cell[3];
        size_t pos[3];
        unsigned int type[3];

#       pragma omp for schedule(dynamic)
        for(size_t i=0; i<n_bends; ++i)
        {
          //-----------------------------DECODE----------------------------------------------------
          uint64_t atom_to_decode_a =  bends_list[i][0]; // atom_from_idmap(bends_list[i][0],id_map,id_map_ghosts); 
          uint64_t atom_to_decode_b =  bends_list[i][1]; // atom_from_idmap(bends_list[i][1],id_map,id_map_ghosts); 
          uint64_t atom_to_decode_c =  bends_list[i][2]; // atom_from_idmap(bends_list[i][2],id_map,id_map_ghosts); 

          assert(atom_to_decode_a != std::numeric_limits<uint64_t>::max());
          assert(atom_to_decode_b != std::numeric_limits<uint64_t>::max());
          assert(atom_to_decode_c != std::numeric_limits<uint64_t>::max());

          decode_cell_particle(atom_to_decode_a, cell[0], pos[0], type[0]);
          decode_cell_particle(atom_to_decode_b, cell[1], pos[1], type[1]);
          decode_cell_particle(atom_to_decode_c, cell[2], pos[2], type[2]);

          // Get potential type from species of atoms triplet
          uint64_t type_a = type[0];
          uint64_t type_b = type[1];
          uint64_t type_c = type[2];
          if( type_a > type_c ) std::swap( type_a, type_c );
          assert( type_a < 65536 );
          assert( type_b < 65536 );
          assert( type_c < 65536 );
          uint64_t key = (type_a<<32) | (type_b<<16) | type_c ;
          auto it = map_bends_potential.find( key );
          if (it == map_bends_potential.end() ) {
            lerr << "ERROR: Parameters for BENDING ANGLE ["<< species->at(type_a).m_name << " ; " << species->at(type_b).m_name << " ; " << species->at(type_c).m_name << "] are missing!" << std::endl;
            std::abort();
          }

        }
      }

    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "check_bend_potentials", make_simple_operator< CheckBendPotentials > );
  }

}

