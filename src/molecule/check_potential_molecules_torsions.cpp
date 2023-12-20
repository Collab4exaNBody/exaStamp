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

#include <exaStamp/molecule/torsions_potentials_parameters.h>

namespace exaStamp
{

  class CheckTorsionPotentials : public OperatorNode
  {
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( ChemicalTorsions            , chemical_torsions       , INPUT, OPTIONAL );
    ADD_SLOT( TorsionsPotentialParameters , potentials_for_torsions , INPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies             , species                 , INPUT, REQUIRED );
    
  public:
    inline void execute ()  override final
    {


      // Initial checks
      if( ! chemical_torsions.has_value() )
      {
        lerr << "chemical_torsions input missing" << std::endl;
        std::abort();
      }

      ParticleSpecies& species = *(this->species);

      // Re-build the map containing the torsion potential parameters from potentials_for_torsions->m_potentials
      // (read automatically by a YAML routine) in case the block ReadTorsionPotentials was not called before
      if( potentials_for_torsions->m_type_to_potential.empty() )
      {
        // Convert name (string) to type (uint_8)
        std::map<std::string, uint8_t> map_species_name_id;
        for(size_t i=0;i<species.size(); ++i)
        {
          map_species_name_id[species.at(i).m_name] = i;
        }

        // Creation map torsions potentials
        //std::map<std::array<unsigned int, 4>, std::shared_ptr<IntraMolecularPotentialFunctional> > map_torsions_potential;
        for(const auto& b_type : potentials_for_torsions->m_potentials )
        {
          uint64_t t0 = map_species_name_id.at(b_type.species.at(0));
          uint64_t t1 = map_species_name_id.at(b_type.species.at(1));
          uint64_t t2 = map_species_name_id.at(b_type.species.at(2));
          uint64_t t3 = map_species_name_id.at(b_type.species.at(3));
          uint64_t key = (t0<<48) | (t1<<32) | (t2<<16) | t3;
          uint64_t alt_key = (t3<<48) | (t2<<32) | (t1<<16) | t0;
          if( alt_key < key ) key = alt_key;
          potentials_for_torsions->m_type_to_potential[key] = b_type.m_potential_function;
        }
      }

      // XXX
      // Charge the torsion angles
      ChemicalTorsions&     torsions_list         = *chemical_torsions;
      size_t n_torsion = torsions_list.size();

#     pragma omp parallel 
      {
        size_t cell[4];
        size_t pos[4];
        unsigned int type[4];
        
#       pragma omp for schedule(static)
        for(size_t i=0; i<n_torsion; ++i)
        {
          //-----------------------------DECODE----------------------------------------------------
          uint64_t atom_to_decode_a = torsions_list[i][0]; // atom_from_idmap(id1,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_b = torsions_list[i][1]; // atom_from_idmap(id2,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_c = torsions_list[i][2]; // atom_from_idmap(id3,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_d = torsions_list[i][3]; // atom_from_idmap(id4,id_map,id_map_ghosts); 

          decode_cell_particle(atom_to_decode_a, cell[0], pos[0], type[0]);
          decode_cell_particle(atom_to_decode_b, cell[1], pos[1], type[1]);
          decode_cell_particle(atom_to_decode_c, cell[2], pos[2], type[2]);
          decode_cell_particle(atom_to_decode_d, cell[3], pos[3], type[3]);
          
          uint64_t t0 = type[0];
          uint64_t t1 = type[1];
          uint64_t t2 = type[2];
          uint64_t t3 = type[3];
          uint64_t key = (t0<<48) | (t1<<32) | (t2<<16) | t3;
          uint64_t alt_key = (t3<<48) | (t2<<32) | (t1<<16) | t0;
          if( alt_key < key ) key = alt_key;
          auto it = potentials_for_torsions->m_type_to_potential.find(key);
          if (it == potentials_for_torsions->m_type_to_potential.end()) {
            lerr << "ERROR: Parameters for TORSION ANGLE ["<< species.at(t0).m_name << " ; " << species.at(t1).m_name << " ; " << species.at(t2).m_name <<  " ; " << species.at(t3).m_name << "] are missing!" << std::endl;
            std::abort();
          }

        }
      }
    }

  };



  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "check_torsion_potentials", make_simple_operator< CheckTorsionPotentials> );
  }

}


