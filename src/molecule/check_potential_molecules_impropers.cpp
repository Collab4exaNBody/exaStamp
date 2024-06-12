#include <yaml-cpp/yaml.h>
#include <memory>
#include <utility>// std::pair
#include <omp.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>

#include <exanb/core/basic_types_yaml.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/log.h>

#include <exaStamp/molecule/impropers_potentials_parameters.h>

namespace exaStamp
{

  class CheckImproperPotentials : public OperatorNode
  {    
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( ChemicalImpropers            , chemical_impropers       , INPUT, OPTIONAL );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies              , species                  , INPUT, REQUIRED );

  public:
    inline void execute ()  override final
    {


      // Initial checks
      if( ! chemical_impropers.has_value() )
      {
        lerr << "chemical_impropers input missing" << std::endl;
        std::abort();
      }

      ParticleSpecies& species = *(this->species);

      // Re-build the map containing the improper potential parameters from potentials_for_impropers->m_potentials
      // (read automatically by a YAML routine) in case the block ReadImproperPotentials was not called before
      if( potentials_for_impropers->m_type_to_potential.empty() )
      {
        // Convert name (string) to type (uint_8)
        std::map<std::string, uint8_t> map_species_name_id;
        for(size_t i=0;i<species.size(); ++i)
        {
          map_species_name_id[species.at(i).m_name] = i;
        }

        // Creation map impropers potentials

        for(const auto& b_type : potentials_for_impropers->m_potentials)
        {
          uint64_t t0 = map_species_name_id.at(b_type.species.at(0));
          uint64_t t1 = map_species_name_id.at(b_type.species.at(1));
          uint64_t t2 = map_species_name_id.at(b_type.species.at(2));
          uint64_t t3 = map_species_name_id.at(b_type.species.at(3));
          // Atom 0 is the central atom
          // We consider that the 6 permutations 0 1 2 3
          //                                     0 2 3 1
          //                                     0 3 1 2
          //                                and  0 2 1 3
          //                                     0 1 3 2
          //                                     0 3 2 1
          // are equivalent
          uint64_t key     = (t0<<48) | (t1<<32) | (t2<<16) | t3;
          uint64_t alt_key = (t0<<48) | (t2<<32) | (t3<<16) | t1;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t3<<32) | (t1<<16) | t2;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t2<<32) | (t1<<16) | t3;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t1<<32) | (t3<<16) | t2;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t3<<32) | (t2<<16) | t1;
          if( alt_key < key ) key = alt_key;
          potentials_for_impropers->m_type_to_potential[key] = b_type.m_potential_function;

          const auto gp = b_type.m_potential_function->generic_parameters();
          ldbg << "Improper potential for "<<b_type.species.at(0)<<"/"<<b_type.species.at(1)<<"/"<<b_type.species.at(2)<<"/"<<b_type.species.at(3)<< " : p0="<<gp.p0<<", p1="<<gp.p1<<", p2="<<gp.p2<<", x0="<<gp.x0<<", coeff="<<gp.coeff<<std::endl;
        }
      }

      // XXX
      // Charge the improper torsion angles
      ChemicalImpropers&     impropers_list         = *chemical_impropers;
      size_t n_improper = impropers_list.size();

      unsigned int max_nt = omp_get_max_threads();
      std::vector<std::vector<uint64_t> > alreadyFound;
      alreadyFound.resize(max_nt);

#     pragma omp parallel
      {
        size_t cell[4];
        size_t pos[4];
        unsigned int type[4];
        int thread_id = omp_get_thread_num();

#       pragma omp for schedule(static)
        for(size_t i=0; i<n_improper; ++i)
        {

          //-----------------------------DECODE----------------------------------------------------
          uint64_t atom_to_decode_a = impropers_list[i][0]; // atom_from_idmap(id1,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_b = impropers_list[i][1]; // atom_from_idmap(id2,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_c = impropers_list[i][2]; // atom_from_idmap(id3,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_d = impropers_list[i][3]; // atom_from_idmap(id4,id_map,id_map_ghosts); 

          decode_cell_particle(atom_to_decode_a, cell[0], pos[0], type[0]);
          decode_cell_particle(atom_to_decode_b, cell[1], pos[1], type[1]);
          decode_cell_particle(atom_to_decode_c, cell[2], pos[2], type[2]);
          decode_cell_particle(atom_to_decode_d, cell[3], pos[3], type[3]);

          uint64_t t0 = type[0];
          uint64_t t1 = type[1];
          uint64_t t2 = type[2];
          uint64_t t3 = type[3];
          uint64_t key     = (t0<<48) | (t1<<32) | (t2<<16) | t3;
          uint64_t alt_key = (t0<<48) | (t2<<32) | (t3<<16) | t1;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t3<<32) | (t1<<16) | t2;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t2<<32) | (t1<<16) | t3;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t1<<32) | (t3<<16) | t2;
          if( alt_key < key ) key = alt_key;
          alt_key          = (t0<<48) | (t3<<32) | (t2<<16) | t1;
          if( alt_key < key ) key = alt_key;
          auto it = potentials_for_impropers->m_type_to_potential.find(key);

          // if the improper type is not found in the potential
          if (it == potentials_for_impropers->m_type_to_potential.end()) {
            // If the improper type is not either in the list "alreadyFound"
            if (std::find(alreadyFound.at(thread_id).begin(), alreadyFound.at(thread_id).end(), key) == alreadyFound.at(thread_id).end()) 
              //then the improper type is added to the list
              alreadyFound.at(thread_id).push_back(key);
            //lerr << "NOTE: Parameters for IMPROPER ANGLE ["<< species.at(t0).m_name << " ; " << species.at(t1).m_name <<  " ; " << species.at(t2).m_name << " ; " << species.at(t3).m_name << " ] were not found!" << std::endl;
            //lerr << "The corresponding energy and force are set to zero." << std::endl;
          }

        }
      }

      //Reduction of the improper angles not found in the potential for each thread
      for (size_t i=1; i < max_nt; i++) {
        for (size_t j=0; j < alreadyFound.at(i).size(); j++) {
          if (std::find(alreadyFound.at(0).begin(), alreadyFound.at(0).end(), alreadyFound.at(i).at(j))  == alreadyFound.at(0).end() )
            alreadyFound.at(0).push_back(alreadyFound.at(i).at(j));
        }

      }

      for (size_t j=0; j < alreadyFound.at(0).size(); j++) {
        //Extract the particle type index (t0, t1, t2, and t3) from  alreadyFound.at(0).at(j)
        size_t MAX_PARTICLE_TYPE = (1ull<<16) - 1;
        uint64_t t0, t1, t2, t3;
        uint64_t key = alreadyFound.at(0).at(j);

        t3 = key & MAX_PARTICLE_TYPE;
        key = key >> 16;
        t2 = key & MAX_PARTICLE_TYPE;
        key = key >> 16;
        t1 = key & MAX_PARTICLE_TYPE;
        key = key >> 16;
        t0 = key & MAX_PARTICLE_TYPE;

        //write warning message (not an error)
        lerr << "NOTE: Parameters for IMPROPER ANGLE ["<< species.at(t0).m_name << " ; " << species.at(t1).m_name <<  " ; " << species.at(t2).m_name << " ; " << species.at(t3).m_name << " ] were not found!" << std::endl;
        lerr << "The corresponding energy and force are set to zero." << std::endl;
        
      }

    }

  };



  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "check_improper_potentials", make_simple_operator< CheckImproperPotentials > );
  }

}


