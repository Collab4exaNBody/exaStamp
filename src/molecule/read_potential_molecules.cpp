#include <yaml-cpp/yaml.h>
#include <memory>
#include <utility>// std::pair
#include <iomanip>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>

#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/log.h>

#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/impropers_potentials_parameters.h>

namespace exaStamp
{
  using namespace exanb;

  class ReadBondPotentials : public OperatorNode
  {

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( BondsPotentialParameters, potentials_for_bonds  , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies         , species               , INPUT, REQUIRED );
  public:
    inline void execute() override final
    {

      // Convert the strings automatically read from the input file containing the force field parameters
      // in potentials_for_bonds->m_bond_desc by a YAML routine
      // into the unordered map potentials_for_bonds->m_type_pair_to_potential
      // This map contains the same info but is more quickly accessed during the simulation.
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

    }

  };

  class ReadBendPotentials : public OperatorNode
  {

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( BendsPotentialParameters , potentials_for_angles , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies          , species               , INPUT, REQUIRED );

  public:
    inline void execute ()  override final
    {

      // Convert the strings automatically read from the input file containing the force field parameters
      // in potentials_for_angles->m_potentials by a YAML routine
      // into the unordered map potentials_for_angles->m_type_to_potential
      // This map contains the same info but is more quickly accessed during the simulation.
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

    }

  };

  class ReadTorsionPotentials : public OperatorNode
  {
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( TorsionsPotentialParameters , potentials_for_torsions , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies             , species                 , INPUT, REQUIRED );

    inline void execute ()  override final
    {

      ParticleSpecies& species = *(this->species);

      // Convert the strings automatically read from the input file containing the force field parameters
      // in potentials_for_torsions->m_potentials by a YAML routine
      // into the unordered map potentials_for_torsions->m_type_to_potential
      // This map contains the same info but is more quickly accessed during the simulation.
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
    }
  };

  class ReadImproperPotentials : public OperatorNode
  {
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies              , species                  , INPUT, REQUIRED );

    inline void execute ()  override final
    {

      ParticleSpecies& species = *(this->species);

      // Convert the strings automatically read from the input file containing the force field parameters
      // in potentials_for_impropers->m_potentials by a YAML routine
      // into the unordered map potentials_for_impropers->m_type_to_potential
      // This map contains the same info but is more quickly accessed during the simulation.


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
        }
      }
    }
  };


  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "read_bond_potentials", make_simple_operator< ReadBondPotentials > );
    OperatorNodeFactory::instance()->register_factory( "read_bend_potentials", make_simple_operator< ReadBendPotentials > );
    OperatorNodeFactory::instance()->register_factory( "read_torsion_potentials", make_simple_operator< ReadTorsionPotentials > );
    OperatorNodeFactory::instance()->register_factory( "read_improper_potentials", make_simple_operator< ReadImproperPotentials > );
  }


}

