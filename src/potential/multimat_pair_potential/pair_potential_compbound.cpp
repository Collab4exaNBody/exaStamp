#include <exaStamp/potential_factory/pair_potential_factory.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <onika/physics/units.h>
#include <onika/cpp_utils.h>

#include <yaml-cpp/yaml.h>

#undef USTAMP_POTENTIAL_WITH_VIRIAL
#include "pair_potential_compbound_operator.hxx"

#define USTAMP_POTENTIAL_WITH_VIRIAL 1
#include "pair_potential_compbound_operator.hxx"

namespace exaStamp
{
  using namespace exanb;

  class PairPotentialCompbound : public PairPotential
  {
  public:

    // constructor
    PairPotentialCompbound(const std::vector< std::pair< double , std::shared_ptr<PairPotential> > >& pots)
      : m_potentials(pots)
      {}
    
    // factory
    static std::shared_ptr<PairPotential> make_potential(YAML::Node node)
    {
      using onika::physics::Quantity;
      
      std::vector< std::pair< double , std::shared_ptr<PairPotential> > > potentials;

      // here, node refers to the "parameters" key
      // in our case, this should be a list of potential descriptors
      if( ! node.IsSequence() )
      {
        lerr << "'parameters' value of a compbound potential should be a list" << std::endl;
        return nullptr;
      }
      
      for( auto pot : node )
      {
        if( ! pot.IsMap() )
        {
          lerr << "potential item is not a map as expected" << std::endl;
          return nullptr;
        }
        if( pot["rcut"] && pot["potential"] && pot["parameters"] )
        {
          potentials.push_back( { pot["rcut"].as<Quantity>().convert() , PairPotentialFactory::make_instance(pot) } );
        }
        else
        {
          lerr << "invalid potential description. must contain at least rcut, potential and parameters" << std::endl;
          return nullptr;
        }
      }
      
      return std::make_shared<PairPotentialCompbound>( potentials );
    }

    // specialized for multimaterial, called from traversal method for each atom pair type
    std::shared_ptr<PairPotentialComputeOperator> force_op() override final
    {
      std::vector< std::shared_ptr<PairPotentialComputeOperator> > ops;
      std::vector< double > rcuts;
      for(auto& pot:m_potentials)
      {
        std::shared_ptr<PairPotentialComputeOperator> op = pot.second->force_op();
        op->set_rcut( pot.first );
        ops.push_back( op );
        rcuts.push_back( pot.first );
      }
      return std::make_shared< PairPotentialCompboundOperator >( ops , rcuts );
    }

    std::shared_ptr<PairPotentialComputeVirialOperator> force_virial_op() override final
    {
      std::vector< std::shared_ptr<PairPotentialComputeVirialOperator> > ops;
      std::vector< double > rcuts;
      for(auto& pot:m_potentials)
      {
        std::shared_ptr<PairPotentialComputeVirialOperator> op = pot.second->force_virial_op();
        op->set_rcut( pot.first );
        ops.push_back( op );
        rcuts.push_back( pot.first );
      }
      return std::make_shared< PairPotentialCompboundVirialOperator >( ops , rcuts );
    }

  private:
    std::vector< std::pair< double , std::shared_ptr<PairPotential> > > m_potentials;
  };

  // === register potential factory ===  
  ONIKA_AUTORUN_INIT(pair_potential_compbound)
  {
    PairPotentialFactory::register_factory( "compbound" , PairPotentialCompbound::make_potential );
  }

}


