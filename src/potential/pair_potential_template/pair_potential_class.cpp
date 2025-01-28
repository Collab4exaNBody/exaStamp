#include "pair_potential_template.h"

#include <exaStamp/potential_factory/pair_potential_factory.h>
#include <exaStamp/potential_factory/pair_potential.h>

#undef USTAMP_POTENTIAL_WITH_VIRIAL
#include "pair_potential_force_op.h"

#define USTAMP_POTENTIAL_WITH_VIRIAL 1
#include "pair_potential_force_op.h"

namespace exaStamp
{
  using namespace exanb;


  class USTAMP_POTENTIAL_CLASS final : public PairPotential
  {
  public:

    // constructor
    USTAMP_POTENTIAL_CLASS(const USTAMP_POTENTIAL_PARAMS& p) : m_p(p) {}
    
    // factory
    static std::shared_ptr<PairPotential> make_potential(YAML::Node node)
    {
      return std::make_shared<USTAMP_POTENTIAL_CLASS>( node.as<USTAMP_POTENTIAL_PARAMS>() );
    }

    // specialized for multimaterial, called from traversal method for each atom pair type
    std::shared_ptr<PairPotentialComputeOperator> force_op() override final
    {
      return std::make_shared< USTAMP_POTENTIAL_FORCE_OP_NAME >(m_p);
    }

    std::shared_ptr<PairPotentialComputeVirialOperator> force_virial_op() override final
    {
      return std::make_shared< USTAMP_POTENTIAL_FORCE_VIRIAL_OP_NAME >(m_p);
    }

  private:
    const USTAMP_POTENTIAL_PARAMS m_p;
  };

  // === register potential factory ===  
  ONIKA_AUTORUN_INIT(pair_potential_class)
  {
    PairPotentialFactory::register_factory( USTAMP_POTENTIAL_STRING , USTAMP_POTENTIAL_CLASS::make_potential );
  }

}


