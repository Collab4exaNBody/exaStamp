#include <exaStamp/potential_factory/pair_potential.h>

#include <onika/log.h>

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  std::shared_ptr<PairPotentialComputeOperator> PairPotential::force_op()
  {
    lerr << "potential does not have an operator for force/no-symetry/no-virial computation" <<std::endl;
    std::abort();
  }
  
  std::shared_ptr<PairPotentialComputeVirialOperator> PairPotential::force_virial_op()
  {
    lerr << "potential does not have an operator for force/no-symetry/virial computation" <<std::endl;
    std::abort();
  }

}


