/// @file 
/// @brief Implementations of LangevinSplitting


#include "domain/domainInterface.hpp"

#include "time/timeIntegration.hpp"


/// @brief Do one step of the Langevin Splitting integration scheme for a domain
/// @param [in,out] domain Domain where the step is done
/// @param [in] time Time step size
void LangevinSplitting::oneStep(DomainInterface* domain, double time) {

  domain->pushPositions2ndOrder(time);
  domain->pushVelocities1stOrder(0.5*time);

  forceComputer->doComputeForces(domain);
  
  // Using the method of Verlet lists, neighbour lists are not cleared at each iteration.
  if(!Global::reference.isVerlet())
    forceComputer->clearNeighborLists(domain);

  domain->pushVelocities1stOrder(0.5*time);

  domain->pushDissipationLangevin(time, gamma, beta);

}


