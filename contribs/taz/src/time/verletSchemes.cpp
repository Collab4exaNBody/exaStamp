/// @file 
/// @brief Implementations of NewtonVerletLeapfrog, NewtonVerletVelocity and NewtonVerlet algorithms


#include "domain/domainInterface.hpp"

#include "time/timeIntegration.hpp"


/// @brief Do one step of the Leapfrog integration scheme for a domain
/// @param [in,out] domain Domain where the step is done
/// @param [in] time Time step size
void NewtonVerletLeapfrog::oneStep(DomainInterface* domain, double time) {

  domain->pushPositions1stOrder(0.5*time);

  forceComputer->doComputeForces(domain);
  forceComputer->clearNeighborLists(domain);

  domain->pushVelocities1stOrder(time);
  domain->pushPositions1stOrder(0.5*time);
  
}

/// @brief Do one step of the Verlet Velocity integration scheme
/// @param [in,out] domain Domain where the step is done
/// @param [in] time Time step size
void NewtonVerletVelocity::oneStep(DomainInterface* domain, double time) {

  domain->pushPositions2ndOrder(time);
  domain->pushVelocities1stOrder(0.5*time);

  forceComputer->doComputeForces(domain);
  
  // Using the method of Verlet lists, neighbour lists are not cleared at each iteration.
  if(!Global::reference.isVerlet())
    forceComputer->clearNeighborLists(domain);

  domain->pushVelocities1stOrder(0.5*time);

}
