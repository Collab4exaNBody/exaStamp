/// @file
/// @brief Implementation of force computation strategies

#ifndef __COMPUTE_FORCE_STRATEGY_HPP_INCLUDED
#define __COMPUTE_FORCE_STRATEGY_HPP_INCLUDED


#include "domain/domainInterface.hpp"

#include <parallel/metrics.hpp>

extern Metrics * ptM; ///< pointer on metrics (in node class)
extern bool updateVerletLists; ///< boolean value to need if the Verlet list will be updated ot not

/*!
  \def __toc(x)
  get a time point and store this value in the metrics structure.
*/
#define __toc(x) { ptM->toc(Metrics::x); \
                 ptM->store(Metrics::x);}

/// @brief Strategy Pattern to compute force with or without communication 
/// recovering within a Domain
class ComputeForceStrategy {

public:

  /// @brief Default constructor
  ComputeForceStrategy() {}
  /// @brief Destructor (nothing to do)
  virtual ~ComputeForceStrategy() {}

  /// @brief Compute forces with only cell lists
  virtual void doComputeForces(DomainInterface* domain) = 0;
  
  /// @brief Clear neighbor lists
  virtual void clearNeighborLists(DomainInterface* domain) = 0;
};


/// @brief Strategy pattern to compute force with communication recovering within a Domain
class ComputeForceRecovering : public ComputeForceStrategy {

public:

  /// @brief Constructor
  ComputeForceRecovering() {}

  /// @brief Destructor (nothing to do)
  virtual ~ComputeForceRecovering() {}

  /// @brief Compute forces with communication recovering
  virtual void doComputeForces(DomainInterface* domain) {
  
    domain->clearGhost();
    
      if(updateVerletLists)
  	  {
        domain->updateCellsInside();
        domain->refineCells();
      }
      domain->updateGhost();

      if(updateVerletLists)
      {
        clearNeighborLists(domain);
        domain->makeNeighborListsInside();
      }

      domain->doComputeForcesInsideVerlet();
      domain->collectGhost();
      
      if(updateVerletLists)
      {
        domain->makeNeighborListsOnEdges();
        updateVerletLists=false;
      }
      domain->doComputeForcesOnEdgesVerlet  ();

  }

  /// @brief Clear neighbor lists
  virtual void clearNeighborLists(DomainInterface* domain) {

    domain->clearNeighborLists();
    
  }

};




/// @brief Strategy pattern to compute force without communication recovering within a Domain
class ComputeForceNoRecovering : public ComputeForceStrategy {

public:

  /// @brief Constructor
  ComputeForceNoRecovering() {}

  /// @brief Destructor (nothing to do)
  virtual ~ComputeForceNoRecovering() {}

  /// @brief Compute forces with communication recovering
  virtual void doComputeForces(DomainInterface* domain) {
  
    // Boolean retrieving (verlet list method is used) 
    // if neighbour lists are updated during this iteration. 
    bool updateParticlesAndNeighbors = updateVerletLists;
    
    domain->clearGhost();

    if(updateParticlesAndNeighbors)
    {
      domain->updateCellsInside();

      if(ptM->data.isTimerRefine())
        ptM->tic(Metrics::REFINE);
        
      domain->refineCells();

      if(ptM->data.isTimerRefine())
      {
        ptM->toc(Metrics::REFINE);
        ptM->store(Metrics::REFINE);
      }
    }

    /* timer for ghost part */
    if(ptM->data.isTimerGhost())
      ptM->tic(Metrics::GHOST);
    
    domain->updateGhost();
    domain->collectGhost();

    if(ptM->data.isTimerGhost())
    {
      ptM->toc(Metrics::GHOST);
      ptM->store(Metrics::GHOST);
    }
    
    /* timer for neighbours part */
     if(ptM->data.isTimerNeighbours())
      ptM->tic(Metrics::NEIGHBOURS);


    if(updateParticlesAndNeighbors)
    {
      clearNeighborLists(domain);
      domain->makeNeighborLists();
      updateVerletLists=false;
    }  

    if(ptM->data.isTimerNeighbours())
    {
      ptM->toc(Metrics::NEIGHBOURS);
      ptM->store(Metrics::NEIGHBOURS);
    } 

    /* timer for potential part */
    if(ptM->data.isTimerPotential())
      ptM->tic(Metrics::POTENTIAL);
    
    domain->doComputeForcesVerlet();

    if(ptM->data.isTimerPotential())
    {
      ptM->toc(Metrics::POTENTIAL);
      ptM->store(Metrics::POTENTIAL);
    }
  }


  /// @brief Clear neighbor lists
  virtual void clearNeighborLists(DomainInterface* domain) {
    
    domain->clearNeighborLists();
    
  }

};

#endif // __COMPUTE_FORCE_STRATEGY_HPP_INCLUDED
