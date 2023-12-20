#if __use_orchestrator

/// @file
/// @brief Header to include shared variables between the threads

#ifndef __SHARED_HPP_INCLUDED
#define __SHARED_HPP_INCLUDED


#include <atomic>

#include "io/particleInSitu.hpp"


/// @brief Structure to share data between master threads
///
/// The goal of this structure is to hold data that need to be shared between the simulation master thread and the orchestrator
struct SharedVariables
{
  
  std::atomic<bool> dataReady;         ///< boolean to tell when simulation wrote the data for analytics
  std::atomic<bool> simulationWake;    ///< boolean to tell if the simulation can resume
  std::atomic<bool> stopOrchestrator;  ///< boolean to tell when the simulation is over
  std::atomic<bool> dataRefreshable;   ///< boolean to tell the data can be refreshed by simulation
  bool isAsynchronous;                 ///< boolean to tell whether the analytics are asynchronous
  bool copyGhost;

  ParticleInSitu* particles;           ///< structure to hold the particles attributes that are copied
  
  SharedVariables() :
    dataReady(false),
    simulationWake(false),
    stopOrchestrator(false),
    dataRefreshable(true), // for the first output iteration
    isAsynchronous(false),
    copyGhost(true)
  {}

  ~SharedVariables() { }
  
};

#endif /* __SHARED_HPP_INCLUDED */

#endif /* __use_orchestrator */
