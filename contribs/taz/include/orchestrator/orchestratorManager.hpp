#if __use_orchestrator

/// @file
/// @brief Definition of the OrchestratorManager class
///
///

#ifndef __ORCHESTRATOR_HPP_INCLUDED
#define __ORCHESTRATOR_HPP_INCLUDED


#include <mpi.h>

#include "globals.hpp"

#include "io/input.hpp"
#include "io/particleInSitu.hpp"

#include "orchestrator/inSituPlugins.hpp"
#include "orchestrator/shared.hpp"

#include "tbb/tbb.h"


class Node;

/// @brief Class that defines the orchestrator
///
/// WARNING: for the moment the orchestrator works only with tbb!!!!
class OrchestratorManager
{
public:

  /// @brief Default constructor
  /// @param[in] input A reference to the input parameters set by the user
  /// @param[in] shared A reference to the shared variables between simulation and orchestrator
  /// @param[in] node A reference to the node
  OrchestratorManager(Input& input, SharedVariables& shared, Node& node);

  ~OrchestratorManager()
  {
#if __use_tbb
    // Free global arenas memory
    delete Global::analyticsArena;
    delete Global::simulationArena;
    delete Global::analyticsGroup;
    delete Global::simulationGroup;
#endif /* __use_tbb */
    delete availablePlugins;
    delete simulationDataBuffer;
  }

  /// @brief Main function of the orchestrator
  ///
  /// This function defines the behavior of the orchestrator.
  /// It is composed of an infinite loop and a set of actions that need to be executed
  void run();

  /// @brief Launch analytics
  ///
  /// This function launches the analytics on the particles data set for a given step (flow graph case: one analytics at a time with a flow graph)
  /// @param[in] particles A reference to the particles copied by the simulation
  /// @param[in] step The step number
  void launchGraph(ParticleInSitu* particles, int step);

  /// @brief Accessor to the simulation arena size
  int getNbSimulationWorkers(){ return nbSimulationWorkers; }

  /// @brief Accessor to the analytics arena size
  int getNbAnalyticsWorkers(){ return nbAnalyticsWorkers; }

  /// @brief Print orchestrator informations
  void printInfo();

  /// @brief Construct the analytics graph
  void constructGraph();
  
private:

  Node* m_node;                            ///< Pointer to the node
  SharedVariables* m_shared;               ///< Pointer to shared data with simulation master thread
  MPI_Comm analyticsComm;                  ///< Analytics communicator (duplicate of the simulation communicator)
  int nbAnalyticsWorkers;                  ///< Size of the analytics arena
  int nbSimulationWorkers;                 ///< Size of the simulation arena
  Input::SynchronousPolicy syncPolicy;     ///< Synchronous policy
  InSituPluginArray* availablePlugins;     ///< Pointer to the available plugins in the plugin directory
  std::vector<std::string> analyticsList;  ///<  List of the analytics asked by the user
  tbb::flow::graph analyticsGraph;         ///< TBB graph that describes the analytics
  tbb::flow::queue_node<ParticleInSitu*>* simulationDataBuffer; ///< TBB node that will take the simulation data as input
  
};

#endif /* __ORCHESTRATOR_HPP_INCLUDED */

#endif /* __use_orchestrator */
