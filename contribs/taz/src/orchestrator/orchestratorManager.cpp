#if __use_orchestrator

/// @file
/// @brief Implementation of the OrchestratorManager class


#include <thread>
#include <vector>
#include <dlfcn.h>
#include "io/particleInSitu.hpp"
#include "io/particleOutputDAT.hpp"
#include "io/particleOutputVTK.hpp"

#include "orchestrator/shared.hpp"
#include "orchestrator/orchestratorManager.hpp"

#include "tbb/tbb.h"
#include "tbb/task_arena.h"

#include "parallel/node.hpp"
#include "parallel/thread/tbb.hxx"


OrchestratorManager::OrchestratorManager(Input& input, SharedVariables& shared, Node& node)
{
  // Initialize member parameters
  this->syncPolicy = input.syncPolicy;
  this->m_node = &node;
  this->m_shared = &shared;
  for (uint i=0; i<input.analyticsGraph.size(); ++i)
    this->analyticsList.push_back(input.analyticsGraph[i]);
  std::string repname = __xsp_compute_rep;
  this->availablePlugins = new InSituPluginArray(repname + "/deps/analytics/plugins");
  this->simulationDataBuffer = new tbb::flow::queue_node<ParticleInSitu*>(this->analyticsGraph);

  // Duplicate the simulation communicator for analytics
  MPI_Comm_dup(MPI_COMM_WORLD, &this->analyticsComm); 

  // Initialize isAsynchronous
  if (this->syncPolicy == Input::ASYNCHRONOUS_WITHOUT_ARENAS || this->syncPolicy == Input::ASYNCHRONOUS_WITH_ARENAS)
    this->m_shared->isAsynchronous = true;

  // Case when arena sizes are not provided by the user
  if (this->syncPolicy == Input::SYNCHRONOUS_WITHOUT_ARENAS || this->syncPolicy == Input::ASYNCHRONOUS_WITHOUT_ARENAS)
    {
      nbSimulationWorkers = Global::maxNumberOfThreads;
      nbAnalyticsWorkers = Global::maxNumberOfThreads;
    }
  // Case when arena sizes are provided by the user
  else if (this->syncPolicy == Input::SYNCHRONOUS_WITH_ARENAS || this->syncPolicy == Input::ASYNCHRONOUS_WITH_ARENAS)
    {
      nbAnalyticsWorkers = input.maxNbAnalyticsWorkers;
      nbSimulationWorkers = input.maxNbSimulationWorkers;
    }
  // Case when no synchronous policy is set
  else
    {
      std::cerr << "Please provide a synchronous policy" << std::endl;
      exit(-1);
    }
  
  // Allocate arenas and groups according to sizes
#if __use_tbb

  Global::analyticsArena = new tbb::task_arena(nbAnalyticsWorkers);
  Global::simulationGroup = new tbb::task_group();
  Global::analyticsGroup = new tbb::task_group();
  Global::simulationArena = new tbb::task_arena(nbSimulationWorkers); // no extra parameter means that we will have 1 master thread joining the arena and hence TBB will create N-1 worker threads
  
#endif /* __use_tbb */
    
}

void OrchestratorManager::constructGraph()
{
  
  int nbNodes=this->analyticsList.size(); // number of function_nodes in the graph 
  std::vector<tbb::flow::function_node<ParticleInSitu*, ParticleInSitu*>* > analyticsNodes; // vector of pointers to function_nodes to store the nodes of the graph
  analyticsNodes.reserve(nbNodes);

  // Create the nodes based on the asked analytics
  for (int a=0; a<nbNodes; ++a)
    {
      analyticsNodes[a] = new tbb::flow::function_node<ParticleInSitu*, ParticleInSitu*>
        (this->analyticsGraph,
         tbb::flow::unlimited,
         [this, a] (ParticleInSitu* particles) -> ParticleInSitu*
        {

          InSituPlugin* plug = this->availablePlugins->getByName(this->analyticsList[a]);
          plug->IP_run(particles, this->analyticsComm);
              
          return particles;
          
        }
         );
    }
    
  make_edge(*this->simulationDataBuffer, *analyticsNodes[0]);

  for (uint a=0; a<analyticsList.size()-1; ++a)
    make_edge(*analyticsNodes[a], *analyticsNodes[a+1]);
  
}

void OrchestratorManager::run()
{
  
  // Call TBB task_scheduler to instantiate the orchestrator as a TBB master thread
  thread_scheduler scheduler;

  // Create TBB Flow graph
  constructGraph();
            
  while (1)
    {

      this->m_node->getMetrics()->tic(Metrics::ORCHESTRATORSLEEP); 
      
      // Wait for data to be available
      while (!this->m_shared->dataReady)
        {
          std::this_thread::sleep_for(std::chrono::microseconds(10)); // arbitrary time for active wait
          
          // Condition to break the infinite loop: simulation master thread sets the stopOrchestrator flag
          if (this->m_shared->stopOrchestrator)
            {
              this->m_node->getMetrics()->toc(Metrics::ORCHESTRATORSLEEP);
              this->m_node->getMetrics()->store(Metrics::ORCHESTRATORSLEEP);
              return;
            }
          
        }
      
      this->m_node->getMetrics()->toc(Metrics::ORCHESTRATORSLEEP);
      this->m_node->getMetrics()->store(Metrics::ORCHESTRATORSLEEP);
      
      this->m_shared->dataReady = false;

      this->m_node->getMetrics()->tic(Metrics::ANALYTICS);

      // Launch graph execution
      Global::analyticsArena->execute( [&] {
          this->simulationDataBuffer->try_put(this->m_shared->particles);
        });

      // Wait for graph execution
      Global::analyticsArena->execute( [&] {
          this->analyticsGraph.wait_for_all();
        });
          
      this->m_node->getMetrics()->toc(Metrics::ANALYTICS);
      this->m_node->getMetrics()->store(Metrics::ANALYTICS);
      
      if (this->syncPolicy == Input::SYNCHRONOUS_WITHOUT_ARENAS || this->syncPolicy == Input::SYNCHRONOUS_WITH_ARENAS)
        {     
          // Unblock simulation
          this->m_shared->simulationWake = true; 
        }
      
      delete this->m_shared->particles; // clean particles data structure after analytics are over
      this->m_shared->dataRefreshable = true; // tell the simulation it can overwrite the data
     
    }
  
}

void OrchestratorManager::printInfo()
{
  std::cout << "ORCHESTRATOR INFO" << std::endl;
  std::cout << "  Synchronous policy      : ";
  if (this->syncPolicy == Input::SYNCHRONOUS_WITH_ARENAS || this->syncPolicy == Input::SYNCHRONOUS_WITHOUT_ARENAS)
    std::cout << "synchronous" << std::endl;
  else if (this->syncPolicy == Input::ASYNCHRONOUS_WITH_ARENAS || this->syncPolicy == Input::ASYNCHRONOUS_WITHOUT_ARENAS)
    std::cout << "asynchronous" << std::endl;
  std::cout << "  Simulation threads      : " << this->nbSimulationWorkers << std::endl;
  std::cout << "  Analytics threads       : " << this->nbAnalyticsWorkers << std::endl;
  std::cout << "  Analytics               : ";
  for (uint a=0; a<this->analyticsList.size(); ++a)
    std::cout << this->analyticsList[a] << " ";
  std::cout << std::endl;

  std::cout << std::endl;
}

#endif /* __use_orchestrator */
