/// @file 
/// @brief Definition of all global variables


#include <map>
#include <typeindex>

#include "domain/domainInfo.hpp"

#include "forceField/forceField.hpp"

#include "parallel/mympi.hpp"

#include "referenceMap.hpp"

#if __use_orchestrator
#if __use_tbb
#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#endif /* __use_tbb */
#endif /* __use_orchestrator */

#if __debug_tasks
#include "tbb/enumerable_thread_specific.h"
#include <utility>
#include <chrono>
#endif /* __debug_tasks */

// Note : DO NOT INCLUDE globals.hpp


class ForceField;


/// @brief Namespace containing all global variables in the code
namespace Global {

  /// @brief Global variable : rank of the "Master" Node
  int masterNode = 0;

  /// @brief Global variable : max number of threads
  uint maxNumberOfThreads = 0;

  /// @brief Global variable : seed for random generator initialization
  int seed = 0;
  
  /// @brief Global variable to store particles types and potentials
  ReferenceMap reference = ReferenceMap();

  /// @brief Global variable to store everything about the domain
  GlobalInfo domainInfo = GlobalInfo();


#if __use_orchestrator
#if __use_tbb
  /// @brief Global variables to store objects dedicated to arenas
  tbb::task_arena* simulationArena;  ///< Arena dedicated to simulation tasks
  tbb::task_arena* analyticsArena;   ///< Arena dedicated to analytics tasks
  tbb::task_group* simulationGroup;  ///< Group of threads for the simulation
  tbb::task_group* analyticsGroup;   ///< Group of threads for the analytics
#endif /* __use_tbb */
#endif /* __use_orchestrator */

#if __debug_tasks
  tbb::enumerable_thread_specific<std::pair<double, double> > threadTiming(std::make_pair(0.0, 0.0)); ///< TLS object to measure the time spent by each thread to tasks execution
#endif /* __debug_tasks */

  /// @brief Global pointer to the force field
	ForceField* ffield;

	/// @brief Things to do at the end of the simulation to end the program correctly
	void terminate() {
		delete ffield;
	}

}


/// @brief Map to store MPI datatypes
std::map<std::type_index, MPI_Datatype> MPI__TYPES_MAP = std::map<std::type_index, MPI_Datatype>();
