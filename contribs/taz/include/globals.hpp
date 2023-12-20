/// @file 
/// @brief Header to include all global variables in the code

#ifndef __GLOBALS_HPP_INCLUDED
#define __GLOBALS_HPP_INCLUDED


#include <cstdlib>


class GlobalInfo;
class ForceField;
class ReferenceMap;

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

namespace Global {

  extern int masterNode;
  extern uint maxNumberOfThreads;
  extern int seed;
  extern ReferenceMap reference;
  extern GlobalInfo domainInfo;
	extern ForceField* ffield;
	extern void terminate();

#if __use_orchestrator
#if __use_tbb
  extern tbb::task_arena* simulationArena;
  extern tbb::task_arena* analyticsArena;
  extern tbb::task_group* simulationGroup;
  extern tbb::task_group* analyticsGroup;
#endif /* __use_tbb */
#endif /* __use_orchestrator */

#if __debug_tasks
  extern tbb::enumerable_thread_specific<std::pair<double, double> > threadTiming;
#endif /* __debug_tasks */

}

#endif // __GLOBALS_HPP_INCLUDED
