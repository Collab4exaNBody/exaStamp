/// @file
/// @brief Implementation of the functions that handle TBB threads

#if __use_lib_tbb
#if __use_lib_hwl


#include <iostream>

#include "hwloc.h"

#include "parallel/thread/thread.hpp"


/// @brief Initialize thread_observer
///
///
void thread_observer::init() {

  int err = hwloc_topology_init(& m_topology);
  if (err!=0) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'hwloc_topology_init' : r=" << err 
	     << std::endl;
    return;
  }

  ulong flags = 0;
  hwloc_topology_set_flags(m_topology, flags);
  err = hwloc_topology_load(m_topology);
  if (err!=0) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'hwloc_topology_load' : r=" << err 
	     << std::endl;
    return;
  }
  
  hwloc_bitmap_t origCpuSet = hwloc_bitmap_alloc();
  hwloc_get_cpubind(m_topology, origCpuSet, HWLOC_CPUBIND_THREAD);
  if (err!=0) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'hwloc_get_cpubind' : r=" << err 
	     << std::endl;
    return;
  }
  
  char* result; // TODO : free result
  hwloc_bitmap_asprintf(&result, origCpuSet);
  hwloc_bitmap_free(origCpuSet);

}


/// @brief Bind a thread
/// @param [in] cpu Index of the thread to bind
void thread_observer::bindThread(int cpu) {

  ulong flags = 0;
  hwloc_bitmap_t newCpuSet = hwloc_bitmap_alloc();
  int err = hwloc_get_cpubind(m_topology, newCpuSet, 0 | HWLOC_CPUBIND_THREAD);
  if (err!=0) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'hwloc_get_cpubind' : r=" << err 
	     << std::endl;
    return;
  }
  
  {
    int  cpuIndex = 0;
    uint cpuIter  = 0;
    hwloc_bitmap_foreach_begin(cpuIter, newCpuSet) {
      if (cpuIndex==cpu)
	hwloc_bitmap_only(newCpuSet, cpuIter);
      ++cpuIndex;
    }
    hwloc_bitmap_foreach_end();
  }
  
  hwloc_set_cpubind(m_topology, newCpuSet, flags | HWLOC_CPUBIND_THREAD);
  
  // this is a check:
  
  // hwloc_get_cpubind(m_topology, newCpuSet, 0 | HWLOC_CPUBIND_THREAD);
  // char* result; // TODO : free result
  // hwloc_bitmap_asprintf(&result, newCpuSet);
  // hwloc_bitmap_free(newCpuSet);
  
}

#endif // __use_lib_hwl

#endif // __use_lib_tbb
