/// @file
/// @brief Definition of the functions that handle TBB threads

#include <atomic>
#include <vector>

#include "tbb/tbb.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/mutex.h"
#include "tbb/spin_mutex.h"

#if __use_lib_hwl

#include "hwloc.h"

#endif

#if __debug_tasks
#include <chrono>
#endif /* __debug_tasks */
#include "globals.hpp"



/// @brief Provides thread local storage (TLS) for elements of type T
/// @version TBB
/// @tparam T Type to store
template <class T> using Th_ = tbb::enumerable_thread_specific<T, tbb::cache_aligned_allocator<T>, tbb::ets_key_usage_type::ets_key_per_instance>;


/// @brief Shortcut for TBB affinity partitioner :
/// tool to hint the same partitioning in related parallel regions
/// @version TBB
/// @brief Shortcut for TBB affinity partitioner :
/// tool to hint the same partitioning in related parallel regions
/// @version TBB
typedef tbb::affinity_partitioner affinity_hint;
/// @brief Shortcut for TBB spin mutex :
/// mutex that make the thread spin will waiting for the lock to be available,
/// adapted to locks with short instructions
/// @version TBB
typedef tbb::spin_mutex           spin_mutex;
/// @brief Shortcut for TBB spin read/write mutex :
/// mutex that make the thread spin will waiting for the lock to be available,
/// adapted to locks with short instructions,
/// allow multiple lock while reading but only one when writing
/// @version TBB
typedef tbb::spin_rw_mutex        spin_rw_mutex;
/// @brief Shortcut for task_scheduler_init ;
/// tool to meddle in task scheduling
/// @version TBB
typedef tbb::task_scheduler_init  thread_scheduler;


/// @brief Accessor to default number of threads
/// @version TBB
/// @brief Accessor to default number of threads
/// @version TBB
inline int default_num_threads() {
  return tbb::task_scheduler_init::default_num_threads();
}



/// @brief Creation of a parallel region (simplest case)
/// @version TBB
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end parameter
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] lambda Function to apply
#if __use_orchestrator

template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda) {

  Global::simulationArena->execute( [&] {
      Global::simulationGroup->run_and_wait ( [&] {
          tbb::parallel_for( tbb::blocked_range<J>(begin, end), 
                             [&](const tbb::blocked_range<J>& r) {
                               
#if __debug_tasks
                               // Get starting time of the task
                               std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */

                               lambda(r.begin(), r.end());
                               
#if __debug_tasks
                               // Get ending time of the task
                               std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                               tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                               my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
                               
                             }
                             );
        });
    });
}

#else /* __use_orchestrator */

template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda) {
  
  tbb::parallel_for( tbb::blocked_range<J>(begin, end), 
                     [&](const tbb::blocked_range<J>& r) {
#if __debug_tasks
                       // Get starting time of the task
                       std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */
                       
                       lambda(r.begin(), r.end());
                       
#if __debug_tasks
                       // Get ending time of the task
                       std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                       tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                       my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
                       
                     }
                     );
  
}

#endif /* __use_orchestrator */



/// @brief Creation of a parallel region (with specified grain)
/// @version TBB
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end and grain parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] grain Grain parameter
/// @param [in] lambda Function to apply
#if __use_orchestrator

template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda) {

  Global::simulationArena->execute( [&] {
      Global::simulationGroup->run_and_wait ( [&] {
          tbb::parallel_for( tbb::blocked_range<J>(begin, end, grain),
                             [&](const tbb::blocked_range<J>& r) {
                               
#if __debug_tasks
                               // Get starting time of the task
                               std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */

                               lambda(r.begin(), r.end());

#if __debug_tasks
                               // Get ending time of the task
                               std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                               tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                               my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
                               
                             }
                             );
        });
    });
}

#else /* __use_orchestrator */

template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda) {
  
  tbb::parallel_for( tbb::blocked_range<J>(begin, end, grain),
                     [&](const tbb::blocked_range<J>& r) {

#if __debug_tasks
                       // Get starting time of the task
                       std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */
                       
                       lambda(r.begin(), r.end());

#if __debug_tasks
                       // Get ending time of the task
                       std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                       tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                       my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */

                     }
                     );
  
}

#endif /* __use_orchestrator */


/// @brief Creation of a parallel region (with specified partitioner)
/// @version TBB
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end parameter
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] lambda Function to apply
/// @param [in] partitioner Partitioner
#if __use_orchestrator

template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda, affinity_hint& partitioner) {
  
#ifdef __use_tbb_affinity

  Global::simulationArena->execute( [&] {
      Global::simulationGroup->run_and_wait ( [&] {
          tbb::parallel_for( tbb::blocked_range<J>(begin, end), 
                             [&](const tbb::blocked_range<J>& r) {
                               
#if __debug_tasks
                               // Get starting time of the task
                               std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */

                               lambda(r.begin(), r.end());

#if __debug_tasks
                               // Get ending time of the task
                               std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                               tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                               my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
                               
                             },
                             partitioner
                             );
        });
    });

#else /* __use_tbb_affinity */

  Global::simulationArena->execute( [&] {
      Global::simulationGroup->run_and_wait ( [&] {
          tbb::parallel_for( tbb::blocked_range<J>(begin, end), 
                             [&](const tbb::blocked_range<J>& r) {
                               
#if __debug_tasks
                               // Get starting time of the task
                               std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();     
#endif /* __debug_tasks */

                               lambda(r.begin(), r.end());

#if __debug_tasks
                               // Get ending time of the task
                               std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                               tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                               my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
                               
                             }
                             );
        });
    });
  
#endif /* __use_tbb_affinity */

}

#else /* __use_orchestrator */

template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda, affinity_hint& partitioner) {

#ifdef __use_tbb_affinity

  tbb::parallel_for( tbb::blocked_range<J>(begin, end), 
                     [&](const tbb::blocked_range<J>& r) {

#if __debug_tasks
                       // Get starting time of the task
                       std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */
                       
                       lambda(r.begin(), r.end());

#if __debug_tasks
                       // Get ending time of the task
                       std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                       tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                       my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */

                     },
                     partitioner
                     );

#else

  tbb::parallel_for( tbb::blocked_range<J>(begin, end), 
                     [&](const tbb::blocked_range<J>& r) {
                       
#if __debug_tasks
                       // Get starting time of the task
                       std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */
                       
                       lambda(r.begin(), r.end());

#if __debug_tasks
                       // Get ending time of the task
                       std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                       tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                       my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */

                     }
                     );
  
#endif /* __use_tbb_affinity */

}
  
#endif /* __use_orchestrator */



/// @brief Creation of a parallel region (with specified grain and partitioner)
/// @version TBB
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end and grain parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] grain Grain parameter
/// @param [in] lambda Function to apply
/// @param [in] partitioner Partitioner
#if __use_orchestrator
  
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda, affinity_hint& partitioner) {
  
#ifdef __use_tbb_affinity
  
  Global::simulationArena->execute( [&] {
      Global::simulationGroup->run_and_wait ( [&] {
          tbb::parallel_for( tbb::blocked_range<J>(begin, end, grain),
                             [&](const tbb::blocked_range<J>& r) {
                               
#if __debug_tasks
                               // Get starting time of the task
                               std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */

                               lambda(r.begin(), r.end());

#if __debug_tasks
                               // Get ending time of the task
                               std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                               tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                               my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
                               
                             },
                             partitioner
                             );
        });
    });
            
#else /* __use_tbb_affinity */

  Global::simulationArena->execute( [&] {
      Global::simulationGroup->run_and_wait ( [&] {
          tbb::parallel_for( tbb::blocked_range<J>(begin, end, grain),
                             [&](const tbb::blocked_range<J>& r) {
                               
#if __debug_tasks
                               // Get starting time of the task
                               std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */

                               lambda(r.begin(), r.end());

#if __debug_tasks
                               // Get ending time of the task
                               std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                               tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                               my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
                               
                             }
                             );
        });
    });
            
#endif /* __use_tbb_affinity */

}
 
#else /* __use_orchestrator */

template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda, affinity_hint& partitioner) {
  
#ifdef __use_tbb_affinity
  
  tbb::parallel_for( tbb::blocked_range<J>(begin, end, grain),
                     [&](const tbb::blocked_range<J>& r) {

#if __debug_tasks
                       // Get starting time of the task
                       std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */
                       
                       lambda(r.begin(), r.end());

#if __debug_tasks
                       // Get ending time of the task
                       std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                       tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                       my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */

                     },
                     partitioner
                     );
            
#else /* __use_tbb_affinity */

  tbb::parallel_for( tbb::blocked_range<J>(begin, end, grain),
                     [&](const tbb::blocked_range<J>& r) {

#if __debug_tasks
                       // Get starting time of the task
                       std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
#endif /* __debug_tasks */
                       
                       lambda(r.begin(), r.end());

#if __debug_tasks
                       // Get ending time of the task
                       std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
                       tbb::enumerable_thread_specific<std::pair<double, double> >::reference my_timer = Global::threadTiming.local();
                       my_timer.first += std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */

                     }
                     );
            
#endif /* __use_tbb_affinity */

}

 
#endif /* __use_orchestrator */



/// @brief Class which provides an array of mutexes
/// @version TBB
class MMutex {

public:

  /// @brief Default constructor
  MMutex() : m_size(0), m_mtx(nullptr) {}

  /// @brief Destructor (nothing to do)
  ~MMutex() {
    if (m_size>0) delete [] m_mtx;
  }

  /// @brief copy constructor
  MMutex(const MMutex& mtx) : m_size(mtx.m_size), m_mtx(nullptr) {
    if (m_size>0) m_mtx = new spin_mutex [m_size];
  }

  /// @brief Resize the array of mutexes
  /// @param [in] n New size
  inline void check(const uint n) {
    if (n>m_size || 4*n<m_size) {
      m_size=2*n;
      delete [] m_mtx;
      m_mtx = new spin_mutex [m_size];
    }
  }

  /// @brief Lock specified mutex
  /// @param [in] i Index of the mutex
  inline void lock(uint i) {
    m_mtx[i].lock();
  }

  /// @brief Unlock specified mutex
  /// @param [in] i Index of the mutex
  inline void unlock(uint i) {
    m_mtx[i].unlock();
  }

  /// @brief Test if specified mutex is locked
  /// @param [in] i Index of the mutex
  inline bool try_lock(uint i) {
    return m_mtx[i].try_lock();
  }

  /// @brief Accessor to the size
  /// @return size
  inline uint size() const {
  	return m_size;
  }

private:

  // prevent assignment operator by putting private operators
  /// @brief Assignement operator (do not use)
  /// @param [in] mtx Mutex to copy
  MMutex& operator = (const MMutex& mtx);

  uint m_size;       ///< size of the array
  spin_mutex* m_mtx; ///< array of mutexes

};


/// @brief Class which provide a read/write mutex
class MrwMutex : public spin_rw_mutex {

public :

	/// @brief Switch from read mode to write mode
	///
	///
	inline void switch_to_write() {
		unlock();
		lock();
	}

	/// @brief Switch from write mode to read mode
	///
	///
	inline void switch_to_read() {
		unlock();
		lock_read();
	}

};


#if __use_lib_hwl

/// @brief Class for the thread observer
/// @version TBB, lib_hwl
class thread_observer : public tbb::task_scheduler_observer {

public:
  
  /// @brief Default constructor
  thread_observer() : m_threadId(0) {}

  void init();
  void bindThread(int cpu);

  /// @brief Reimplementation of on_scheduler_entry
  /// @param [in] isWorker Not used
  virtual void on_scheduler_entry(bool isWorker) {
    uint i = ++m_threadId;
    bindThread(i-1);
  }

  /// @brief Reimplementation of on_scheduler_exit
  /// @param [in] isWorker Not used
  virtual void on_scheduler_exit(bool isWorker) {
  }

private:

  std::atomic<uint> m_threadId; ///< Thread ID
  hwloc_topology_t  m_topology; ///< Topology

};

#else

/// @brief Class for the thread observer
/// @version TBB, default (empty containers)
class thread_observer {
public:
  
  /// @brief Empty
  void init() {}

  /// @brief Empty
  void bindThread(int) {}

  /// @brief Empty with unused parameters
  void on_scheduler_entry(bool) {}

  /// @brief Empty with unused parameters
  void on_scheduler_exit (bool) {}

  /// @brief Empty with unused parameters
  void observe(bool) {}

};



#endif // __use_lib_hwl

