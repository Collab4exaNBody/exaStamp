/// @file
/// @brief Definition of the functions that handle OpenMP threads

#include <atomic>
#include <vector>

#include "omp.h"

#include <mutex>

#if __use_lib_hwl

#include "hwloc.h"

#endif



typedef int affinity_hint;

// 3 mutex options : tbb, standard, openMP or handwritten by (someone)

// Available with openmp thread.
/*#include "tbb/tbb.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/mutex.h"
#include "tbb/spin_mutex.h"
typedef tbb::spin_mutex           spin_mutex;*/

//#define STANDARD

#ifdef STANDARD
typedef std::mutex spin_mutex;
#endif

#define ATOMIC_MUTEX

#ifdef ATOMIC_MUTEX

#include <atomic>
#include <thread>

class spin_mutex {
private:
    std::atomic_flag _lock = ATOMIC_FLAG_INIT;
    std::atomic<std::size_t> _spin_pred{0};

public:

    spin_mutex() {}
  
    ~spin_mutex() {}

    inline bool try_lock() {
        return !_lock.test_and_set(std::memory_order_acquire);
    }

    inline void lock() {
        std::size_t spin_count{0};

        while (!try_lock()) { 
  
            if (spin_count < _spin_pred * 2)
                continue;

            std::this_thread::sleep_for(std::chrono::microseconds(1));
            _spin_pred += (spin_count - _spin_pred) / 8;
          }

    }

    inline void unlock() {
        _lock.clear(std::memory_order_release);
    }
};
#endif


//#define OMP_MUTEX

#ifdef OMP_MUTEX
/// @brief Dummy spin_mutex class
/// @version OMP mutex 
class spin_mutex : public omp_lock_t {
public: 

  spin_mutex() {omp_init_lock(this);}
  
  ~spin_mutex() {omp_destroy_lock(this);}
  

  inline void lock    () {omp_set_lock(this);}

  inline void unlock  () {omp_unset_lock(this);}

  inline bool try_lock() { return omp_test_lock(this); }

};
#endif


/// @brief Provides thread local storage (TLS) for elements of type T
///
/// There is no thread so this just provide a T
/// @version OpenMP
template <class T> class Th_ {

public:

  /// @brief default constructor
  explicit Th_() { 
  
    m_t.resize(omp_get_num_procs());

    #pragma omp parallel for
    for(int i = 0; i < omp_get_num_procs(); ++i)
    {
      m_t[i] = new T;
    }
  }

	/// @brief Return the local version
	/// @return Local and only version of T
  inline T& local() { return (*m_t[omp_get_thread_num()]); }

private:

  std::vector<T*> m_t; ///< Local and only storage

};


/// @brief Accessor to default number of threads
/// @version OpenMP
/// @brief Accessor to default number of threads
/// @version OpenMP
inline int default_num_threads() {
  return omp_get_max_threads();
}


/// @brief Creation of a parallel region (simplest case)
/// @version OpenMP
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end parameter
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] lambda Function to apply
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda)
{
    #pragma omp parallel for schedule(runtime)
    for(auto i = begin; i < end; ++i)
       lambda(i,i+1); 

}



/// @brief Creation of a parallel region (with specified grain)
/// @version OpenMP
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end and grain parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] grain Grain parameter
/// @param [in] lambda Function to apply
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda) {
  
  /*  int end_t;
    #pragma omp parallel for private(end_t) schedule(runtime)
    for(int i = begin; i < end; i+=grain)
    {
       end_t=i+grain;
       end_t = end_t>end ? end : end_t;
       lambda(i,end_t);  
    }*/

     #pragma omp parallel for schedule(runtime)
    for(auto i = begin; i < end; ++i)
       lambda(i,i+1); 
  
}


/// @brief Creation of a parallel region (with specified partitioner)
/// @version OpenMP
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end parameter
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] lambda Function to apply
/// @param [in] partitioner Partitioner
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda, affinity_hint& partitioner) {

     #pragma omp parallel for schedule(runtime)
    for(auto i = begin; i < end; ++i)
       lambda(i,i+1); 
 

} 



/// @brief Creation of a parallel region (with specified grain and partitioner)
/// @version OpenMP
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end and grain parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] grain Grain parameter
/// @param [in] lambda Function to apply
/// @param [in] partitioner Partitioner
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda, affinity_hint& partitioner) {
  
   /* int end_t;
    #pragma omp parallel for private(end_t) schedule(runtime)
    for(int i = begin; i < end; i+=grain)
    {
       end_t=i+grain;
       end_t = end_t>end ? end : end_t;
       lambda(i,end_t);  
    }*/

     #pragma omp parallel for schedule(runtime)
    for(auto i = begin; i < end; ++i)
       lambda(i,i+1); 
}

 
/// @brief Class which provides an array of mutexes
/// @version OpenMP
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
    if (n>m_size) {
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



///// WARNING MAYBE FALSE
///// DON'T USE OpenMP with MOLECULES

/// @brief Class which provide a read/write mutex
class MrwMutex : public spin_mutex {

public :

	/// @brief Switch from read mode to write mode
	///
	///
	inline void switch_to_write() {
		unlock();
		lock();
	}

	/// @brief Switch from write mode to read mode
	inline void switch_to_read() {
		unlock();
		lock();
	}
	
		/// @brief Switch from write mode to read mode
	inline void lock_read() {
		lock();
	}

};


/// @brief Class for the thread observer
/// @version OpenMP, default (empty containers)
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



