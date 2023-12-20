/// @file
/// @brief Empty containers to handle no thread case

/// @brief Dummy affinity_hint shortcut
/// @version no_threads
typedef int affinity_hint;


/// @brief Provides thread local storage (TLS) for elements of type T
///
/// There is no thread so this just provide a T
/// @version no_threads
template <class T> class Th_ {

public:

	/// @brief Return the local version
	/// @return Local and only version of T
  inline T& local() { return m_t; }

private:

  T m_t; ///< Local and only storage

};


/// @brief Dummy spin_mutex class
/// @version no_threads
class spin_mutex {
public: 

	/// @brief Empty
  inline void lock    () {}
	/// @brief Empty
  inline void unlock  () {}
	/// @brief Empty with meaningless return
  inline bool try_lock() { return true; }

};


/// @brief Dummy spin_rw_mutex class
/// @version no_threads
class MrwMutex {

public:

	/// @brief Empty
  inline void lock    () {}
	/// @brief Empty
  inline void lock_read    () {}
	/// @brief Empty
  inline void unlock  () {}
	/// @brief Empty
  inline void switch_to_write  () {}
	/// @brief Empty
  inline void switch_to_read  () {}
	/// @brief Empty with meaningless return
  inline bool try_lock() { return true; }
	/// @brief Empty with meaningless return
  inline bool try_lock_read() { return true; }

};


/// @brief Dummy thread_scheduler class
/// @version no_threads
class thread_scheduler {
public:

	/// @brief Empty with unused parameters
  thread_scheduler(const uint) {}

};


/// @brief Return default number of threads if no parallelization
/// @version no_threads
inline int default_num_threads() {
  return 1;
}


/// @brief Dummy parallel_region (simplest case) : just launch the function
/// @version no_threads
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] lambda Function to apply
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda) {
  lambda(begin, end);
}


/// @brief Dummy parallel_region (with specified grain) : just launch the function
/// @version no_threads
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end and grain parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] grain Grain parameter (not used)
/// @param [in] lambda Function to apply
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda) {
  lambda(begin, end);
}


/// @brief Dummy parallel_region (with specified partitioner) : just launch the function
/// @version no_threads
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] lambda Function to apply
/// @param [in] partitioner Partitioner (not used)
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, Lambda lambda, affinity_hint& partitioner) {
  lambda(begin, end);
}

/// @brief Dummy parallel_region (with specified grain and partitioner) : just launch the function
/// @version no_threads
/// @tparam I Type of the begin parameter
/// @tparam J Type of the end and grain parameters
/// @tparam Lambda Type of the function to apply
/// @param [in] begin Begin parameter
/// @param [in] end End parameter
/// @param [in] grain Grain parameter (not used)
/// @param [in] lambda Function to apply
/// @param [in] partitioner Partitioner (not used)
template <class I, class J, typename Lambda>
inline void parallel_region(const I& begin, const J& end, const J& grain, Lambda lambda, affinity_hint& partitioner) {
  lambda(begin, end);
}


/// @brief Dummy MMutex class
/// @version no_threads
class MMutex {

public:

	/// @brief Default constructor
  MMutex() {}

  /// @brief Destructor (nothing to do)
  ~MMutex() {}

	/// @brief Empty
  inline void check   (uint) {}
	/// @brief Empty
  inline void lock    (uint) {}
	/// @brief Empty
  inline void unlock  (uint) {}
	/// @brief Empty with meaningless return
  inline bool try_lock(uint) { return true; }
  /// @brief Empty accessor to the size
  inline uint size() const {return 1;}
};


/// @brief Dummy thread_observer class
/// @version no_threads
class thread_observer {
public:
  
	/// @brief Empty
  void init() {}
  /// @brief Empty with unused parameters
  void bindThread(int) {}

  /// @brief Empty with unused parameters
  void on_scheduler_entry(bool) {}
  /// @brief Empty with unused parameters
  void on_scheduler_exit (bool) {}

  /// @brief Empty with unused parameters
  void observe(bool) {}

};
