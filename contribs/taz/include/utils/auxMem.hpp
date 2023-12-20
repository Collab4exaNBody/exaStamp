/// @file 
/// @brief A set of basic memory functions

#ifndef __AUX_MEM_HPP_INCLUDED
#define __AUX_MEM_HPP_INCLUDED


#include <sys/resource.h>

#include <cstring>
#include <iostream>
#include <string>


/// @brief Copy a block of memory
///
/// Template wrapper of memcpy
/// @tparam T Class encased in the memory
/// @param [out] to Where to copy
/// @param [in] from What to copy
/// @param [in] n Number of element to copy
template <class T> 
inline void auxMemCpy(T* to, const T* from, const uint n) {
  memcpy(to, from, n*sizeof(T));
}


/// @brief Set values in a block of memory (not used)
///
/// Template wrapper of memset
/// @tparam T Class encased in the memory
/// @param [out] ptr Pointer to the block
/// @param [in] value Value to assign
/// @param [in] n Number of element to set
template <class T> 
inline void auxMemSet(T* ptr, int value, const uint n) {
  memset(ptr, value, n*sizeof(T));
}


/// @brief Copy a block of memory
///
/// Specialization for string (has to be copied by hand)
/// @param [out] to Where to copy
/// @param [in] from What to copy
/// @param [in] n Number of element to copy
template <>
inline void auxMemCpy(std::string* to, const std::string* from, const uint n) { 
  for (uint i=0; i<n; ++i) 
    to[i] = from[i];
}


/// @brief Print a rusage struct [used for debug]
/// @param [in] _rusage Structure to print
inline void print(const rusage* _rusage) {

  // struct timeval ru_utime; /* user time used */
  // struct timeval ru_stime; /* system time used */
 
  std::cout<< "maximum resident set size    : " << _rusage->ru_maxrss   << std::endl;
  std::cout<< "integral shared memory size  : " << _rusage->ru_ixrss    << std::endl;
  std::cout<< "integral unshared data size  : " << _rusage->ru_idrss    << std::endl;
  std::cout<< "integral unshared stack size : " << _rusage->ru_isrss    << std::endl;
  std::cout<< "page reclaims                : " << _rusage->ru_minflt   << std::endl;
  std::cout<< "page faults                  : " << _rusage->ru_majflt   << std::endl;
  std::cout<< "swaps                        : " << _rusage->ru_nswap    << std::endl;
  std::cout<< "block input operations       : " << _rusage->ru_inblock  << std::endl;
  std::cout<< "block output operations      : " << _rusage->ru_oublock  << std::endl;
  std::cout<< "messages sent                : " << _rusage->ru_msgsnd   << std::endl;
  std::cout<< "messages received            : " << _rusage->ru_msgrcv   << std::endl;
  std::cout<< "signals received             : " << _rusage->ru_nsignals << std::endl;
  std::cout<< "voluntary context switches   : " << _rusage->ru_nvcsw    << std::endl;
  std::cout<< "involuntary context switches : " << _rusage->ru_nivcsw   << std::endl;

}

#endif // __AUX_MEM_HPP_INCLUDED
