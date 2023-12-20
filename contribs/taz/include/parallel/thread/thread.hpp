/// @file 
/// @brief Include for function to handle TBB threads or equivalent empty functions if there is no threads

#ifndef __THREADS_HPP_INCLUDED
#define __THREADS_HPP_INCLUDED


#ifdef __use_lib_omp

#include "parallel/thread/omp.hxx"

#else

#ifdef __use_lib_tbb

#include "parallel/thread/tbb.hxx"

#else

#include "parallel/thread/no_threads.hxx"

#endif

#endif

#endif // __THREADS_HPP_INCLUDED
