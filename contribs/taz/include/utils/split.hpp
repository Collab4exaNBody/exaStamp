/// @file 
/// @brief Tool to split work between threads

#ifndef __SPLIT_HPP_INCLUDED
#define __SPLIT_HPP_INCLUDED


/// @brief Split tasks between workers
/// @tparam T Integral type used
/// @param [in] numTasks number of tasks to share
/// @param [in] numWorkers number of workers
/// @param [in] id identifier of the current worker
/// @param [out] start start index for the worker designed by id
/// @param [out] size number of tasks for worker id
template <typename T> inline void split(T numTasks, T numWorkers, T id, T& start, T& size) {

    T q=numTasks/numWorkers;
    T r=numTasks%numWorkers;
    
    if (id<r){
      start = id*(q+1);
      size  = q+1;
    }
    else{
      start = id*q+r;
      size  = q;
    }

}

#endif // __SPLIT_HPP_INCLUDED
