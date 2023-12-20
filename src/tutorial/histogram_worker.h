#pragma once

#include <exanb/core/histogram.h>

#include <memory>
#include <vector>
#include <type_traits>
#include <limits>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <mpi.h>

namespace exaStamp
{

  struct HistogramWorker
  {
    // thread stuff
    std::shared_ptr<std::thread> m_thread;
    std::mutex m;
    std::condition_variable cv;    
    bool ready = false;
    bool processed = false;  

    // consumer input
    std::vector<double> input;
    
    // consumer output
    exanb::Histogram<> output;

    // thread executes this function
    void run( MPI_Comm comm );
  }; 

}

