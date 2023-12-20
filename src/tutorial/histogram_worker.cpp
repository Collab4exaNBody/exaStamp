#include "histogram_worker.h"
#include <iostream>

namespace exaStamp
{

  void HistogramWorker::run( MPI_Comm comm )
  {
    int nprocs = 1;
    int rank = 0;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&rank);

    std::cout << "WORKER: thread is running (rank="<<rank<<", np="<<nprocs<<", thread="<<std::this_thread::get_id()<<")"<<std::endl;

    while(true)
    {
      // Wait until main() sends data
      std::cout << "WORKER: thread is waiting input data to be ready\n";
      std::unique_lock<std::mutex> lk(m);
      cv.wait(lk, [this]{return ready;});
      ready = false;
        
      // after the wait, we own the lock.
      std::cout << "WORKER: thread is processing data\n";

      // ---------------- processing start -------------------
      size_t nsamples = output.m_data.size();
      size_t ndata = input.size();
              
      // min max computation
      double min_val = std::numeric_limits<double>::max();
      double max_val = std::numeric_limits<double>::lowest();
      for(size_t i=0;i<ndata;i++)
      {
        min_val = std::min( min_val , input[i] );
        max_val = std::max( max_val , input[i] );
      }

      // MPI min/max
      {
        double tmp[2] = { -min_val , max_val };
        MPI_Allreduce(MPI_IN_PLACE,tmp,2, MPI_DOUBLE , MPI_MAX , comm);
        min_val = - tmp[0];
        max_val = tmp[1];
      }

      size_t hist_size = nsamples;
      output.m_min_val = min_val;
      output.m_max_val = max_val;
      output.m_data.resize( hist_size );

      for(size_t i=0;i<ndata;i++)
      {
        double x = input[i];
        ssize_t bin = static_cast<ssize_t>( ( (x-min_val) * hist_size ) / ( max_val - min_val ) );
        if( bin < 0 ) { bin=0; }
        if( bin >= static_cast<ssize_t>(hist_size) ) { bin = hist_size-1; }
        output.m_data[bin] += 1.0;
      }

      MPI_Allreduce(MPI_IN_PLACE,output.m_data.data(),hist_size,MPI_DOUBLE,MPI_SUM,comm);
      // ---------------- processing end -------------------

      // Send data back to main()
      processed = true;
      std::cout << "WORKER: thread signals data processing completed\n";
   
      // Manual unlocking is done before notifying, to avoid waking up
      // the waiting thread only to block again (see notify_one for details)
      // TODO unlock and notify
      // lk. ...
      // cv. ...
    }
  }

}

