#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/histogram.h>
#include <exanb/mpi/data_types.h>

#include <mpi.h>

#include "histogram_worker.h"

namespace exaStamp
{

  struct TutorialHistoThreadFinish : public OperatorNode
  {      
    ADD_SLOT( HistogramWorker, insitu_worker , INPUT_OUTPUT );
    ADD_SLOT( Histogram<> , histogram , OUTPUT );

    inline void execute () override final
    {      
      // wait for the worker
      {
          std::unique_lock<std::mutex> lk(insitu_worker->m);
          insitu_worker->cv.wait(lk, [this]{return insitu_worker->processed;});
          insitu_worker->processed = false;
      }
      std::cout << "MAIN: Worker done" << std::endl;
      
      histogram->m_min_val = insitu_worker->output.m_min_val;
      histogram->m_max_val = insitu_worker->output.m_max_val;
      histogram->m_data = std::move( insitu_worker->output.m_data );
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "histoseq_thread_finish" , make_simple_operator< TutorialHistoThreadFinish > );
  }

}

