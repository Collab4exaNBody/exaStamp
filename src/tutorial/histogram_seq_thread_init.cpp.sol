#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/histogram.h>

#include "histogram_worker.h"

#include <mpi.h>

namespace exaStamp
{

  struct TutorialHistoThreadInit : public OperatorNode
  {      
    ADD_SLOT( MPI_Comm        , mpi_dup       , INPUT , REQUIRED);
    ADD_SLOT( HistogramWorker , insitu_worker , OUTPUT );

    inline void execute () override final
    {
      // create and start worker thread
      insitu_worker->m_thread = std::make_shared<std::thread> ( &HistogramWorker::run , & (*insitu_worker) , *mpi_dup );
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "histoseq_thread_init" , make_simple_operator< TutorialHistoThreadInit > );
  }

}

