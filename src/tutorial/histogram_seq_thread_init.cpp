#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/histogram.h>

#include "histogram_worker.h"

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct TutorialHistoThreadInit : public OperatorNode
  {      
    ADD_SLOT( MPI_Comm        , mpi_dup       , INPUT , REQUIRED);
    ADD_SLOT( HistogramWorker , insitu_worker , OUTPUT );

    inline void execute () override final
    {
      // create and start worker thread
      // TODO 
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "histoseq_thread_init" , make_simple_operator< TutorialHistoThreadInit > );
  }

}

