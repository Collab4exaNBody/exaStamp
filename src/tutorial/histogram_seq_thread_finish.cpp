/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/histogram.h>
#include <onika/mpi/data_types.h>

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
  ONIKA_AUTORUN_INIT(histogram_seq_thread_finish)
  {
   OperatorNodeFactory::instance()->register_factory( "histoseq_thread_finish" , make_simple_operator< TutorialHistoThreadFinish > );
  }

}

