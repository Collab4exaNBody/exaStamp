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

  template<
    class GridT,
    class HistField,
    class = std::enable_if_t< GridHasField<GridT,HistField>::value >
    >
  struct TutorialHistoThreadRun : public OperatorNode
  {      
    ADD_SLOT( GridT     , grid      , INPUT , REQUIRED );
    ADD_SLOT( long      , samples   , INPUT , 1000 );
    ADD_SLOT( HistogramWorker, insitu_worker , INPUT_OUTPUT );

    inline void execute () override final
    {
      static constexpr onika::soatl::FieldId<HistField> hist_field{};

      // prepare (copy) input data for worker thread
      {
        std::lock_guard<std::mutex> lk(insitu_worker->m);
        lout << "MAIN: copy input data" << std::endl;
                
        size_t n_cells = grid->number_of_cells();
        auto cells = grid->cells();
        IJK dims = grid->dimension();
        int gl = grid->ghost_layers();
        
        insitu_worker->input.clear();
        for(size_t i=0;i<n_cells;i++)
        {
          IJK loc = grid_index_to_ijk(dims,i);
          if( loc.i>=gl && loc.i<(dims.i-gl) && loc.j>=gl && loc.j<(dims.j-gl) && loc.k>=gl && loc.k<(dims.k-gl) )
          {
            size_t n = cells[i].size();          
            const auto * __restrict__ value_ptr = cells[i][hist_field];
            insitu_worker->input.insert( insitu_worker->input.end() , value_ptr , value_ptr+n );
          }
        }
        insitu_worker->output.m_data.resize( *samples );
        insitu_worker->ready = true;
      }
      
      // notify worker thread data is ready
      lout << "MAIN: notify input data is ready" << std::endl;
      insitu_worker->cv.notify_one();      
    }

  };

  template<typename GridT> using TutorialHistoThreadRunEnergy = TutorialHistoThreadRun<GridT,field::_ep>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(histogram_seq_thread_run)
  {
   OperatorNodeFactory::instance()->register_factory( "histoseq_thread_run" , make_grid_variant_operator< TutorialHistoThreadRunEnergy > );
  }

}

