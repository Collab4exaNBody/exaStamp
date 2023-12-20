#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/histogram.h>
#include <exanb/mpi/data_types.h>

#include <mpi.h>

// The following code will use the struct/class and donctions defined in this include file.
#include "histogram_worker.h"

namespace exaStamp
{

  template<
    class GridT,
    class HistField,
    class = std::enable_if_t< GridHasField<GridT,HistField>::value >
    >
  struct TutorialHistoNTOperator : public OperatorNode
  {
    // System allowing to use ExaStamp data and connect it to other functions
    ADD_SLOT( MPI_Comm  , mpi_dup   , INPUT , REQUIRED);
    ADD_SLOT( GridT     , grid      , INPUT , REQUIRED );
    ADD_SLOT( long      , samples   , INPUT , 1000 );
    ADD_SLOT( long      , timestep  , INPUT , 1 );
    ADD_SLOT( Histogram<> , histogram , OUTPUT );
    ADD_SLOT( HistogramWorker, insitu_worker , OUTPUT ); // insitu_worker is an object of class HistogramWorker

    inline void execute () override final
    {
      static constexpr onika::soatl::FieldId<HistField> hist_field{};

      // create and start worker thread
      if( insitu_worker->m_thread == nullptr ) // Creation of the worker thread during the first iteration (it should not be created each time)
      {	
        // TODO : initialization of the insitu_worker->m_thread (using HistogramWorker::run, ... and the dupliacte MPI communicator mpi_dup)
	// insitu_worker->m_thread = 
      }      

      // Prepare (copy) input data for worker thread
      {
	// TODO : A creation of a lock_guard  mutex "lk" is required here using the insitu_worker mutex.
        // ...
        lout << "MAIN: copy input data" << std::endl;
                
        size_t n_cells = grid->number_of_cells();
        auto cells = grid->cells();
        IJK dims = grid->dimension();
        int gl = grid->ghost_layers();
        
        insitu_worker->input.clear(); // Clearing the input data after a previous iteration
        for(size_t i=0;i<n_cells;i++)
        {
          IJK loc = grid_index_to_ijk(dims,i);
          if( loc.i>=gl && loc.i<(dims.i-gl) && loc.j>=gl && loc.j<(dims.j-gl) && loc.k>=gl && loc.k<(dims.k-gl) )
          {
            // TODO: uncomment the following 2 lines
            size_t n = cells[i].size();          
            const auto * __restrict__ value_ptr = cells[i][hist_field];
            
	    std::ignore = value_ptr[n-1]; // avoid unnecessary warnings
	    // TODO : Fill the vector input de insitu_worker
            // insitu_worker->input.insert(...
          }
        }
        insitu_worker->output.m_data.resize( *samples );

	// TODO : ready ?
        // ...
      }
      
      // notify worker thread data is ready
      lout << "MAIN: notify input data is ready" << std::endl;

      // TODO : action with the condition variable
      // ...
      
      // wait for the worker
      {
          std::unique_lock<std::mutex> lk(insitu_worker->m);
          insitu_worker->cv.wait(lk, [this]{return insitu_worker->processed;});
          insitu_worker->processed = false;
      }
      std::cout << "MAIN: Worker done" << std::endl;

      // TODO : fill  histogram->m_min_val,  histogram->m_max_val and histogram->m_data
      // histogram->m_min_val = ...
      // histogram->m_max_val = ...
      // histogram->m_data = ...
    }

  };

  template<typename GridT> using HistoEnergyNT = TutorialHistoNTOperator<GridT,field::_ep>;
  template<typename GridT> using HistoVxNT = TutorialHistoNTOperator<GridT,field::_vx>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "histoseqnt_energy" , make_grid_variant_operator< HistoEnergyNT > );
   OperatorNodeFactory::instance()->register_factory( "histoseqnt_vx" , make_grid_variant_operator< HistoVxNT > );
  }

}

