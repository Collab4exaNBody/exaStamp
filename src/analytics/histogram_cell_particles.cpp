#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/histogram.h>

#include <mpi.h>

#include <memory>
#include <vector>
#include <type_traits>
#include <limits>

namespace exaStamp
{

  template<class GridT>
  class HistogramCellParticles : public OperatorNode
  {
    ADD_SLOT( MPI_Comm  , mpi       , INPUT , REQUIRED );
    ADD_SLOT( GridT     , grid      , INPUT , REQUIRED );
    ADD_SLOT( long      , samples   , INPUT , 25 );
    ADD_SLOT( bool      , ghost     , INPUT , false );
    ADD_SLOT( Histogram<> , histogram , OUTPUT );

  public:
    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      int nprocs = 1;
      int rank = 0;
      MPI_Comm_size(comm,&nprocs);
      MPI_Comm_rank(comm,&rank);

      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();
      if( *ghost ) { gl = 0; }

      // min max computation
      long min_val = std::numeric_limits<long>::max();
      long max_val = std::numeric_limits<long>::lowest();

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, reduction(min:min_val) reduction(max:max_val) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          long n = cells[i].size();
          min_val = std::min( min_val , n );
          max_val = std::max( max_val , n );
        }
        GRID_OMP_FOR_END
      }

      // MPI min/max
      if( nprocs > 1 )
      {
        long tmp[2] = { -min_val , max_val };
        MPI_Allreduce(MPI_IN_PLACE,tmp,2, MPI_LONG ,MPI_MAX,comm);
        min_val = - tmp[0];
        max_val = tmp[1];
      }

      // hitogram counting
      size_t hist_size = *samples;
      histogram->m_min_val = min_val;
      histogram->m_max_val = max_val;
      histogram->m_data.resize( hist_size );

      int max_nt = omp_get_max_threads();
      double* per_thread_histogram[max_nt];

#     pragma omp parallel 
      {
        int nt = omp_get_num_threads();
        int tid = omp_get_thread_num();
        assert( tid<max_nt && nt<=max_nt );
        
        // stack allocated, per thread histogram
        double local_hist[hist_size];
        for(size_t i=0;i<hist_size;i++) { local_hist[i] = 0.0; }

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc)
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          ssize_t bin = ( (n-min_val) * hist_size ) / ( max_val - min_val );
          if( bin < 0 ) { bin=0; }
          if( bin >= static_cast<ssize_t>(hist_size) ) { bin = hist_size-1; }
          local_hist[bin] += 1.0;
        }
        GRID_OMP_FOR_END

        per_thread_histogram[tid] = local_hist;
        size_t start = ( hist_size * tid ) / nt;
        size_t end = ( hist_size * (tid+1) ) / nt;

#       pragma omp barrier

        double* h = per_thread_histogram[0];        
        for(size_t i=start;i<end;i++)
        {
          histogram->m_data[i] = h[i];
        }
        for(int t=1;t<nt;t++)
        {
          h = per_thread_histogram[t]; 
          for(size_t i=start;i<end;i++)
          {
            histogram->m_data[i] += h[i];
          }
        }

      }

      if( nprocs > 1 )
      {
        MPI_Allreduce(MPI_IN_PLACE,histogram->m_data.data(),hist_size,MPI_DOUBLE,MPI_SUM,comm);
      }

    }
  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(histogram_cell_particles)
  {
   OperatorNodeFactory::instance()->register_factory( "histogram_cell_particles" , make_grid_variant_operator< HistogramCellParticles > );
  }

}

