#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/histogram.h>

#include <mpi.h>
#include <exanb/mpi/data_types.h>

#include <memory>
#include <vector>
#include <type_traits>
#include <limits>

namespace exaStamp
{

  template<
    class GridT,
    class HistField,
    class = std::enable_if_t< GridHasField<GridT,HistField>::value >
    >
  struct HistogramOperator : public OperatorNode
  {
    using ValueType = typename onika::soatl::FieldId<HistField>::value_type ;
      
    ADD_SLOT( MPI_Comm  , mpi       , INPUT , REQUIRED );
    ADD_SLOT( GridT     , grid      , INPUT , REQUIRED );
    ADD_SLOT( long      , samples   , INPUT , 1000 );
    ADD_SLOT( bool      , ghost     , INPUT , false );
    ADD_SLOT( Histogram<> , histogram , OUTPUT );

    inline void execute () override final
    {
      static constexpr onika::soatl::FieldId<HistField> hist_field{};

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
      ValueType min_val = std::numeric_limits<ValueType>::max();
      ValueType max_val = std::numeric_limits<ValueType>::lowest();

#     pragma omp parallel
      {
        ValueType local_min_val = std::numeric_limits<ValueType>::max();
        ValueType local_max_val = std::numeric_limits<ValueType>::lowest();
        
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc /*, reduction(min:min_val) reduction(max:max_val) */ )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();          
          const ValueType * __restrict__ value_ptr = cells[i][hist_field];

//#         pragma omp simd reduction(min:min_val) reduction(max:max_val)
          for(size_t j=0;j<n;j++)
          {
            ValueType x = value_ptr[j];
            local_min_val = std::min( local_min_val , x );
            local_max_val = std::max( local_max_val , x );
          }
        }
        GRID_OMP_FOR_END
#       pragma omp critical
        {
          min_val = std::min( local_min_val , min_val );
          max_val = std::max( local_max_val , max_val );
        }
      }

      // MPI min/max
      if( nprocs > 1 )
      {
        ValueType tmp[2] = { -min_val , max_val };
        MPI_Allreduce(MPI_IN_PLACE,tmp,2, exanb::mpi_datatype<ValueType>() ,MPI_MAX,comm);
        min_val = - tmp[0];
        max_val = tmp[1];
      }

      // std::cout << "min="<<min_val<<", max="<<max_val<<std::endl;

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

          const ValueType * __restrict__ value_ptr = cells[i][hist_field];
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            ssize_t bin = static_cast<size_t>( ( (value_ptr[j]-min_val) * hist_size ) / ( max_val - min_val ) );
            if( bin < 0 ) { bin=0; }
            if( bin >= static_cast<ssize_t>(hist_size) ) { bin = hist_size-1; }
            local_hist[bin] += 1.0;
          }
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

  template<typename GridT> using HistogramEnergy = HistogramOperator<GridT,field::_ep>;
  template<typename GridT> using HistogramCharge = HistogramOperator<GridT,field::_charge>;
  template<typename GridT> using HistogramVx = HistogramOperator<GridT,field::_vx>;
  template<typename GridT> using HistogramRx = HistogramOperator<GridT,field::_rx>;
  template<typename GridT> using HistogramRy = HistogramOperator<GridT,field::_ry>;
  template<typename GridT> using HistogramRz = HistogramOperator<GridT,field::_rz>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "histogram_energy" , make_grid_variant_operator< HistogramEnergy > );
   OperatorNodeFactory::instance()->register_factory( "histogram_charge" , make_grid_variant_operator< HistogramCharge > );
   OperatorNodeFactory::instance()->register_factory( "histogram_vx" , make_grid_variant_operator< HistogramVx > );
   OperatorNodeFactory::instance()->register_factory( "histogram_rx" , make_grid_variant_operator< HistogramRx > );
   OperatorNodeFactory::instance()->register_factory( "histogram_ry" , make_grid_variant_operator< HistogramRy > );
   OperatorNodeFactory::instance()->register_factory( "histogram_rz" , make_grid_variant_operator< HistogramRz > );
  }

}

