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

#include <memory>
#include <vector>
#include <type_traits>
#include <limits>

#include <mpi.h>

namespace exaStamp
{

  template<
    class GridT,
    class HistField,
    class = std::enable_if_t< GridHasField<GridT,HistField>::value >
    >
  struct TutorialHistoOperator : public OperatorNode
  {
    using ValueType = typename onika::soatl::FieldId<HistField>::value_type ;
      
    ADD_SLOT( MPI_Comm  , mpi       , INPUT , REQUIRED);
    ADD_SLOT( GridT     , grid      , INPUT , REQUIRED );
    ADD_SLOT( long      , samples   , INPUT , 1000 );
    ADD_SLOT( long      , timestep  , INPUT , 1 );
    ADD_SLOT( Histogram<> , histogram , OUTPUT );

    inline void execute () override final
    {
      static constexpr onika::soatl::FieldId<HistField> hist_field{};
      MPI_Comm comm = *mpi;
      size_t nsamples = *samples;
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      int gl = grid->ghost_layers();
      size_t n_cells = grid->number_of_cells();
            
      int nprocs = 1;
      int rank = 0;
      MPI_Comm_size(comm,&nprocs);
      MPI_Comm_rank(comm,&rank);

//      lout << "timestep = " << *timestep << std::endl;
      // ============ modification area ============================
      
      // min max computation
      ValueType min_val = std::numeric_limits<ValueType>::max();
      ValueType max_val = std::numeric_limits<ValueType>::lowest();
//      lout << "A. min="<<min_val<<", max="<<max_val<<std::endl;

      for(size_t i=0;i<n_cells;i++)
      {
        IJK loc = grid_index_to_ijk(dims,i);
        if( loc.i>=gl && loc.i<(dims.i-gl) && loc.j>=gl && loc.j<(dims.j-gl) && loc.k>=gl && loc.k<(dims.k-gl) )
        {
          size_t n = cells[i].size();          
          const ValueType * __restrict__ value_ptr = cells[i][hist_field];
          for(size_t j=0;j<n;j++)
          {
            ValueType x = value_ptr[j];
            min_val = std::min( min_val , x );
            max_val = std::max( max_val , x );
          }
        }
      }

//      lout << "B. min="<<min_val<<", max="<<max_val<<std::endl;

      // MPI min/max
      if( nprocs > 1 )
      {
        ValueType tmp[2] = { -min_val , max_val };
        MPI_Allreduce(MPI_IN_PLACE,tmp,2, onika::mpi::mpi_datatype<ValueType>() ,MPI_MAX,comm);
        min_val = - tmp[0];
        max_val = tmp[1];
      }
      //lout << "C. min="<<min_val<<", max="<<max_val<<std::endl;

      // hitogram counting
      size_t hist_size = nsamples;
      histogram->m_min_val = min_val;
      histogram->m_max_val = max_val;
      histogram->m_data.resize( hist_size );

      for(size_t i=0;i<n_cells;i++)
      {
        IJK loc = grid_index_to_ijk(dims,i);
        if( loc.i>=gl && loc.i<(dims.i-gl) && loc.j>=gl && loc.j<(dims.j-gl) && loc.k>=gl && loc.k<(dims.k-gl) )
        {
          size_t n = cells[i].size();          
          const ValueType * __restrict__ value_ptr = cells[i][hist_field];
          for(size_t j=0;j<n;j++)
          {
            ValueType x = value_ptr[j];
            for(int k=0;k<10;k++)
            {
              x = std::exp( std::sin( std::sqrt(x) ) );
            }
            if(x!=666.0) x = value_ptr[j];
            ssize_t bin = static_cast<ssize_t>( ( (x-min_val) * hist_size ) / ( max_val - min_val ) );
            if( bin < 0 ) { bin=0; }
            if( bin >= static_cast<ssize_t>(hist_size) ) { bin = hist_size-1; }
            histogram->m_data[bin] += 1.0;
          }
        }
      }

      if( nprocs > 1 )
      {
        MPI_Allreduce(MPI_IN_PLACE,histogram->m_data.data(),hist_size,MPI_DOUBLE,MPI_SUM,comm);
      }

      // ============ end of modification area ============================
    }
  };

  template<typename GridT> using HistoEnergy = TutorialHistoOperator<GridT,field::_ep>;
  template<typename GridT> using HistoVx = TutorialHistoOperator<GridT,field::_vx>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(histogram_seq)
  {
   OperatorNodeFactory::instance()->register_factory( "histoseq_energy" , make_grid_variant_operator< HistoEnergy > );
   OperatorNodeFactory::instance()->register_factory( "histoseq_vx" , make_grid_variant_operator< HistoVx > );
  }

}

