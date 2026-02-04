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
#include <exanb/core/grid.h>
#include <onika/math/basic_types.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/position_long_term_backup.h>

#include <onika/soatl/packed_field_arrays.h>
#include <onika/mpi/xs_data_move.h>

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT, class = AssertGridHasFields< GridT, field::_id >, class = AssertGridHasFields< GridT, field::_rxf, field::_ryf, field::_rzf> >
  struct FilteredPositionBackupLongTermMoveData : public OperatorNode
  {
    ADD_SLOT( MPI_Comm               , mpi         , INPUT , REQUIRED );
    ADD_SLOT( GridT                  , grid        , INPUT , REQUIRED );
    ADD_SLOT( PositionLongTermBackup , backup_r_rf_lt , INPUT_OUTPUT);

    inline void execute ()  override final
    {
      GridT& grid = *(this->grid);
      size_t n_cells = grid.number_of_cells();
      IJK dims = grid.dimension();
      auto cells = grid.cells();
      ssize_t gl = grid.ghost_layers();
      
      size_t backup_n_ids = backup_r_rf_lt->m_ids.size();

      uint64_t idMin = backup_r_rf_lt->m_idMin;
      uint64_t idMax = backup_r_rf_lt->m_idMax;

      backup_r_rf_lt->m_cell_offset.assign( n_cells+1 , 0 );
      size_t total_particles = 0;
      GRID_FOR_BEGIN(dims,i,loc)
      {
        size_t n_particles = 0;
        if( ! grid.is_ghost_cell(loc) ) { n_particles = cells[i].size(); }
        backup_r_rf_lt->m_cell_offset[i] = total_particles;
        total_particles += n_particles;
      }
      GRID_FOR_END
      backup_r_rf_lt->m_cell_offset[n_cells] = total_particles;

      std::vector< uint64_t > current_ids( total_particles , std::numeric_limits<uint64_t>::max() );
#     pragma omp parallel
      {      
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
	        const size_t n_particles = cells[i].size();
	        size_t start = backup_r_rf_lt->m_cell_offset[i];
          const auto* __restrict__ id = cells[i][field::id];
#         pragma omp simd 
          for(size_t j=0;j<n_particles;j++)
          {
            current_ids[start+j] = id[j];
          }
        }
        GRID_OMP_FOR_END
      }

#     ifndef NDEBUG
      for(auto x : current_ids)
      {
        assert( x != std::numeric_limits<uint64_t>::max() );
        assert( x>=idMin && x<idMax );
      }
      for(auto x : backup_r_rf_lt->m_ids)
      {
        assert( x != std::numeric_limits<uint64_t>::max() );
        assert( x>=idMin && x<idMax );
      }
#     endif

      assert( backup_n_ids == backup_r_rf_lt->m_ids.size() );

      std::vector<int> send_indices;     // resized to localElementCountBefore, contain indices in 'before' array to pack into send buffer
      std::vector<int> send_count;       // resized to number of processors in comm, unit data element count (not byte size)
      std::vector<int> send_displ;       // resized to number of processors in comm, unit data element count (not byte size)
      std::vector<int> recv_indices;     // resized to localElementCountAfter
      std::vector<int> recv_count;       // resized to number of processors in comm, unit data element count (not byte size)
      std::vector<int> recv_displ;       // resized to number of processors in comm, unit data element count (not byte size)
      onika::mpi::communication_scheme_from_ids(*mpi, idMin, idMax /*exclusive max(ids)+1*/, backup_n_ids, backup_r_rf_lt->m_ids.data(), total_particles, current_ids.data(),
                                                send_indices, send_count, send_displ, recv_indices, recv_count, recv_displ );

#     ifndef NDEBUG
      backup_r_rf_lt->m_ids.resize( std::max( backup_r_rf_lt->m_ids.size(), current_ids.size() ) );
      onika::mpi::data_move(*mpi, send_indices, send_count, send_displ, recv_indices, recv_count, recv_displ, backup_r_rf_lt->m_ids.data(), backup_r_rf_lt->m_ids.data() );
      backup_r_rf_lt->m_ids.resize(total_particles);
      for(size_t i=0;i<total_particles;i++)
      {
        assert( backup_r_rf_lt->m_ids[i] == current_ids[i] );
      }
      current_ids.clear();
#     else
      backup_r_rf_lt->m_ids = std::move( current_ids );      
#     endif

      // resize position buffer to move data in place
      backup_r_rf_lt->m_positions.resize( std::max(backup_n_ids,total_particles) );
      backup_r_rf_lt->m_filtered_positions.resize( std::max(backup_n_ids,total_particles) );      

      // move position data across processors
      onika::mpi::data_move(*mpi, send_indices, send_count, send_displ, recv_indices, recv_count, recv_displ, backup_r_rf_lt->m_positions.data(), backup_r_rf_lt->m_positions.data() );
      onika::mpi::data_move(*mpi, send_indices, send_count, send_displ, recv_indices, recv_count, recv_displ, backup_r_rf_lt->m_filtered_positions.data(), backup_r_rf_lt->m_filtered_positions.data() );      
      backup_r_rf_lt->m_positions.resize( total_particles );
      backup_r_rf_lt->m_filtered_positions.resize( total_particles );      

      assert( backup_r_rf_lt->m_positions.size() == backup_r_rf_lt->m_ids.size() );
      assert( backup_r_rf_lt->m_filtered_positions.size() == backup_r_rf_lt->m_ids.size() );      

      backup_r_rf_lt->m_positions.shrink_to_fit();
      backup_r_rf_lt->m_filtered_positions.shrink_to_fit();      
      backup_r_rf_lt->m_ids.shrink_to_fit();      

#     ifndef NDEBUG
#     pragma omp parallel
      {      
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
	        const size_t n_particles = cells[i].size();
	        size_t start = backup_r_rf_lt->m_cell_offset[i];
          const auto* __restrict__ ids = cells[i][field::id];
          for(size_t j=0;j<n_particles;j++)
          {
            assert( ids[j] == backup_r_rf_lt->m_ids[ start + j ] );
            /*
            assert( backup_r_rf_lt->m_positions[start+j].x == cells[i][field::rx0][j] );
            assert( backup_r_rf_lt->m_positions[start+j].y == cells[i][field::ry0][j] );
            assert( backup_r_rf_lt->m_positions[start+j].z == cells[i][field::rz0][j] );
            */
          }
        }
        GRID_OMP_FOR_END
      }
#     endif

    }

  };

  template<class GridT> using FilteredPositionBackupLongTermMoveDataTmpl = FilteredPositionBackupLongTermMoveData<GridT>;  

 // === register factories ===  
  ONIKA_AUTORUN_INIT(backup_r_rf_lt_move_data)
  {
   OperatorNodeFactory::instance()->register_factory( "backup_r_rf_lt_move_data", make_grid_variant_operator< FilteredPositionBackupLongTermMoveDataTmpl > );
  }

}

