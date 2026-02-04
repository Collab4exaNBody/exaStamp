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
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/position_long_term_backup.h>

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT,  class = AssertGridHasFields< GridT, field::_id> >
  class PositionRestoreLongTerm : public OperatorNode
  {
    ADD_SLOT( GridT                  , grid        , INPUT_OUTPUT );
    ADD_SLOT( Domain                 , domain      , INPUT , REQUIRED );  // 
    ADD_SLOT( PositionLongTermBackup , backup_r_lt , INPUT , REQUIRED );

  public:
    inline void execute ()  override final
    {    
      IJK dims = grid->dimension();
      auto cells = grid->cells();
      ssize_t gl = grid->ghost_layers();

      ldbg << "xform backup : " << backup_r_lt->m_xform << std::endl;
      ldbg << "xform        : " << domain->xform() << std::endl;      

      Mat3d mat;
      if( domain->xform() == backup_r_lt->m_xform )
      {
        mat = make_identity_matrix();
      }
      else
      {
        mat = domain->inv_xform() * backup_r_lt->m_xform;
      }

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
	        const size_t n_particles = cells[i].size();
	        size_t start = backup_r_lt->m_cell_offset[i];
#         ifndef NDEBUG
          const auto* __restrict__ ids = cells[i][field::id];
#         endif
          auto* __restrict__ rx = cells[i][field::rx];
          auto* __restrict__ ry = cells[i][field::ry];
          auto* __restrict__ rz = cells[i][field::rz];
#         pragma omp simd
          for(size_t j=0;j<n_particles;j++)
          {
            assert( backup_r_lt->m_ids[start+j] == ids[j] );
	    //            Vec3d r = inverse(backup_r_lt->m_xform) * backup_r_lt->m_positions[start+j];
            Vec3d r = backup_r_lt->m_positions[start+j];	    
            rx[j] = r.x;
            ry[j] = r.y;
            rz[j] = r.z;
          }
        }
        GRID_OMP_FOR_END
      }
    }

  };

  template<class GridT> using PositionRestoreLongTermTmpl = PositionRestoreLongTerm<GridT>;

 // === register factories ===  
  ONIKA_AUTORUN_INIT(restore_r_lt)
  {
   OperatorNodeFactory::instance()->register_factory( "restore_r_lt", make_grid_variant_operator< PositionRestoreLongTermTmpl > );
  }

}

