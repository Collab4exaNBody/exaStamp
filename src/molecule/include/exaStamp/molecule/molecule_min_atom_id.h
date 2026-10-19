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

#pragma once

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/mpi/grid_update_ghosts.h>
#include <exanb/mpi/update_ghost_functors.h>
#include <exaStamp/molecule/id_map.h>
#include <exaStamp/fields.h>
#include <onika/log.h>
#include <mpi.h>
#include <limits>
#include <algorithm>

namespace exaStamp
{
  using namespace exanb;

  // Sets the idmol field of every particle to the smallest atom id of its molecule (connected component of the bond
  // graph given by cmol), across MPI processes. idmol is used as work storage : its previous content is lost.
  // pecfunc / peqfunc are the calling operator's parallel execution context and queue accessors.
  template<class LDBG, class GridT, class PecFuncT, class PeqFuncT>
  inline void propagate_molecule_min_atom_id( LDBG& ldbg, MPI_Comm comm, GhostCommunicationScheme& ghost_comm_scheme,
                                              GridT* gridp, const Domain& domain, const UpdateGhostConfig& update_ghost_config,
                                              const PecFuncT& pecfunc, const PeqFuncT& peqfunc,
                                              const IdMap& id_map, const IdMapGhosts& id_map_ghosts )
  {
    GridT& grid_ref = *gridp;
    auto * grid = &grid_ref;
    auto cells = grid->cells_accessor();
    size_t n_cells = grid->number_of_cells();
    const auto idmol_field = grid->field_accessor( field::idmol );
    const auto cmol_field = grid->field_accessor( field::cmol );
    static constexpr uint64_t null_loc = std::numeric_limits<uint64_t>::max();
    static constexpr size_t null_index = std::numeric_limits<size_t>::max();
    static constexpr GridCellValues * null_cell_values = nullptr;

      // preamble : assign particle's molid particles' ids, so that we start with individual molecules
      // of one atom each, and we will then merge those bound by intramolecular bond
      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++) cells[cell_i][idmol_field][p_i] = cells[cell_i][field::id][p_i];
      }
      ldbg << "molecule ids: id_map size = "<<id_map.size()<<" , id_map_ghosts size = "<<id_map_ghosts.size() <<std::endl;

      // first multi pass algorithm : propagate molecule ids through bonds, keeping only the minimum id a the molecule id
      long update_count = 0;
      int pass = 0;
      do
      {
        const auto ghost_update_fields = onika::make_flat_tuple( idmol_field ); 
        auto ghost_comm_buffers = UpdateGhostsUtils::make_forward_ghost_update_manager( cells, ghost_update_fields );
        grid_update_ghosts( exanb::ldbg, comm, ghost_comm_scheme, gridp, domain, null_cell_values,
                            ghost_comm_buffers, pecfunc, peqfunc, ghost_update_fields, update_ghost_config );

        update_count = 0;
        for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
        {
#         ifndef NDEBUG
          const auto * __restrict__ ids = cells[cell_i][field::id];
#         endif
          const auto *  __restrict__ cmol = cells[cell_i][cmol_field];
          auto * __restrict__ mol_ids   = cells[cell_i][idmol_field];
          size_t n = cells[cell_i].size();
          for(size_t p_i=0;p_i<n;p_i++)
          {
            assert( mol_ids[p_i] <= ids[p_i] );
            uint64_t molidmin = mol_ids[p_i];
            for(size_t j=0;j<cmol[p_i].size();j++)
            {
              const auto conloc = atom_from_idmap_if_found( cmol[p_i][j] , id_map , id_map_ghosts , null_loc );
              if( conloc != null_loc )
              {
                size_t cell=null_index, pos=null_index;
                decode_cell_particle( conloc , cell, pos );
                assert( cell!=null_index && pos!=null_index );
                const auto con_mol_id = cells[cell][idmol_field][pos];
                molidmin = std::min( molidmin , con_mol_id );
                if( molidmin < con_mol_id )
                {
                  cells[cell][idmol_field][pos] = molidmin;
                  ++ update_count;
                }
              }
            }
            if( molidmin < mol_ids[p_i] )
            {
              mol_ids[p_i] = molidmin;
              ++ update_count;
              // ldbg << "cell #"<<cell_i<<" particle #"<<p_i<<" molid "<<mol_ids[p_i]<<" -> "<<molidmin<<std::endl;
            }
          }
        }
        MPI_Allreduce( MPI_IN_PLACE, &update_count , 1 , MPI_LONG , MPI_SUM , comm );
        ++ pass;
        ldbg << "molecule ids: connect molecules, pass="<<pass<<" , updates="<<update_count<<std::endl;
      } while( update_count > 0 );
  }

}
