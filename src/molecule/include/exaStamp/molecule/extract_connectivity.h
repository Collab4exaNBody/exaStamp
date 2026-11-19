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
#include <onika/log.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/grid_fields.h>

#include <exaStamp/molecule/mol_connectivity.h>
#include <exaStamp/molecule/id_map.h>

#include <exanb/core/parallel_grid_algorithm.h>

namespace exaStamp
{

  namespace extract_connectivity_details
  {
    struct PerThreadCount
    {
      size_t bond_count=0;
      size_t angle_count=0;
      size_t torsion_count=0;
      size_t improper_count=0;
    };
    
    template<typename GridT>
    static inline bool is_leader_cell(bool cell_is_inner, const GridT& grid, const IJK& cell_loc, const IJK& domain_cells)
    {
      if( cell_is_inner )
      {
        return true;
      }
      else
      {
        bool is_leader = true;
        for(int di=-1;di<=1;di++) for(int dj=-1;dj<=1;dj++) for(int dk=-1;dk<=1;dk++) if( di!=0 || dj!=0 || dk!=0 )
        {
          IJK loc = cell_loc + IJK{di*domain_cells.i,dj*domain_cells.j,dk*domain_cells.k};
          if( grid.contains(loc) )
          {
            bool loc_is_here = !grid.is_ghost_cell(loc);
            bool lex_order = lexicographic_order(loc,cell_loc);
            if( loc_is_here || lex_order )
            {
              is_leader = false;
            }
          }
        }
        return is_leader;
      }
    }
    
  }
  // Sets the connectivity of the molecules and adds an id to the molecules
  template<typename GridT> inline  void extract_connectivity(
    IJK domain_cells,
    GridT& grid,
    ChemicalBonds& result_bonds,
    ChemicalAngles& result_angles,
    ChemicalTorsions& result_torsions,
    ChemicalImpropers& result_impropers,
    IdMap& id_map,
    IdMapGhosts& id_map_ghosts)
  {
    using namespace extract_connectivity_details;
    static constexpr uint64_t atom_not_found = std::numeric_limits<uint64_t>::max();

    assert( ! omp_in_parallel() );

    auto cells = grid.cells_accessor();

    IJK dims = grid.dimension();
    
    int max_nt = omp_get_max_threads();
    PerThreadCount per_thread_counts [max_nt];
    ChemicalBonds per_thread_bonds [max_nt];
    ChemicalAngles per_thread_angles [max_nt];
    ChemicalTorsions per_thread_torsions [max_nt];
    ChemicalImpropers per_thread_impropers [max_nt];

    // reuse previously allocated storage if any
    per_thread_bonds[0] = std::move( result_bonds );
    per_thread_bonds[0].clear();
    
    per_thread_angles[0] = std::move( result_angles );
    per_thread_angles[0].clear();
    
    per_thread_torsions[0] = std::move( result_torsions );
    per_thread_torsions[0].clear();

    per_thread_impropers[0] = std::move( result_impropers );
    per_thread_impropers[0].clear();

    std::set<int64_t> atomIds;
    std::set<int64_t> atomGhostIds;
    const size_t n_cells = grid.number_of_cells();
    for(size_t cell_i=0;cell_i<n_cells;cell_i++)
    {
      size_t n_particles = cells[cell_i].size();
      const auto * __restrict__ ids = cells[cell_i][field::id];
      if( grid.is_ghost_cell(cell_i) ) for(size_t p=0; p<n_particles; ++p) atomGhostIds.insert( ids[p] );
      else for(size_t p=0; p<n_particles; ++p) atomIds.insert( ids[p] );
    }

    auto field_cmol = grid.field_accessor( field::cmol );

#   pragma omp parallel
    {
      int nt = omp_get_num_threads();
      int tid = omp_get_thread_num();
      if(tid>=max_nt || nt>max_nt)
      {
        fatal_error() << "Internal error: bad number of threads" << std::endl;
      }
      auto& bonds = per_thread_bonds[tid];
      auto& angles = per_thread_angles[tid];
      auto& torsions = per_thread_torsions[tid];
      auto& impropers = per_thread_impropers[tid];

      GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(dynamic) nowait )
      {
        bool cell_is_here = ! grid.is_ghost_cell(cell_loc); // 'here' meaning in the central area (not ghost) of the local processor
        size_t n_particles = cells[cell_i].size();
        const auto * __restrict__ ids = cells[cell_i][field::id];
        const auto * __restrict__ cmols = cells[cell_i].field_pointer_or_null(field_cmol);
        const auto * __restrict__ types = cells[cell_i][field::type];
        
        for(size_t p=0; p<n_particles; ++p)
        {
          uint64_t p_id = ids[p];
          const auto& neigh = cmols[p];
          
          if( cell_is_here )
          {
            // first check connectivity consistency (all connex atoms must be found, at least in the ghost area)
            for(int nbh=0;nbh<4;nbh++)
            {
              auto nbh_id = neigh[nbh];
              if(nbh_id != atom_not_found )
              {
                if( atomIds.find(nbh_id)==atomIds.end() && atomGhostIds.find(nbh_id)==atomGhostIds.end() )
                {
#                 pragma omp critical(dbg_mesg)
                  {
                    fatal_error() << "Bad atom connectivity: cell @"<<cell_loc<<" atom #"<<p<<" (id="<<ids[p]<<") : cmol["<<nbh<<"]="<<nbh_id<<" not found in grid"<<std::endl;
                  }
                }
              }
            }
          }
          
          const bool p_is_leader = ( encode_cell_particle(cell_i,p,types[p]) == atom_from_idmap(p_id,id_map,id_map_ghosts) );
          if( p_is_leader )
          {        
            // IMPROPERS CASE ---------------------------------------------------------------
            // There are exactly 3 neighbors. In the (common) case where an atom has four neighbors
            // we do not consider (for now) that an improper angle potential can be applied.
            if(neigh[2]!=std::numeric_limits<uint64_t>::max() && neigh[3]==std::numeric_limits<uint64_t>::max())
            {
              // if ghost and a no ghost version exist, ignore this
              // otherwise, detect if we are the preffered ghost cell
              bool n0_is_here = ( id_map.find(neigh[0]) != id_map.end() );
              bool n1_is_here = ( id_map.find(neigh[1]) != id_map.end() );
              bool n2_is_here = ( id_map.find(neigh[2]) != id_map.end() );
              if( cell_is_here || n0_is_here || n1_is_here || n2_is_here )
              {
                impropers.push_back( ChemicalImproper{p_id, neigh[0], neigh[1], neigh[2]} );
              }
            }

            // first neighbor step
            for(auto id_a:neigh)
            {
              assert( id_a != p_id );
              if(id_a==std::numeric_limits<uint64_t>::max()) break;
              bool id_a_is_here = ( id_map.find(id_a) != id_map.end() );
              if( cell_is_here )
              {
                if( p_id < id_a )         bonds.push_back( ChemicalBond{p_id,id_a} );
                else if( ! id_a_is_here ) bonds.push_back( ChemicalBond{id_a,p_id} );
              }
              
              // second neighbor step
              uint64_t enc_atom_a = atom_from_idmap_if_found(id_a,id_map,id_map_ghosts,atom_not_found);
              // can be outside of processor's box if and only if none of preceding atoms in the chain are in the inner box (all are ghosts)
              assert( enc_atom_a!=atom_not_found || !cell_is_here );
              if( enc_atom_a!=atom_not_found )
              {
                size_t cell_a=-1, pos_a=-1;
                decode_cell_particle( enc_atom_a , cell_a, pos_a);
                const auto& neigh_a = cells[cell_a][field_cmol][pos_a];
                for(auto id_b:neigh_a) if(id_b!=p_id)
                {
                  assert( id_b != id_a );
                  if(id_b==std::numeric_limits<uint64_t>::max()) break;                
                  bool id_b_is_here = ( id_map.find(id_b) != id_map.end() );                
                  bool angle_must_be_added = (cell_is_here || id_a_is_here || id_b_is_here) && (p_id!=id_b);
                  bool angle_owner = (cell_is_here && !id_b_is_here) || (cell_is_here==id_b_is_here && p_id<id_b);
                  if( angle_must_be_added && angle_owner )
                  {
                    if( p_id < id_b ) angles.push_back( ChemicalAngle{ p_id , id_a , id_b } );
                    else              angles.push_back( ChemicalAngle{ id_b , id_a , p_id } );
                  }
                  
                  // third neighbor step
                  uint64_t enc_atom_b = atom_from_idmap_if_found(id_b,id_map,id_map_ghosts,atom_not_found);
                  assert( enc_atom_b!=atom_not_found || ( !id_b_is_here && !cell_is_here ) );
                  if( enc_atom_b!=atom_not_found )
                  {
                    size_t cell_b=-1, pos_b=-1;
                    decode_cell_particle( enc_atom_b, cell_b, pos_b);
                    const auto& neigh_b = cells[cell_b][field_cmol][pos_b];
                    for(auto id_c:neigh_b) if(id_c!=id_a && id_c!=p_id)
                    {
                      assert( id_c != id_b );
                      if(id_c==std::numeric_limits<uint64_t>::max()) break;
                      bool id_c_is_here = ( id_map.find(id_c) != id_map.end() );
                      bool add_torsion = (cell_is_here || id_a_is_here || id_b_is_here || id_c_is_here) && (p_id!=id_b) && (id_c!=id_a) && (p_id!=id_c) ;
                      bool torsion_owner = (cell_is_here && !id_c_is_here) || (cell_is_here==id_c_is_here && p_id<id_c);
                      if( add_torsion && torsion_owner )
                      {
                        if( p_id < id_c ) torsions.push_back( ChemicalTorsion{ p_id , id_a , id_b , id_c } );
                        else              torsions.push_back( ChemicalTorsion{ id_c , id_b , id_a , p_id } );

                      }
                    }// for(auto id_c:neigh_b)
                  }// if enc_atom_b!=atom_not_found
                  
                }// for(auto id_b:neigh_a)
              }// if enc_atom_pos_a!=atom_not_found
              
            }

          }
        
        } // if p_is_leader         
      }
      GRID_OMP_FOR_END

      size_t bond_count = bonds.size();
      size_t angle_count = angles.size();
      size_t torsion_count = torsions.size();
      size_t improper_count = impropers.size();
      per_thread_counts[tid] = { bond_count , angle_count , torsion_count , improper_count };
      
#     pragma omp barrier
      
#     pragma omp single nowait
      {
        for(int i=1;i<nt;i++)
        {
          per_thread_counts[i].bond_count     += per_thread_counts[i-1].bond_count;
          per_thread_counts[i].angle_count    += per_thread_counts[i-1].angle_count;
          per_thread_counts[i].torsion_count  += per_thread_counts[i-1].torsion_count;
          per_thread_counts[i].improper_count += per_thread_counts[i-1].improper_count;
        }
        per_thread_bonds    [0].resize( per_thread_counts[nt-1].bond_count );
        per_thread_angles   [0].resize( per_thread_counts[nt-1].angle_count );
        per_thread_torsions [0].resize( per_thread_counts[nt-1].torsion_count );
        per_thread_impropers[0].resize( per_thread_counts[nt-1].improper_count );
      }

#     pragma omp barrier
      
      int bond_offset = 0;
      int angle_offset = 0;
      int torsion_offset = 0;
      int improper_offset = 0;
      if( tid > 0 )
      {
        bond_offset     = per_thread_counts[tid-1].bond_count;
        angle_offset    = per_thread_counts[tid-1].angle_count;
        torsion_offset  = per_thread_counts[tid-1].torsion_count;
        improper_offset = per_thread_counts[tid-1].improper_count;       
      }
      
      for(size_t i=0;i<bond_count    ;i++) per_thread_bonds    [0][bond_offset    +i] = per_thread_bonds    [tid][i];
      for(size_t i=0;i<angle_count   ;i++) per_thread_angles   [0][angle_offset   +i] = per_thread_angles   [tid][i];
      for(size_t i=0;i<torsion_count ;i++) per_thread_torsions [0][torsion_offset +i] = per_thread_torsions [tid][i];
      for(size_t i=0;i<improper_count;i++) per_thread_impropers[0][improper_offset+i] = per_thread_impropers[tid][i];
    }

    result_bonds     = std::move( per_thread_bonds    [0] );
    result_angles    = std::move( per_thread_angles   [0] );
    result_torsions  = std::move( per_thread_torsions [0] );
    result_impropers = std::move( per_thread_impropers[0] );

  } // end extract_connectivity(...)
  
  
} // namespace exaStamp

