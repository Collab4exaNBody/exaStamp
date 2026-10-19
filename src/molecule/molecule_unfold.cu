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

//#include <chrono>
#include <memory>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/molecule/id_map.h>
#include <exaStamp/molecule/intramolecular_image.h>
#include <exaStamp/molecule/molecule_min_atom_id.h>
#include <exaStamp/molecule/molecule_species.h>

#include <exanb/mpi/grid_update_ghosts.h>
#include <exanb/mpi/update_ghost_functors.h>
#include <mpi.h>

#include <exanb/core/grid_fields.h>
XNB_DECLARE_FIELD( exanb::Vec3d, molufpos, "molecule unfolded particle position" );

namespace exaStamp
{

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id >
    >
  class MoleculeUnfold : public OperatorNode
  {
    using UpdateGhostsScratch = typename UpdateGhostsUtils::UpdateGhostsScratch;

    ADD_SLOT( MPI_Comm                 , mpi               , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( long                     , mpi_tag           , INPUT , 0 );

    ADD_SLOT( Domain                   , domain            , INPUT );
    ADD_SLOT( GridT       , grid       , INPUT_OUTPUT );

    ADD_SLOT( IdMap       , id_map        , INPUT , REQUIRED );
    ADD_SLOT( IdMapGhosts , id_map_ghosts , INPUT , REQUIRED );

    ADD_SLOT( UpdateGhostConfig        , update_ghost_config , INPUT_OUTPUT , UpdateGhostConfig{} );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( MoleculeFoldedPositions  , molecule_folded_positions , INPUT_OUTPUT , DocString{"positions before unfolding, for molecule_restore_positions"} );
    ADD_SLOT( MoleculeSpeciesVector    , molecules         , INPUT_OUTPUT , OPTIONAL , DocString{"Molecule descriptions, keeps the unfolded anchor position of each molecule"} );
    ADD_SLOT( bool                     , continuous        , INPUT , true , DocString{"if true, each molecule is unfolded at the periodic image closest to its previous unfolding, so that successive unfolded snapshots form a continuous trajectory"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Unfolds molecules : bonded atoms are moved to their periodic image closest to their bond neighbors, so that each
molecule is whole. Unfolding starts from the anchor atom of each molecule (its smallest atom id).
With continuous: true (default), the anchor is placed at its periodic image closest to its position at the previous
unfolding (kept in molecules, written to restart files), so that molecules do not jump by a box vector from one
unfolded snapshot to the next. This requires dense molecule ids (molecule_numbering) and anchor atoms moving less
than half of the cell width between two unfoldings. The first unfolding keeps anchors inside the cell.
Use molecule_restore_positions after an output to restore folded positions.
)EOF";
    }

    inline void execute ()  override final
    {
      auto cells = grid->cells_accessor();
      size_t n_cells = grid->number_of_cells();
      const auto idmol_field = grid->field_accessor( field::idmol );
      const auto cmol_field = grid->field_accessor( field::cmol );
      static constexpr uint64_t null_id = std::numeric_limits<uint64_t>::max();
      static constexpr uint64_t null_loc = std::numeric_limits<uint64_t>::max();
      static constexpr size_t null_index = std::numeric_limits<size_t>::max();
      static constexpr GridCellValues * null_cell_values = nullptr;

      // basic assumption
      assert( domain->origin() == grid->origin() );

/*
      const auto ghost_update_fields = onika::make_flat_tuple( idmol_field ); 
      auto ghost_comm_buffers = UpdateGhostsUtils::make_forward_ghost_update_manager( cells, ghost_update_fields );
        grid_update_ghosts( exanb::ldbg, *mpi, *ghost_comm_scheme, *grid, *domain, nullptr,
                            ghost_comm_buffers, pecfunc, peqfunc, ghost_update_fields, *update_ghost_config );

*/

      // those two lambdas are used to launch parallel kernels within a function exterior to this OperatorNode's implementation
      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args ...); };
      auto peqfunc = [self=this]() -> onika::parallel::ParallelExecutionQueue& { return self->parallel_execution_queue(); };

      // map of molecule id (identical to owner's particle id) to owner particle location ( encoded cell / position in cell )
      std::unordered_map<uint64_t , uint64_t> molecule_owner;
      
      // maximum molecule bond distance found in dataset
      double bond_max_dist = 0.0;

      // save molecule ids, used as work storage below and restored at the end
      std::vector<uint64_t> saved_idmol;
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        const size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++) saved_idmol.push_back( cells[cell_i][idmol_field][p_i] );
      }

      // previous unfolded anchor positions, valid if molecules are densely numbered
      const IntramolecularImageBuilder image_builder( *domain );
      const size_t n_mol_instances = molecules.has_value() ? molecules->number_of_molecule_instances() : 0;
      const bool use_anchor_history = *continuous && n_mol_instances > 0;
      if( *continuous && ! use_anchor_history )
      {
        lout << "WARNING: molecule_unfold: no dense molecule numbering, continuous unfolding disabled (use molecule_numbering)" << std::endl;
      }
      if( use_anchor_history && ( molecules->m_unfold_anchor_id.size() != n_mol_instances || molecules->m_unfold_anchor_pos.size() != n_mol_instances ) )
      {
        molecules->m_unfold_anchor_id.assign( n_mol_instances , 0 );
        molecules->m_unfold_anchor_pos.assign( n_mol_instances , Vec3d{0.,0.,0.} );
      }
      std::vector<uint64_t> anchor_id( use_anchor_history ? n_mol_instances : 0 , 0 );
      std::vector<Vec3d> anchor_pos( use_anchor_history ? n_mol_instances : 0 , Vec3d{0.,0.,0.} );
      size_t n_continued = 0;

      // first pass : molecule id of each particle = smallest atom id of its molecule
      long update_count = 0;
      int pass = 0;
      propagate_molecule_min_atom_id( ldbg, *mpi, *ghost_comm_scheme, grid.get_pointer(), *domain, *update_ghost_config,
                                      pecfunc, peqfunc, *id_map, *id_map_ghosts );

      // this temporary position field holds particle's unfolded position with respect
      // to molecule owner particle.
      const auto ufpos_field = grid->field_accessor( field::molufpos );

      // a first pre-process pass initialize molufpos to unmodified particle's postiotion,
      // except for anchors of continuously unfolded molecules, placed at their image closest to their previous position
      size_t k_idmol = 0;
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        const auto * __restrict__ ids = cells[cell_i][field::id];
        const auto * __restrict__ rx = cells[cell_i][field::rx];
        const auto * __restrict__ ry = cells[cell_i][field::ry];
        const auto * __restrict__ rz = cells[cell_i][field::rz];
        auto * __restrict__ mol_ids = cells[cell_i][idmol_field];
        auto * __restrict__ ufpos = cells[cell_i][ufpos_field];
        const size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++,k_idmol++)
        {
          ufpos[p_i] = Vec3d{ rx[p_i] , ry[p_i] , rz[p_i] };
          const auto mid = mol_ids[p_i];
          if( ids[p_i] == mid )
          {
            if( use_anchor_history )
            {
              const uint64_t idmol = saved_idmol[k_idmol];
              if( idmol >= n_mol_instances )
              {
                fatal_error() << "molecule_unfold: molecule id "<<idmol<<" out of range (number of molecules = "<<n_mol_instances<<")" << std::endl;
              }
              if( molecules->m_unfold_anchor_id[idmol] == ids[p_i] + 1 )
              {
                const Vec3d prev = molecules->m_unfold_anchor_pos[idmol];
                ufpos[p_i] = prev + image_builder.closest_delta( ufpos[p_i] - prev );
                ++ n_continued;
              }
              anchor_id[idmol] = ids[p_i] + 1;
              anchor_pos[idmol] = ufpos[p_i];
            }
            const auto ploc = encode_cell_particle(cell_i,p_i,0);
            if( molecule_owner.find(mid) != molecule_owner.end() )
            {
              fatal_error() << "duplicate owner particle for molecule id #"<<mol_ids[p_i]<<std::endl;
            }
            molecule_owner.insert( { mid , ploc } );
            // ldbg << "cell #"<<cell_i<<" particle #"<<p_i<<" owner of molid #"<<mid<<std::endl;
          }
          else
          {
            mol_ids[p_i] = null_id;
          }
        }
      }

      // bonded atoms are unfolded to their closest periodic image (physical distance, works for triclinic cells),
      // any number of box lengths away from the already unfolded position
      const double max_bond_dist = image_builder.max_unique_distance();
      ldbg << "Number of molecules = "<< molecule_owner.size()<<" , max. unique bond distance = "<<max_bond_dist<<std::endl;

      // debug map to keep track of duplicated positions
      std::unordered_map<uint64_t,Vec3d> particle_pos_map;

      // second multipass algorithm : unfold positions by applying boundary conditions shift
      // whenever a bound particle has an excessive distance, taking the one with the
      pass = 0;
      do
      {                  
        int sub_pass = 0;        
        do
        {
          update_count = 0;
          for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
          {
            const auto *  __restrict__ cmol = cells[cell_i][cmol_field];
            auto * __restrict__ mol_ids = cells[cell_i][idmol_field];
            auto * __restrict__ ufpos = cells[cell_i][ufpos_field];
            size_t n = cells[cell_i].size();
            for(size_t p_i=0;p_i<n;p_i++)
            {
              if( mol_ids[p_i] != null_id ) // unfolded central particle
              {
                const Vec3d ri = ufpos[p_i];
                for(size_t j=0;j<cmol[p_i].size();j++)
                {
                  const auto conloc = atom_from_idmap_if_found( cmol[p_i][j] , *id_map , *id_map_ghosts , null_loc );
                  if( conloc != null_loc )
                  {
                    size_t cell=null_index, pos=null_index;
                    decode_cell_particle( conloc , cell, pos );
                    assert( cell!=null_index && pos!=null_index );
                    assert( cell!=cell_i || pos!=p_i );
                    const Vec3d rj = cells[cell][ufpos_field][pos] ;
                    if( cells[cell][idmol_field][pos] == null_id ) // not unfolded bond neighbor
                    {
                      const Vec3d rij = image_builder.closest_delta( rj - ri );
                      const double norm_rij = image_builder.physical_norm(rij);
                      if( norm_rij >= max_bond_dist || norm_rij == 0.0 )
                      {
                        fatal_error() << "pre-check: in C#"<<cell_i<<"P#"<<p_i<<" (id="<<cells[cell_i][field::id][p_i]<<") bond#"<<j
                        <<" with C#"<<cell<<"P#"<<pos<<" (id="<<cells[cell][field::id][pos]<<",ghost="<<grid->is_ghost_cell(cell)<<") bad dist. "<<norm_rij<<" not in ] 0 ; "<<max_bond_dist<<" ["
                        <<" , ri="<<ri<<" , rj="<<rj<< std::endl;
                      }
                      bond_max_dist = std::max( bond_max_dist , norm_rij );   
                      const Vec3d rj_unfolded = ri + rij;
                      auto it = particle_pos_map.find( cells[cell][field::id][pos] );
                      if( it != particle_pos_map.end() )
                      {
                        if( rj_unfolded != it->second )
                        {
                          fatal_error() << "different ghost pos : rj_unfolded="<<rj_unfolded<<" , stored="<<it->second<<" dist="<<norm(rj_unfolded-it->second)<<std::endl;
                        }
                      }
                      else
                      {
                        particle_pos_map.insert( { cells[cell][field::id][pos] , rj_unfolded } );
                      }
                      cells[cell][ufpos_field][pos] = rj_unfolded;
                      cells[cell][idmol_field][pos] = mol_ids[p_i];
                      ++ update_count;
                    } // if bond neighbor nor unfolded yet
                  } // bond neighbor is present
                } // for each bond neighbor
              } // if central is unfolded
            } // for each cell's particle
          } // for each non ghost cell

          ++ sub_pass;
          ldbg << "molecule_unfold: unfold central positions, sub_pass="<<sub_pass<<" , updates="<<update_count <<std::endl;
        } while( update_count > 0 );

        // propagate ghost updates
        const auto ghost_update_fields = onika::make_flat_tuple( idmol_field, ufpos_field );
        auto ghost_comm_buffers = UpdateGhostsUtils::make_forward_ghost_update_manager( cells, ghost_update_fields );
        grid_update_ghosts( exanb::ldbg, *mpi, *ghost_comm_scheme, grid.get_pointer(), *domain, null_cell_values,
                            ghost_comm_buffers, pecfunc, peqfunc, ghost_update_fields, *update_ghost_config );      

        // gather updates from ghosts (from neighbor sub-domains)
        update_count = 0;
        for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
        {
          const auto *  __restrict__ cmol = cells[cell_i][cmol_field];
          auto * __restrict__ mol_ids = cells[cell_i][idmol_field];
          auto * __restrict__ ufpos = cells[cell_i][ufpos_field];
          size_t n = cells[cell_i].size();
          for(size_t p_i=0;p_i<n;p_i++)
          {
            if( mol_ids[p_i] == null_id ) // central particle not unfolded yet
            {
              const Vec3d ri = ufpos[p_i];
              int ri_updates = 0;
              Vec3d last_ri_unfolded = { 0, 0, 0 };
              for(size_t j=0;j<cmol[p_i].size();j++)
              {
                const auto conloc = atom_from_idmap_if_found( cmol[p_i][j] , *id_map , *id_map_ghosts , null_loc );
                if( conloc != null_loc )
                {
                  size_t cell=null_index, pos=null_index;
                  decode_cell_particle( conloc , cell, pos );
                  assert( cell!=null_index && pos!=null_index );
                  assert( cell!=cell_i || pos!=p_i );
                  const Vec3d rj = cells[cell][ufpos_field][pos] ;
                  if( cells[cell][idmol_field][pos] != null_id ) // bond neighbor already unfolded
                  {
                    const Vec3d rij = image_builder.closest_delta( rj - ri );
                    const Vec3d ri_unfolded = rj - rij;
                    if( ri_updates > 0 )
                    {
                      if( ri_unfolded != last_ri_unfolded )
                      {
                        fatal_error() << "different update of central particle : last_ri_unfolded="<<last_ri_unfolded<<" ri_unfolded="<<ri_unfolded<<" dist="<<norm(ri_unfolded-last_ri_unfolded)<<" ri_updates="<<ri_updates<<" j="<<j<<std::endl;
                      }
                    }
                    else
                    {
                      last_ri_unfolded = ri_unfolded;
                    }
                    ufpos[p_i] = ri_unfolded;
                    mol_ids[p_i] = cells[cell][idmol_field][pos];
                    ++ ri_updates;
                  } // if bond neighbor nor unfolded yet
                } // bond neighbor is present
              } // for each bond neighbor
              update_count += ri_updates;
            } // if central is unfolded
          } // for each cell's particle
        } // for each non ghost cell

        MPI_Allreduce( MPI_IN_PLACE, &update_count , 1 , MPI_LONG , MPI_SUM , *mpi );
        ++ pass;
        ldbg << "molecule_unfold: unfold from bond neighbor positions, pass="<<pass<<" , updates="<<update_count <<std::endl;
        
      } while( update_count > 0 );

      // keep folded positions, so that molecule_restore_positions can undo the unfolding (e.g. after an output)
      molecule_folded_positions->m_rx.clear();
      molecule_folded_positions->m_ry.clear();
      molecule_folded_positions->m_rz.clear();
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        const size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++)
        {
          molecule_folded_positions->m_rx.push_back( cells[cell_i][field::rx][p_i] );
          molecule_folded_positions->m_ry.push_back( cells[cell_i][field::ry][p_i] );
          molecule_folded_positions->m_rz.push_back( cells[cell_i][field::rz][p_i] );
        }
      }

      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        const auto * __restrict__ ufpos = cells[cell_i][ufpos_field];
        const auto * __restrict__ mol_ids = cells[cell_i][idmol_field];
        auto * __restrict__ rx   = cells[cell_i][field::rx];
        auto * __restrict__ ry   = cells[cell_i][field::ry];
        auto * __restrict__ rz   = cells[cell_i][field::rz];
        const size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++)
        {
          if( mol_ids[p_i] != null_id )
          {
            rx[p_i] = ufpos[p_i].x;
            ry[p_i] = ufpos[p_i].y;
            rz[p_i] = ufpos[p_i].z;
          }
          else
          {
            fatal_error() << "Cell #"<<cell_i<<" , particle #"<<p_i<<" not unfolded"<<std::endl;
          }
        }
      }

      grid->remove_flat_array( "molufpos" );

      // restore molecule ids
      size_t k = 0;
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        const size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++) cells[cell_i][idmol_field][p_i] = saved_idmol[k++];
      }
      
      // each anchor is owned by a single process : gather unfolded anchor positions of all molecules
      static_assert( sizeof(Vec3d) == 3 * sizeof(double) , "Vec3d reduced as 3 doubles" );
      if( use_anchor_history )
      {
        MPI_Allreduce( MPI_IN_PLACE, anchor_id.data() , n_mol_instances , MPI_UINT64_T , MPI_MAX , *mpi );
        MPI_Allreduce( MPI_IN_PLACE, anchor_pos.data() , 3 * n_mol_instances , MPI_DOUBLE , MPI_SUM , *mpi );
        MPI_Allreduce( MPI_IN_PLACE, &n_continued , 1 , MPI_UNSIGNED_LONG , MPI_SUM , *mpi );
        molecules->m_unfold_anchor_id = std::move( anchor_id );
        molecules->m_unfold_anchor_pos = std::move( anchor_pos );
        ldbg << "molecule_unfold: "<<n_continued<<" / "<<n_mol_instances<<" molecules unfolded continuously" << std::endl;
      }

      MPI_Allreduce( MPI_IN_PLACE, &bond_max_dist , 1 , MPI_DOUBLE , MPI_MAX , *mpi );
      ldbg << "bond_max_dist = " << bond_max_dist << std::endl;
    }

  };

  template<class GridT> using MoleculeUnfoldTmpl = MoleculeUnfold<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_unfold)
  {
   OperatorNodeFactory::instance()->register_factory( "molecule_unfold", make_grid_variant_operator<MoleculeUnfoldTmpl> );
  }

}
