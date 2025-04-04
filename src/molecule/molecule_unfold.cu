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
#include <exaStamp/molecule/periodic_r_delta.h>

#include <mpi.h>
#include <exanb/mpi/grid_update_ghosts.h>

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

    ADD_SLOT( UpdateGhostsScratch      , ghost_comm_buffers, PRIVATE );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme , INPUT_OUTPUT , REQUIRED );

  public:
    inline void execute ()  override final
    {
      auto cells = grid->cells_accessor();
      size_t n_cells = grid->number_of_cells();
      const auto idmol_field = grid->field_accessor( field::idmol );
      const auto cmol_field = grid->field_accessor( field::cmol );
      static constexpr uint64_t null_id = std::numeric_limits<uint64_t>::max();
      static constexpr uint64_t null_loc = std::numeric_limits<uint64_t>::max();
      static constexpr size_t null_index = std::numeric_limits<size_t>::max();

      // basic assumption
      assert( domain->origin() == grid->origin() );

      // those two lambdas are used to launch parallel kernels within a function exterior to this OperatorNode's implementation
      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args ...); };
      //      auto pesfunc = [self=this](unsigned int i) { return self->parallel_execution_stream(i); };
      auto peqfunc = [self=this]() -> onika::parallel::ParallelExecutionQueue& { return self->parallel_execution_queue(); };

      // map of molecule id (identical to owner's particle id) to owner particle location ( encoded cell / position in cell )
      std::unordered_map<uint64_t , uint64_t> molecule_owner;
      
      // maximum molecule bond distance found in dataset
      double bond_max_dist = 0.0;

      // preamble : assign particle's molid particles' ids, so that we start with individual molecules
      // of one atom each, and we will then merge those bound by intramolecular bond
      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++) cells[cell_i][idmol_field][p_i] = cells[cell_i][field::id][p_i];
      }
      ldbg << "molecule_unfold: id_map @"<<id_map.get_pointer()<<" size = "<<id_map->size()<<" , id_map_ghosts @ "<<id_map_ghosts.get_pointer() <<" size = "<<id_map_ghosts->size() <<std::endl;

      // first multi pass algorithm : propagate molecule ids through bonds, keeping only the minimum id a the molecule id
      long update_count = 0;
      int pass = 0;
      do
      {
        const auto ghost_update_fields = onika::make_flat_tuple( idmol_field ); 
        grid_update_ghosts( exanb::ldbg, *mpi, *ghost_comm_scheme, *grid, *domain, nullptr,
                            *ghost_comm_buffers, pecfunc, peqfunc, ghost_update_fields,
                            *mpi_tag, true, true, true, true, false, std::false_type{} );
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
              const auto conloc = atom_from_idmap_if_found( cmol[p_i][j] , *id_map , *id_map_ghosts , null_loc );
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
        MPI_Allreduce( MPI_IN_PLACE, &update_count , 1 , MPI_LONG , MPI_SUM , *mpi );
        ++ pass;
        ldbg << "molecule_unfold: connect molecules, pass="<<pass<<" , updates="<<update_count<<std::endl;
      } while( update_count > 0 );

      // this temporary position field holds particle's unfolded position with respect
      // to molecule owner particle.
      const auto ufpos_field = grid->field_accessor( field::molufpos );

      // a first pre-process pass initialize molufpos to unmodified particle's postiotion
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        const auto * __restrict__ ids = cells[cell_i][field::id];
        const auto * __restrict__ rx = cells[cell_i][field::rx];
        const auto * __restrict__ ry = cells[cell_i][field::ry];
        const auto * __restrict__ rz = cells[cell_i][field::rz];
        auto * __restrict__ mol_ids = cells[cell_i][idmol_field];
        auto * __restrict__ ufpos = cells[cell_i][ufpos_field];
        const size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++)
        {
          ufpos[p_i] = Vec3d{ rx[p_i] , ry[p_i] , rz[p_i] };
          const auto mid = mol_ids[p_i];
          if( ids[p_i] == mid )
          {
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

      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 
      ldbg << "Number of molecules = "<< molecule_owner.size()<<" , size_box = "<<size_box<<" , half_min_size_box = "<<half_min_size_box<<std::endl;

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
                      const Vec3d rij = periodic_r_delta_loop( ri , rj , size_box , half_min_size_box );
                      const double norm_rij = norm(rij);
                      if( norm_rij > half_min_size_box || norm_rij == 0.0 )
                      {
                        fatal_error() << "pre-check: in C#"<<cell_i<<"P#"<<p_i<<" (id="<<cells[cell_i][field::id][p_i]<<") bond#"<<j
                        <<" with C#"<<cell<<"P#"<<pos<<" (id="<<cells[cell][field::id][pos]<<",ghost="<<grid->is_ghost_cell(cell)<<") bad dist. "<<norm_rij<<" not in ] 0 ; "<<half_min_size_box<<" ]"
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
        grid_update_ghosts( exanb::ldbg, *mpi, *ghost_comm_scheme, *grid, *domain, nullptr,
                            *ghost_comm_buffers, pecfunc,peqfunc, ghost_update_fields,
                            *mpi_tag, true, true, true, true, false, std::false_type{} );        

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
                    const Vec3d rij = periodic_r_delta_loop( ri , rj , size_box , half_min_size_box );
                    const Vec3d ri_unfolded = ri - rij;
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
