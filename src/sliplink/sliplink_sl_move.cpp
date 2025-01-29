#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/parallel/random.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/particle_id_translation.h>
#include <onika/memory/allocator.h> // for ONIKA_ASSUME_ALIGNED macro
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/algorithm.h>
#include <exaStamp/sliplink/sliplink.h>

#include <onika/mpi/data_types.h>

#include <random>
#include <mpi.h>

#include <vector>

namespace exaStamp
{

  // apply SL<->anchor spring force and friction,
  // directly accounted to SL's position, as sliplinks method is only valid in an overdemt scheme
  template<
    class GridT
//    , class = AssertGridHasFields< GridT, field::_id > // necessary only in debug mode
    >
  class SlipLinkSLMoveOperator : public OperatorNode
  {
    ADD_SLOT(MPI_Comm           , mpi             , INPUT        , MPI_COMM_WORLD );
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT        , REQUIRED); 
    ADD_SLOT(double             , bond_max_stretch , INPUT , 0.5 ); // fraction of bond_max_dist.
    ADD_SLOT(SLGrid             , sl_grid         , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT(ParticleIdMap      , id_map          , INPUT );
    ADD_SLOT(GridT              , grid            , INPUT        , REQUIRED );
    ADD_SLOT(long               , sl_regen_count  , INPUT_OUTPUT );

  public:

    inline void execute () override final
    {        
      MPI_Comm comm = *mpi;
      int nprocs = 1;
      int rank = 0;
      MPI_Comm_size(comm,&nprocs);
      MPI_Comm_rank(comm,&rank);

      if( ! m_SLRegenInfo_mpitype_initialized )
      {
        MPI_Type_contiguous( sizeof(SLRegenInfo), MPI_BYTE, &m_SLRegenInfo_mpitype );
        MPI_Type_commit( &m_SLRegenInfo_mpitype );
        m_SLRegenInfo_mpitype_initialized = true;
      }

      auto cells = grid->cells();
      auto sl_cells = sl_grid->cells();

      assert( grid->number_of_cells() == sl_grid->number_of_cells() );  
          
      IJK dims = grid->dimension();
      assert( dims == sl_grid->dimension() );

      assert( grid->ghost_layers() == sl_grid->ghost_layers() );
      const ssize_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);
 
      const size_t n_chains = sliplink_config->number_of_chains;
      const size_t n_beads = sliplink_config->beads_per_chain;
       
      const double cte4_h = sliplink_config->cte4 * sliplink_config->h;
      const double sigma4 = sliplink_config->sigma4;
      const double bond_max_search_dist = sliplink_config->bond_max_dist * ( 1. + *bond_max_stretch );

      std::map<uint64_t,SLRegenInfoNoId> update_sl_map;

      size_t n_slipped = 0;
      size_t n_regen = 0;

#     pragma omp parallel
      {
        auto& re = onika::parallel::random_engine();
        // random distributions used for sliplink parameters
        std::uniform_int_distribution<size_t> random_chain( 0 , n_chains-1 );
        std::uniform_int_distribution<size_t> random_bead( 0 , n_beads-2 );
        std::uniform_real_distribution<double> uniform_01( 0.0 , 1.0 );
        std::normal_distribution<double> gaussian_friction(0.0, sigma4);

        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc, reduction(+:n_slipped,n_regen) )
        {
          size_t i = grid_ijk_to_index( dims , loc + ghost_layers);
          assert( ! sl_grid->is_ghost_cell(i) );

          double * __restrict__ rx = sl_cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
          double * __restrict__ ry = sl_cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
          double * __restrict__ rz = sl_cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
          SlipLinkField * __restrict__ sl = sl_cells[i][field::sl]; ONIKA_ASSUME_ALIGNED(sl);
          size_t n = sl_cells[i].size();

          for(size_t j=0;j<n;j++)
          {
            uint64_t left_bead_id = sl[j].left_bead_id;
            uint64_t right_bead_id = left_bead_id+1;
            assert( is_particle_id_valid(left_bead_id) && is_particle_id_valid(right_bead_id) );
            assert( sliplink_left_bead_id(right_bead_id,n_beads) == left_bead_id );
            assert( sliplink_right_bead_id(left_bead_id,n_beads) == right_bead_id );

            Vec3d old_r = { rx[j], ry[j], rz[j] };
                        
            // left bead
            uint64_t local_left_bead = global_to_nearest_local_id( left_bead_id, *id_map, *grid, old_r, bond_max_search_dist );
#           ifndef NDEBUG
#           pragma omp critical
            { assert( nearest_local_id_check(lerr,left_bead_id,local_left_bead,*id_map, *grid, old_r, bond_max_search_dist) ); }
#           endif
            size_t lc = 0;
            size_t lp = 0;
            decode_cell_particle( local_left_bead , lc , lp );
            assert( grid->is_valid_cell_particle(lc,lp) );
            Vec3d lr = { cells[lc][field::rx][lp] , cells[lc][field::ry][lp], cells[lc][field::rz][lp] };

            // right bead
            uint64_t local_right_bead = global_to_nearest_local_id( right_bead_id, *id_map, *grid, old_r, bond_max_search_dist );
#           ifndef NDEBUG
#           pragma omp critical
            { assert( nearest_local_id_check(lerr,right_bead_id,local_right_bead,*id_map, *grid, old_r, bond_max_search_dist) ); }
#           endif
            size_t rc = 0;
            size_t rp = 0;
            decode_cell_particle( local_right_bead , rc , rp );
            assert( grid->is_valid_cell_particle(rc,rp) );
            Vec3d rr = { cells[rc][field::rx][rp] , cells[rc][field::ry][rp], cells[rc][field::rz][rp] };

            // SL move
            double slip = cte4_h * dot( rr-lr , sl[j].anchor_displ );
            double friction = gaussian_friction(re);
            //slip=0.0; friction=0.0; // FIXME: REMOVE THIS LINE
            double xj_f = sl[j].xj_frac + slip + friction;

            if( xj_f < -1.0 || xj_f > (2.0-1.e-6) )
            {
              lerr << "sliplink over slip : xj_frac = "<<xj_f << std::endl;
              xj_f = clamp( xj_f , -1.0 , 2.0-1.e-6 );
            }

            bool local_update = true;
            bool full_regen = false;
            uint64_t new_left_bead_id = left_bead_id;
            uint64_t new_right_bead_id = right_bead_id;

            // shift left in the chain
            if( xj_f<0.0 || xj_f>=1.0 )
            {
              if( xj_f<0.0 )
              {           
                xj_f += 1.0;
                new_left_bead_id = sliplink_left_bead_id( left_bead_id , sliplink_config->beads_per_chain );
              }
              else
              {
                xj_f -= 1.0;
                new_left_bead_id = right_bead_id;
              }
              assert( xj_f >= 0.0 && xj_f <= 1.0 );
              
              new_right_bead_id = sliplink_right_bead_id( new_left_bead_id , sliplink_config->beads_per_chain );
              
              if( new_left_bead_id != PARTICLE_NO_ID && new_right_bead_id != PARTICLE_NO_ID )
              {                
                if( new_right_bead_id == left_bead_id ) // left shitfed
                {
                  uint64_t new_local_left_bead = global_to_nearest_local_id( new_left_bead_id, *id_map, *grid, lr, bond_max_search_dist );
                  if( is_particle_id_valid(new_local_left_bead) )
                  {
                    rc=lc; rp=lp; rr=lr;
                    decode_cell_particle( new_local_left_bead , lc , lp );
                    assert( grid->is_valid_cell_particle(lc,lp) );
                    lr = Vec3d{ cells[lc][field::rx][lp] , cells[lc][field::ry][lp], cells[lc][field::rz][lp] };
                  }
                  else { local_update = false; }
                }
                else // right sifhted
                {
                  assert( new_left_bead_id == right_bead_id );
                  uint64_t new_local_right_bead = global_to_nearest_local_id( new_right_bead_id, *id_map, *grid, rr, bond_max_search_dist );
                  if( is_particle_id_valid(new_local_right_bead) )
                  {
                    lc=rc; lp=rp; lr=rr;
                    decode_cell_particle( new_local_right_bead , rc , rp );
                    assert( grid->is_valid_cell_particle(rc,rp) );
                    rr = { cells[rc][field::rx][rp] , cells[rc][field::ry][rp], cells[rc][field::rz][rp] };
                  }
                  else { local_update = false; }                  
                }
                ++ n_slipped;
              }
              else
              {
                full_regen = true;
                local_update = false;
                ++ n_regen;
              }
            }
            
            if( local_update )  // simple local update
            {
              assert( is_particle_id_valid(new_left_bead_id) );
              sl[j].xj_frac = xj_f;
              sl[j].left_bead_id = new_left_bead_id;
              Vec3d new_r = (1.0-xj_f) * lr + xj_f * rr;
              rx[j] = new_r.x;
              ry[j] = new_r.y;
              rz[j] = new_r.z;
              sl[j].anchor_displ += old_r - new_r;             
            }
            else if( full_regen ) // SL and SL peer regeneration
            {               
              uint64_t sl_id = sl_cells[i][field::id][j];
              uint64_t sl_peer_id = sl_id^1ull;
              
              uint64_t bead_id_1 = random_chain(re)*n_beads + random_bead(re);
              double xj_frac_1 = uniform_01(re);

              uint64_t bead_id_2 = random_chain(re)*n_beads + random_bead(re);
              double xj_frac_2 = uniform_01(re);

#             pragma omp critical
              {
                ldbg << "invalidate SL #" << sl_id << " and #" << sl_peer_id << std::endl << std::flush;
                update_sl_map.insert( std::pair<uint64_t,SLRegenInfoNoId>(sl_id,SLRegenInfoNoId{bead_id_1,xj_frac_1,true}) );
                update_sl_map.insert( std::pair<uint64_t,SLRegenInfoNoId>(sl_peer_id,SLRegenInfoNoId{bead_id_2,xj_frac_2,true}) );
              }
            }
            else // non local slip
            {
              uint64_t sl_id = sl_cells[i][field::id][j];
#             pragma omp critical
              {
                auto it = update_sl_map.find( sl_id );
                if( it == update_sl_map.end() ) // only if not fully regenerated
                {
                  ldbg << "invalidate SL #" << sl_id << std::endl << std::flush;
                  update_sl_map.insert( std::pair<uint64_t,SLRegenInfoNoId>(sl_id,SLRegenInfoNoId{new_left_bead_id,xj_f,false}) );
                }
              }
            }

          }

        }
        GRID_OMP_FOR_END
      }
        
      // ********************************************
      // *** 2. regenerate destroyed SL and peers ***
      // ***    and update non local slips        ***
      // ********************************************

      ldbg << "n_slipped="<< n_slipped<<", n_regen="<<n_regen<<std::endl;
      *sl_regen_count = n_regen; // first approximation

      // exchange all SL destroy accross all processors so that peer SLs
      std::vector<SLRegenInfo> regen_sl;
      regen_sl.reserve(update_sl_map.size());
      for(const auto& p : update_sl_map)
      {
        regen_sl.push_back( SLRegenInfo{ p.first, p.second.left_bead_id, p.second.xj_frac, p.second.full_regen } );
      }
      update_sl_map.clear();

      int regen_count = regen_sl.size();
      int all_regen_count[nprocs];
      int all_regen_displs[nprocs];
      MPI_Allgather( &regen_count, 1, MPI_INT, all_regen_count, 1, MPI_INT, comm );
      assert( all_regen_count[rank] == regen_count );

      size_t total_regen_count = 0;
      for(int i=0;i<nprocs;i++)
      {
        all_regen_displs[i] = total_regen_count;
        total_regen_count += all_regen_count[i];
      }
            
      if( total_regen_count > 0 )
      {
        ldbg << "slmove: total_regen_count="<< total_regen_count<<std::endl;

        std::vector<SLRegenInfo> all_regen_sl(total_regen_count);
        MPI_Allgatherv( regen_sl.data(), regen_count , m_SLRegenInfo_mpitype, all_regen_sl.data(), all_regen_count, all_regen_displs, m_SLRegenInfo_mpitype, comm );
        regen_sl.clear();

        // remove duplicate SL ids, prioritizing destroyed SLs over slipped SLs
        update_sl_map.clear();
        for(const auto& sl : all_regen_sl)
        {
          auto it = update_sl_map.find( sl.sl_id );
          if( it == update_sl_map.end() )
          {
            update_sl_map[ sl.sl_id ] = SLRegenInfoNoId { sl.left_bead_id, sl.xj_frac, sl.full_regen };
          }
          else if( ! it->second.full_regen )
          {
            it->second = SLRegenInfoNoId { sl.left_bead_id, sl.xj_frac, sl.full_regen };
          }
        }
        *sl_regen_count = update_sl_map.size();
        ldbg << "output sl_regen_count = " << *sl_regen_count << std::endl;

        all_regen_sl.clear();
        for(const auto& p : update_sl_map)
        {
          all_regen_sl.push_back( SLRegenInfo{ p.first, p.second.left_bead_id, p.second.xj_frac, p.second.full_regen } );
        }
        update_sl_map.clear();
        total_regen_count = all_regen_sl.size();
        ldbg << "slmove: reduced total_regen_count="<< total_regen_count<<std::endl;

        size_t n_resolved_positions = 0;
        std::vector<double> sl_positions( total_regen_count*3 , std::numeric_limits<double>::lowest() );
        for(size_t i=0;i<total_regen_count;i++)
        {
          assert( is_particle_id_valid(all_regen_sl[i].left_bead_id) );
          uint64_t left_local_id = global_to_own_local_id( all_regen_sl[i].left_bead_id, *id_map, *grid );
          if( is_particle_id_valid(left_local_id) )
          {
            size_t lc=0,lp=0;
            decode_cell_particle( left_local_id, lc, lp );
            assert( grid->is_valid_cell_particle(lc,lp) );
            Vec3d lr = { cells[lc][field::rx][lp] , cells[lc][field::ry][lp], cells[lc][field::rz][lp] };
            
            uint64_t right_local_id = global_to_nearest_local_id( all_regen_sl[i].left_bead_id+1, *id_map, *grid, lr, bond_max_search_dist );
            assert( is_particle_id_valid(right_local_id) );
            size_t rc=0,rp=0;
            decode_cell_particle( right_local_id, rc, rp );
            assert( grid->is_valid_cell_particle(rc,rp) );
            Vec3d rr = { cells[rc][field::rx][rp] , cells[rc][field::ry][rp], cells[rc][field::rz][rp] };
            
            double xj_f = all_regen_sl[i].xj_frac;
            Vec3d slpos = (1.0-xj_f) * lr + xj_f * rr;
            sl_positions[i*3+0] = slpos.x;
            sl_positions[i*3+1] = slpos.y;
            sl_positions[i*3+2] = slpos.z;
            ++ n_resolved_positions;
          }
        }
        MPI_Allreduce(MPI_IN_PLACE,sl_positions.data(),total_regen_count*3,MPI_DOUBLE,MPI_MAX,comm);

#       ifndef NDEBUG
        MPI_Allreduce(MPI_IN_PLACE,&n_resolved_positions,1,exanb::mpi_datatype<size_t>(),MPI_SUM,comm);
        ldbg << "total resolved positions = "<< n_resolved_positions << std::endl;
        assert( n_resolved_positions == total_regen_count );
        for(size_t i=0;i<total_regen_count;i++)
        {
          assert( sl_positions[i*3+0] != std::numeric_limits<double>::lowest() );
          assert( sl_positions[i*3+1] != std::numeric_limits<double>::lowest() );
          assert( sl_positions[i*3+2] != std::numeric_limits<double>::lowest() );
        }
#       endif

        std::unordered_map<uint64_t,SLRegenInfo2> all_regen_sl_map;
        for(size_t i=0;i<total_regen_count;i++)
        {
          all_regen_sl_map[ all_regen_sl[i].sl_id ] = SLRegenInfo2{ all_regen_sl[i].left_bead_id, Vec3d{sl_positions[i*3+0],sl_positions[i*3+1],sl_positions[i*3+2]} , all_regen_sl[i].xj_frac, all_regen_sl[i].full_regen };
        }
        sl_positions.clear();
        all_regen_sl.clear();

        size_t n_regen_sl = all_regen_sl_map.size();
        ldbg << "regenerating "<< n_regen_sl << " sliplinks"<<std::endl;

        // random distributions for sliplinks' anchors
        std::normal_distribution<double> gaussian_anchor_displ( 0.0, sliplink_config->sigma3 );

#       pragma omp parallel
        {
          auto& re = onika::parallel::random_engine();
          GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghost)
          {
            IJK loc = loc_no_ghost + ghost_layers;
            size_t i = grid_ijk_to_index(dims,loc);

            double * __restrict__ rx = sl_cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
            double * __restrict__ ry = sl_cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
            double * __restrict__ rz = sl_cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
            SlipLinkField * __restrict__ sl = sl_cells[i][field::sl]; ONIKA_ASSUME_ALIGNED(sl);
            uint64_t * __restrict__ sl_ids = sl_cells[i][field::id]; ONIKA_ASSUME_ALIGNED(sl_ids);
            size_t n = sl_cells[i].size();
            for(size_t j=0;j<n;j++)
            {
              auto sl_it = all_regen_sl_map.find( sl_ids[j] );
              if( sl_it != all_regen_sl_map.end() )
              {
                sl[j].left_bead_id = sl_it->second.bead_id;
                sl[j].xj_frac = sl_it->second.frac;
                if( sl_it->second.full_regen )
                {
                  // sliplink anchor point follow a gaussian distribution around the sliplink position
                  sl[j].anchor_displ = Vec3d { gaussian_anchor_displ(re),
                                               gaussian_anchor_displ(re),
                                               gaussian_anchor_displ(re) };
                }
                Vec3d r = sl_it->second.pos;
                rx[j] = r.x; ry[j] = r.y; rz[j] = r.z;
              }
            }
          }
          GRID_OMP_FOR_END
        }
      }

    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return
R"EOF(
Written by Michelin & CEA/DIF
Credits to Claire Lemarchand, Ioanis Tanis, Thierry Carrard

computes sliplink movement along curvilinear abcsissa. movement has two origins :
  - spring force between anchor and sliplink
  - sliplink friction
)EOF";
    }

  private:
    MPI_Datatype m_SLRegenInfo_mpitype;
    bool m_SLRegenInfo_mpitype_initialized = false;
  };

  template<class GridT> using SlipLinkSLMoveOperatorTmpl = SlipLinkSLMoveOperator<GridT>;
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(sliplink_sl_move)
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_sl_move", make_grid_variant_operator< SlipLinkSLMoveOperatorTmpl > );
  }

}

