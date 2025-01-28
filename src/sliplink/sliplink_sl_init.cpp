#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/parallel/random.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/particle_id_constants.h>
#include <exaStamp/sliplink/sliplink.h>
#include <exanb/core/particle_id_translation.h>
#include <exanb/core/integer_range_algorithm.h>

#include <mpi.h>

#include <vector>
#include <algorithm>
#include <random>


namespace exaStamp
{

  template<class GridT>
  class SlipLinkInitOperator : public OperatorNode
  {
    ADD_SLOT(MPI_Comm           , mpi             , INPUT , MPI_COMM_WORLD );
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT, REQUIRED ); 
    ADD_SLOT(double             , bond_max_stretch , INPUT , 0.5 ); // fraction of bond_max_dist.
    ADD_SLOT(Domain             , domain          , INPUT, REQUIRED );
    ADD_SLOT(ParticleIdMap      , id_map          , INPUT, REQUIRED );
    ADD_SLOT(GridT              , grid            , INPUT_OUTPUT );
    ADD_SLOT(SLGrid             , sl_grid         , OUTPUT );

  public:

    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      size_t nc = sliplink_config->number_of_chains;
      size_t n_beads = sliplink_config->beads_per_chain;      
      size_t n_sl = sliplink_config->number_of_sliplinks;
      assert( n_sl%2 == 0 );

      const double sigma3 = sliplink_config->sigma3;

#     ifndef NDEBUG
      const double bond_max_search_dist = sliplink_config->bond_max_dist * ( 1. + *bond_max_stretch );
#     endif

      assert( id_map->size() == grid->number_of_particles() );
//      lout << "id_map size = " << id_map->size()<<", grid particles = "<< grid->number_of_particles() << std::endl;

      int nprocs = 1;
      int rank = 0;
      MPI_Comm_size(comm,&nprocs);
      MPI_Comm_rank(comm,&rank);
      
      // ========================================================
      // === 1) assign sliplinks handled by local mpi process ===
      // ========================================================
      size_t sl_start = sub_range_begin( rank, nprocs, 0, n_sl );
      size_t sl_end = sub_range_end( rank, nprocs, 0, n_sl );
      size_t sl_count = sl_end - sl_start;
      assert( sl_start == (rank*n_sl)/nprocs );
      assert( sl_end == ((rank+1)*n_sl)/nprocs );

      std::vector< SLBeadPlacement > locally_assigned_sliplinks;

      {
        // randomly choosen locations (bead id and bond placement)
        std::vector< SLBeadPlacement > sl_positions(sl_count);
        #pragma omp parallel
        {
          auto& re = onika::parallel::random_engine();
          std::uniform_int_distribution<size_t> random_chain( 0 , nc-1 );
          std::uniform_int_distribution<size_t> random_bead( 0 , n_beads-2 );
          std::uniform_real_distribution<double> uniform_01( 0.0 , 1.0 );
          for(size_t sl=0;sl<sl_count;sl++)
          {
            sl_positions[sl].sl_id = sl_start + sl;
            sl_positions[sl].bead_id = random_chain(re)*n_beads + random_bead(re);
            sl_positions[sl].frac = uniform_01(re);
          }
        }

        std::sort( sl_positions.begin(), sl_positions.end() , [](const SLBeadPlacement& a, const SLBeadPlacement& b) ->bool { return a.bead_id < b.bead_id; } );
        std::vector<int> send_bead_count( nprocs, 0 );

#       ifndef NDEBUG
        ssize_t oldp = 0;
#       endif
        for(size_t i=0;i<sl_count;i++)
        {
          size_t bead_id = sl_positions[i].bead_id;
          ssize_t chain_id = bead_id / n_beads;
          ssize_t p = sub_range_index(chain_id,nprocs,0,nc);
#         ifndef NDEBUG
          assert( chain_id>=0 && chain_id < static_cast<ssize_t>(nc) );
          assert( p>=0 && p<nprocs );
          assert( p >= oldp ); // strictly increasing owner, guaranteed by previous sort. 
          oldp=p;
          assert( chain_id>=sub_range_begin(p,nprocs,0,nc) && chain_id<sub_range_end(p,nprocs,0,nc) );
#         endif
          ++ send_bead_count[p];
        }
        std::vector<int> recv_bead_count( nprocs, -1 );
        MPI_Alltoall( send_bead_count.data() , 1 , MPI_INT , recv_bead_count.data() , 1 , MPI_INT , comm );
        assert( send_bead_count[rank] == recv_bead_count[rank] );
        
        std::vector<int> send_bead_displ( nprocs, 0 );
        std::vector<int> recv_bead_displ( nprocs, 0 );
        
        size_t total_received_beads = 0;
        size_t total_sent_beads = 0;
        for(int i=0;i<nprocs;i++)
        {
          send_bead_displ[i] = total_sent_beads;
          recv_bead_displ[i] = total_received_beads;
          total_sent_beads += send_bead_count[i];
          total_received_beads += recv_bead_count[i];
          // we'll work with raw bytes
          send_bead_count[i] *= sizeof(SLBeadPlacement);
          recv_bead_count[i] *= sizeof(SLBeadPlacement);
          send_bead_displ[i] *= sizeof(SLBeadPlacement);
          recv_bead_displ[i] *= sizeof(SLBeadPlacement);
        }
        assert( total_sent_beads == sl_count );
        locally_assigned_sliplinks.resize( total_received_beads );
        MPI_Alltoallv( sl_positions.data() , send_bead_count.data(), send_bead_displ.data(), MPI_BYTE, locally_assigned_sliplinks.data(), recv_bead_count.data(), recv_bead_displ.data(), MPI_BYTE, comm );
      }
  
      // ========================================================
      // === 2) initializes sliplinks in a specific grid      ===
      // ========================================================
  
      // at initialization, all MPI processes cover the whole grid, but they have a different subset of chains
      IJK dims = domain->grid_dimension();

      auto cells = grid->cells();

      // initializes support grid for sliplinks
      sl_grid->clear_particles();
      sl_grid->set_max_neighbor_distance( 0.0 );
      sl_grid->set_offset( IJK{0,0,0} );
      sl_grid->set_origin( domain->bounds().bmin );
      sl_grid->set_cell_size( domain->cell_size() );
      sl_grid->set_dimension( dims );
      auto sl_cells = sl_grid->cells();
  
      // random distributions for sliplinks' anchors
      auto& all_re = onika::parallel::random_engine();
      std::normal_distribution<double> gaussian_anchor_displ(0.0, sigma3);

      for( const auto & random_sl : locally_assigned_sliplinks )
      {
        double frac = random_sl.frac;
        uint64_t sl_id = random_sl.sl_id;
        uint64_t left_bead_id = random_sl.bead_id;

        //uint64_t peer_sl_id = sl_id ^ 1ull;
        uint64_t right_bead_id = left_bead_id + 1;
        
        assert( left_bead_id>=0 && left_bead_id<(nc*n_beads) );
        assert( (left_bead_id%n_beads)>=0 && (left_bead_id%n_beads)<(n_beads-1) );
        assert( (left_bead_id/n_beads)>=0 && (left_bead_id/n_beads)<nc );

        // where in the grid is the selected bead pair
        // note that at() will abort if element is not found.
        // both left and right beads are assumed to be in the local chains set
        uint64_t local_left_bead = global_to_own_local_id( left_bead_id, *id_map, *grid );
        assert( is_particle_id_valid(local_left_bead) );
        size_t lc=0,lp;
        decode_cell_particle( local_left_bead , lc, lp ); // ->at() virtually ok, only one bead can be found at this time. may not compile though
        assert( grid->is_valid_cell_particle(lc,lp) );

        uint64_t local_right_bead = global_to_own_local_id( right_bead_id, *id_map, *grid );
        assert( is_particle_id_valid(local_right_bead) );
        size_t rc=0,rp=0;
        decode_cell_particle( local_right_bead , rc, rp );
        assert( grid->is_valid_cell_particle(rc,rp) );
        assert( lc!=rc || lp!=rp );
        
        Vec3d left_r  = { cells[lc][field::rx][lp] , cells[lc][field::ry][lp] , cells[lc][field::rz][lp] };
        Vec3d right_r = { cells[rc][field::rx][rp] , cells[rc][field::ry][rp] , cells[rc][field::rz][rp] };
        Vec3d periodic_right_r = find_periodic_closest_point( right_r , left_r , domain->bounds() );

#       ifndef NDEBUG
        assert( norm(periodic_right_r-left_r) < bond_max_search_dist );
        Vec3d check_right_r = periodic_right_r;
        domain_periodic_location( *domain , check_right_r );
        assert( norm(check_right_r-right_r) < 1.e-30 );
#       endif        

        Vec3d sl_r = ( left_r * (1.0-frac) ) + ( periodic_right_r * frac );
        IJK loc = domain_periodic_location( *domain, sl_r );

        assert( grid->contains(loc) );
        size_t cell_index = grid_ijk_to_index( dims , loc );

        // sliplink anchor point follow a gaussian distribution around the sliplink position
        Vec3d anchor_displ = { gaussian_anchor_displ(all_re),
                               gaussian_anchor_displ(all_re),
                               gaussian_anchor_displ(all_re) };
    
        SLTupleInput t = { sl_r.x, sl_r.y, sl_r.z, sl_id, SlipLinkField{ anchor_displ, frac, left_bead_id } };
        sl_cells[cell_index].push_back(t);
      }

      sl_grid->rebuild_particle_offsets();
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return
R"EOF(
Written by Michelin & CEA/DIF
Credits to Claire Lemarchand, Ioanis Tanis, Thierry Carrard

Initializes the bead chains and sliplinks.
Number of chains is controlled by number_of_chains input.
Every chain has beads_per_chain beads.
)EOF";
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(sliplink_sl_init)
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_sl_init", make_grid_variant_operator< SlipLinkInitOperator > );
  }

}

