#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_random.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/particle_id_translation.h>
#include <onika/memory/allocator.h> // for ONIKA_ASSUME_ALIGNED macro
#include <exanb/core/thread.h> // for ONIKA_ASSUME_ALIGNED macro
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/sliplink/sliplink.h>

#include <random>
#include <mpi.h>

#include <vector>

namespace exaStamp
{

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az >
    >
  class SlipLinkForceOperator : public OperatorNode
  {
    ADD_SLOT(MPI_Comm           , mpi             , INPUT , MPI_COMM_WORLD );
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT, REQUIRED );
    ADD_SLOT(double             , bond_max_stretch , INPUT , 0.5 ); // fraction of bond_max_dist.
    ADD_SLOT(GridParticleLocks          , particle_locks  , INPUT, REQUIRED );      
    ADD_SLOT(SLGrid             , sl_grid         , INPUT, REQUIRED );
    ADD_SLOT(ParticleIdMap      , id_map          , INPUT, REQUIRED );
    ADD_SLOT(GridT              , grid            , INPUT_OUTPUT );

  public:

    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      int nprocs = 1;
      int rank = 0;
      MPI_Comm_size(comm,&nprocs);
      MPI_Comm_rank(comm,&rank);

      GridParticleLocks& particle_locks = *(this->particle_locks);

      auto cells = grid->cells();
      auto sl_cells = sl_grid->cells();

      assert( grid->number_of_cells() == sl_grid->number_of_cells() );  
          
      IJK dims = grid->dimension();
      assert( dims == sl_grid->dimension() );
//      size_t ghost_layers = grid->ghost_layers();
 
//      const size_t n_chains = sliplink_config->number_of_chains;
//      const size_t n_beads = sliplink_config->beads_per_chain;
      const double cte8_cte2 = sliplink_config->cte8 * sliplink_config->cte2;
      const double sigma4 = sliplink_config->sigma4;
      const double bond_max_search_dist = sliplink_config->bond_max_dist * ( 1. + *bond_max_stretch );
     
#     pragma omp parallel
      {
        std::normal_distribution<double> gaussian_friction(0.0, sigma4);
        GRID_OMP_FOR_BEGIN(dims,i,_)
        {
          SlipLinkField * __restrict__ sl = sl_cells[i][field::sl]; ONIKA_ASSUME_ALIGNED(sl);
          double * __restrict__ sl_rx = sl_cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(sl_rx);
          double * __restrict__ sl_ry = sl_cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(sl_ry);
          double * __restrict__ sl_rz = sl_cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(sl_rz);
          size_t n = sl_cells[i].size();

          for(size_t j=0;j<n;j++)
          {
            Vec3d sl_force = sl[j].anchor_displ;
            Vec3d sl_pos = { sl_rx[j], sl_ry[j], sl_rz[j] };
            double xj_f = sl[j].xj_frac;
            assert( xj_f>=0.0 && xj_f<=1.0 );
            uint64_t left_bead_id = sl[j].left_bead_id;
            uint64_t right_bead_id = left_bead_id+1;
            assert( is_particle_id_valid(left_bead_id) && is_particle_id_valid(right_bead_id) );

            // left bead
            uint64_t local_left_bead_id = global_to_nearest_local_id( left_bead_id, *id_map, *grid, sl_pos, bond_max_search_dist );
            if( is_particle_id_valid(local_left_bead_id) )
            {
              size_t lc = 0;
              size_t lp = 0;
              decode_cell_particle( local_left_bead_id , lc , lp );
              assert( grid->is_valid_cell_particle(lc,lp) );
              if( ! grid->is_ghost_cell(lc) )
              {
                particle_locks[lc][lp].lock();
                cells[lc][field::ax][lp] += sl_force.x * (1.0-xj_f) * cte8_cte2;
                cells[lc][field::ay][lp] += sl_force.y * (1.0-xj_f) * cte8_cte2;
                cells[lc][field::az][lp] += sl_force.z * (1.0-xj_f) * cte8_cte2;
                particle_locks[lc][lp].unlock();
              }
            }
            
            // right bead
            uint64_t local_right_bead_id = global_to_nearest_local_id( right_bead_id, *id_map, *grid, sl_pos, bond_max_search_dist );
            if( is_particle_id_valid(local_right_bead_id) )
            {
              size_t rc = 0;
              size_t rp = 0;
              decode_cell_particle( local_right_bead_id , rc , rp );
              assert( grid->is_valid_cell_particle(rc,rp) );
              if( ! grid->is_ghost_cell(rc) )
              {
                particle_locks[rc][rp].lock();
                cells[rc][field::ax][rp] += sl_force.x * xj_f * cte8_cte2;
                cells[rc][field::ay][rp] += sl_force.y * xj_f * cte8_cte2;
                cells[rc][field::az][rp] += sl_force.z * xj_f * cte8_cte2;
                particle_locks[rc][rp].unlock();
              }
            }
            
          }
        }
        GRID_OMP_FOR_END
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

computes stretch force between successive beads in a chain
)EOF";
    }

  };

  template<class GridT> using SlipLinkForceOperatorTmpl = SlipLinkForceOperator<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_sl_force", make_grid_variant_operator< SlipLinkForceOperatorTmpl > );
  }

}

