#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_random.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/particle_id_codec.h>
#include <onika/memory/allocator.h> // for ONIKA_ASSUME_ALIGNED macro
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/sliplink/sliplink.h>
#include <exanb/core/particle_id_translation.h>

#include <random>
#include <mpi.h>

#include <vector>
#include <iomanip>

namespace exaStamp
{
  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id >
    >
  class SlipLinkCheckSL : public OperatorNode
  {  
    ADD_SLOT(SlipLinkParameters , sliplink_config   , INPUT , REQUIRED );
    ADD_SLOT(double             , bond_max_stretch , INPUT , 0.5 ); // fraction of bond_max_dist.
    ADD_SLOT(ParticleIdMap      , id_map            , INPUT , REQUIRED );
    ADD_SLOT(GridT              , grid              , INPUT , REQUIRED );
    ADD_SLOT(SLGrid             , sl_grid           , INPUT , REQUIRED );

  public:

    inline void execute () override final
    {
      IJK dims = sl_grid->dimension();
      size_t ghost_layers = sl_grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);
      
      assert( sl_grid->dimension() == grid->dimension() );
      assert( sl_grid->ghost_layers() == grid->ghost_layers() );

      const double bond_max_search_dist = sliplink_config->bond_max_dist * ( 1. + *bond_max_stretch );
      const double bond_max_dist2 = bond_max_search_dist * bond_max_search_dist;

      auto cells = grid->cells();
      auto sl_cells = sl_grid->cells();
      
      IJK grid_offset = grid->offset();
      double cell_size = grid->cell_size();
      Vec3d origin = grid->origin();
      
      AABB inner_bounds = { origin + (grid_offset+ghost_layers)*cell_size , origin + (grid_offset+dims-2*ghost_layers)*cell_size };
      AABB outter_bounds = { inner_bounds.bmin - grid->max_neighbor_distance() , inner_bounds.bmax + grid->max_neighbor_distance() };

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc)
        {
          size_t i = grid_ijk_to_index( dims , loc + ghost_layers);
          assert( ! sl_grid->is_ghost_cell(i) );

          uint64_t const * __restrict__ ids = sl_cells[i][field::id]; ONIKA_ASSUME_ALIGNED(ids);
          double const * __restrict__ rx = sl_cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
          double const * __restrict__ ry = sl_cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
          double const * __restrict__ rz = sl_cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
          SlipLinkField const * __restrict__ sl = sl_cells[i][field::sl]; ONIKA_ASSUME_ALIGNED(sl);
          size_t n = sl_cells[i].size();

          for(size_t j=0;j<n;j++)
          {
            Vec3d r = { rx[j], ry[j], rz[j] };
            uint64_t sl_id = ids[j];

            // 0=left bead, 1=right bead
#           ifndef NDEBUG
            if( sliplink_right_bead_id(sl[j].left_bead_id,sliplink_config->beads_per_chain) != (sl[j].left_bead_id+1) )
            {
              lerr << "SL #"<<sl_id<<" left bead ("<<sl[j].left_bead_id<<") inconsistent with right bead" <<std::endl;
              std::abort();
            }
#           endif
            uint64_t siblings[2] = { sl[j].left_bead_id , sl[j].left_bead_id+1 };
            
            for(size_t k=0;k<2;k++)
            {
              uint64_t bead_id = siblings[k];
              if( is_particle_id_valid(bead_id) )
              {
                uint64_t local_bead = global_to_nearest_local_id( bead_id, *id_map, *grid, r, bond_max_search_dist );
                if( ! is_particle_id_valid(local_bead) )
                {
#                 pragma omp critical
                  {
                    lerr << "SL #"<<sl_id<<" bead("<<k<<")="<<bead_id<<" not found in id_map. r = "
                      <<std::setprecision(16)<<r.x<<" , " <<std::setprecision(16)<<r.y<<" , " <<std::setprecision(16)<<r.z <<std::endl;
                    auto range = id_map->equal_range(bead_id);
                    for (auto it=range.first; it!=range.second; ++it)
                    {
                      size_t c=0,p=0;
                      decode_cell_particle( it->second , c, p );
                      lerr << "\tcell #"<<c<<" part #"<<p<<" (id="<<cells[c][field::id][p]<<")";
                      if( grid->is_ghost_cell(c) ) { lerr << " (ghost)"; }
                      Vec3d sr = { cells[c][field::rx][p] , cells[c][field::ry][p] , cells[c][field::rz][p] };
                      lerr << " pos="<<sr << ", dist="<<norm(sr-r)<< std::endl;
                    }
                    lerr << "inner_bounds = "<<inner_bounds<<std::endl;
                    lerr << "outter_bounds = "<<outter_bounds<<std::endl;
                    lerr << std::flush;
                    std::abort();
                  }
                }
                size_t c = 0, p = 0;
                decode_cell_particle( local_bead , c , p );
                if( ! grid->is_valid_cell_particle(c,p) )
                {
#                 pragma omp critical
                  {
                    lerr << "SL #"<<sl_id<<" sibling("<<k<<") not found in grid" << std::endl;
                    std::abort();                  
                  }
                }
                Vec3d sr = { cells[c][field::rx][p] , cells[c][field::ry][p] , cells[c][field::rz][p] };
                double d2 = norm2(r-sr);
                if( d2 > bond_max_dist2 )
                {
#                 pragma omp critical
                  {
                    lerr << "bad distance to SL "<< std::sqrt(d2) << " > "<< sliplink_config->bond_max_dist << std::endl;
                    std::abort();
                  }
                }
                if( cells[c][field::id][p] != bead_id )
                {
#                 pragma omp critical
                  {
                    lerr << "SL #"<<sl_id<<", sibling("<<k<<"), inconsistent id" <<cells[c][field::id][p] << std::endl;
                    std::abort();                  
                  }
                }
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

  template<class GridT> using SlipLinkCheckSLTmpl = SlipLinkCheckSL<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_check_sl", make_grid_variant_operator< SlipLinkCheckSLTmpl > );
  }

}

