#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/parallel/random.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/particle_id_codec.h>
#include <onika/memory/allocator.h> // for ONIKA_ASSUME_ALIGNED macro
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/sliplink/sliplink.h>
#include <exanb/core/particle_id_translation.h>

#include <random>
#include <mpi.h>

#include <vector>

namespace exaStamp
{
  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id >
    >
  class SlipLinkCheckBonds : public OperatorNode
  {  
    ADD_SLOT(SlipLinkParameters , sliplink_config   , INPUT , REQUIRED );
    ADD_SLOT(double             , bond_max_stretch  , INPUT , 0.5 ); // fraction of bond_max_dist.
    ADD_SLOT(ParticleIdMap      , id_map            , INPUT , REQUIRED );
    ADD_SLOT(GridT              , grid              , INPUT , REQUIRED );

  public:

    inline void execute () override final
    {
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);

      const double bond_max_search_dist = sliplink_config->bond_max_dist * ( 1. + *bond_max_stretch );
      const double bond_max_dist2 = bond_max_search_dist * bond_max_search_dist;
      const double bond_min_dist2 = bond_max_dist2*1.e-18;

      // ldbg << "bond_max_search_dist = "<<bond_max_search_dist<<std::endl;
  
      // this simple test ensures that we're connected with the right id_map
      assert( grid->number_of_particles() == id_map->size() );

      double max_d2 = std::numeric_limits<double>::lowest();
      double min_d2 = std::numeric_limits<double>::max();

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc, reduction(min:min_d2) reduction(max:max_d2) )
        {
          size_t i = grid_ijk_to_index( dims , loc + ghost_layers);
          assert( ! grid->is_ghost_cell(i) );

          uint64_t const * __restrict__ ids = cells[i][field::id]; ONIKA_ASSUME_ALIGNED(ids);
          double const * __restrict__ rx = cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
          double const * __restrict__ ry = cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
          double const * __restrict__ rz = cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
          size_t n = cells[i].size();

          for(size_t j=0;j<n;j++)
          {
            Vec3d r = { rx[j], ry[j], rz[j] };

            // 0=left bead, 1=right bead
            uint64_t id = ids[j];
            uint64_t siblings[2] = { sliplink_left_bead_id(id,sliplink_config->beads_per_chain) , sliplink_right_bead_id(id,sliplink_config->beads_per_chain) };
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
                    lerr << "bead #"<<id<<" sibling("<<k<<")="<<bead_id<<" not found in id_map. r=" << r <<std::endl;
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
                    lerr << "bead #"<<id<<" sibling("<<k<<") not found in grid" << std::endl; 
                    std::abort();                  
                  }
                }
                Vec3d sr = { cells[c][field::rx][p] , cells[c][field::ry][p] , cells[c][field::rz][p] };
                double d2 = norm2(r-sr);

                min_d2 = std::min( min_d2 , d2 );
                max_d2 = std::max( max_d2 , d2 );
                
                if( d2 > bond_max_dist2 || d2<bond_min_dist2)
                {
#                 pragma omp critical
                  {
                    lerr << "bad bond distance "<< std::sqrt(d2) << " > "<< sliplink_config->bond_max_dist << std::endl;
                    std::abort();
                  }
                }
                if( cells[c][field::id][p] != bead_id )
                {
#                 pragma omp critical
                  {
                    lerr << "bead #"<<id<<", sibling("<<k<<"), inconsistent id" <<cells[c][field::id][p] << std::endl;
                    std::abort();                  
                  }
                }
              }
            }
          }
        }
        GRID_OMP_FOR_END
      }
      
      ldbg << "bonds dist range = [ "<<std::sqrt(min_d2)<<" ; "<<std::sqrt(max_d2)<<" ]"<<std::endl;
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

  template<class GridT> using SlipLinkCheckBondsTmpl = SlipLinkCheckBonds<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(sliplink_check_bonds)
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_check_bonds", make_grid_variant_operator< SlipLinkCheckBondsTmpl > );
  }

}

