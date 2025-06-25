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
    , class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_id >
    >
  class SlipLinkBondsStretchOperator : public OperatorNode
  {  
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT, REQUIRED );
    ADD_SLOT(double             , bond_max_stretch , INPUT , 0.5 ); // fraction of bond_max_dist.
    ADD_SLOT(ParticleIdMap      , id_map          , INPUT );
    ADD_SLOT(GridT              , grid            , INPUT_OUTPUT );

  public:

    inline void execute () override final
    {
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);

      const double bond_max_search_dist = sliplink_config->bond_max_dist * ( 1. + *bond_max_stretch );
      const double cte7_cte2 = sliplink_config->cte7 * sliplink_config->cte2;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghost)
        {
          IJK loc = loc_no_ghost + ghost_layers;
          size_t i = grid_ijk_to_index(dims,loc);

          uint64_t const * __restrict__ ids = cells[i][field::id]; ONIKA_ASSUME_ALIGNED(ids);
          double const * __restrict__ rx = cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
          double const * __restrict__ ry = cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
          double const * __restrict__ rz = cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
          double * __restrict__ fx = cells[i][field::fx]; ONIKA_ASSUME_ALIGNED(fx);
          double * __restrict__ fy = cells[i][field::fy]; ONIKA_ASSUME_ALIGNED(fy);
          double * __restrict__ fz = cells[i][field::fz]; ONIKA_ASSUME_ALIGNED(fz);
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
                assert( is_particle_id_valid(local_bead) );
                size_t c = 0, p = 0;
                decode_cell_particle( local_bead , c , p );
                assert( grid->is_valid_cell_particle(c,p) );

                // Ajout de la force de ressort
                fx[j] += ( cells[c][field::rx][p] - r.x ) * cte7_cte2 ;
                fy[j] += ( cells[c][field::ry][p] - r.y ) * cte7_cte2 ;
                fz[j] += ( cells[c][field::rz][p] - r.z ) * cte7_cte2 ;
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

  template<class GridT> using SlipLinkBondsStretch = SlipLinkBondsStretchOperator<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(sliplink_bond_stretch)
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_bond_stretch", make_grid_variant_operator< SlipLinkBondsStretch > );
  }

}

