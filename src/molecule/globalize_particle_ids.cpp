#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_random.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/parallel_grid_algorithm.h>

#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/particle_id_constants.h>

#include <random>
#include <mpi.h>

#include <vector>
#include <utility>

namespace exaStamp
{

  template<
      class GridT
    , class = AssertGridHasFields< GridT, field::_id, field::_cmol >
    >
  class GlobalizeParticleIdsOperator : public OperatorNode
  {
    ADD_SLOT(GridT  , grid    , INPUT_OUTPUT );

  public:

    inline void execute () override final
    {
      auto cells = grid->cells();
               
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);
      
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghost)
        {
          IJK loc = loc_no_ghost + ghost_layers;
          size_t i = grid_ijk_to_index(dims,loc);
          MoleculeConnectivity * __restrict__ cmol = cells[i][field::cmol];
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            for(size_t k=0;k<cmol[j].size();k++)
            {
              // can be PARTICLE_NO_ID, but bonded particle must be known as we traverse only non ghost particles
              assert( cmol[j][k] != PARTICLE_MISSING_ID ); 
              if( is_particle_id_valid( cmol[j][k] ) )
              {
                size_t cell_index = 0;
                size_t part_index = 0;
                decode_cell_particle( cmol[j][k] , cell_index , part_index );
		assert( grid->is_valid_cell_particle(cell_index,part_index) );
                cmol[j][k] = cells[cell_index][field::id][part_index];
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
replaces local ids in cmol by real unique ids of particles. Does it only for inner particles (not ghosts)
)EOF";
    }

  };

  template<class GridT> using GlobalizeParticleIds = GlobalizeParticleIdsOperator<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "globalize_particle_ids", make_grid_variant_operator< GlobalizeParticleIds > );
  }

}

