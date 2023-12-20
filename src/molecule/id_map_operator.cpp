//#include <chrono>
#include <memory>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/molecule/id_map.h>

namespace exaStamp
{

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id, field::_type >
    >
  struct IdMapNode : public OperatorNode
  {
    ADD_SLOT( GridT , grid   , INPUT);
    ADD_SLOT( IdMap , id_map , INPUT_OUTPUT);
    ADD_SLOT( IdMapGhosts , id_map_ghosts , INPUT_OUTPUT);

    inline void execute ()  override final
    {
      // clear previous content
      id_map->clear();
      id_map_ghosts->clear();

      auto cells = grid->cells();
      size_t n_cells = grid->number_of_cells();

#     pragma omp parallel for schedule(dynamic)
      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        const uint64_t* __restrict__ ids   = cells[cell_i][field::id];
        const uint8_t*  __restrict__ types = cells[cell_i][field::type];
        size_t n = cells[cell_i].size();
        if(!grid->is_ghost_cell(cell_i))
        {
          for(size_t p_i=0;p_i<n;p_i++)
          {
            id_map->insert( { ids[p_i], encode_cell_particle(cell_i, p_i, types[p_i]) } ); // IdMap is thread safe for concurrent insertion, see exanb/molecule/id_map.h
          }
        }
        else
        {
          for(size_t p_i=0;p_i<n;p_i++)
          {
            id_map_ghosts->insert( { ids[p_i], encode_cell_particle(cell_i, p_i, types[p_i]) } ); // IdMapGhosts is thread safe for concurrent insertion, see exanb/molecule/id_map.h
          }
        }
      }
    }

  };

  template<class GridT> using IdMapNodeTmpl = IdMapNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "id_map", make_grid_variant_operator<IdMapNodeTmpl> );
  }

}
