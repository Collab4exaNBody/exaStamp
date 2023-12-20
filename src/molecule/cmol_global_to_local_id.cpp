#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>

#include <exaStamp/molecule/mol_connectivity.h>

#include <memory>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_cmol >
    >
  class CMolGlobalToLocalIdOperator : public OperatorNode
  {  
    using IdMap = std::map<uint64_t, uint64_t>;

    ADD_SLOT( GridT , grid   , INPUT_OUTPUT);
    ADD_SLOT( IdMap , id_map , INPUT);

  public:
    inline void execute ()  override final
    {
      GridT& grid = *(this->grid);
      const IdMap& id_map = *(this->id_map);

      auto cells = grid.cells();

      size_t n_cells = grid.number_of_cells();
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) // FIXME: WRONG, works only if limited to own cells (without ghosts)
      {
        MoleculeConnectivity* __restrict__ cmol = cells[cell_i][field::cmol]; 
        size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++)
        {
          cmol[p_i][0] = id_map.at( cmol[p_i][0] );
          cmol[p_i][1] = id_map.at( cmol[p_i][1] );
          cmol[p_i][2] = id_map.at( cmol[p_i][2] );
          cmol[p_i][3] = id_map.at( cmol[p_i][3] );
        }
      }
    }

  };

  template<class GridT> using CMolGlobalToLocalIdOperatorTmpl = CMolGlobalToLocalIdOperator<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("cmol_global_to_local", make_grid_variant_operator< CMolGlobalToLocalIdOperatorTmpl > );
  }

}
