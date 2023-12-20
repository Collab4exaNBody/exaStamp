#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/make_grid_variant_operator.h>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ResizeGridCellValues : public OperatorNode
  {
    ADD_SLOT( GridT       , grid            , INPUT , REQUIRED );
    ADD_SLOT( GridCellValues , grid_cell_values , INPUT_OUTPUT );

  public:
    inline void execute () override final
    {
      grid_cell_values->set_grid_dims( grid->dimension() );
      grid_cell_values->set_ghost_layers( grid->ghost_layers() );
      grid_cell_values->set_grid_offset( grid->offset() );
    }

  };

  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "resize_grid_cell_values", make_grid_variant_operator<ResizeGridCellValues> );
  }

}

