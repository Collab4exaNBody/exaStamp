#include <exanb/mpi/simple_cost_model.h>
#include <exanb/core/make_grid_variant_operator.h>

namespace exanb
{

  template<class GridT> using SimpleCostMultimat = SimpleCostModel<GridT,field::_type>;

  // === register factory ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory(
      "simple_cost_multimat",
      make_grid_variant_operator< SimpleCostMultimat > );
  }

}

