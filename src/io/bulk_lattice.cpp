#include <exanb/grid_cell_particles/bulk_generator.h>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT> using BulkLatticeTmpl = exanb::BulkLattice<GridT,field::_type>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory("bulk_lattice", make_grid_variant_operator< BulkLatticeTmpl >);
  }

}
