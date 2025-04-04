#include <exanb/grid_cell_particles/lattice_generator.h>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT> using RegionLatticeTmpl = exanb::RegionLattice<GridT,field::_type>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(region_lattice)
  {
    OperatorNodeFactory::instance()->register_factory("lattice", make_grid_variant_operator< RegionLatticeTmpl >);
  }

}
