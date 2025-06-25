#include <md/snap/snap_force.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <onika/cpp_utils.h>

namespace exaStamp
{
  template<class GridT> using SnapForceXSTmpl = md::SnapNewForce<GridT,field::_ep,field::_virial>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(snap_force)
  {
    OperatorNodeFactory::instance()->register_factory( "snap_force" ,make_grid_variant_operator< SnapForceXSTmpl > );
  }

}


