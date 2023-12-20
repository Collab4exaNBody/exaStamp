#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/log.h>
#include <exaStamp/molecule/extract_connectivity.h>
#include <exaStamp/molecule/get_r.h>
#include <exaStamp/molecule/id_map.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_id, field::_idmol, field::_cmol>
    >
  class MolExtractConnectivity : public OperatorNode
  {
    ADD_SLOT( Domain            , domain             , INPUT );
    ADD_SLOT( GridT             , grid               , INPUT );
    ADD_SLOT( IdMap             , id_map             , INPUT  );
    ADD_SLOT( IdMapGhosts       , id_map_ghosts      , INPUT  );
    ADD_SLOT( ChemicalBonds     , chemical_bonds     , INPUT_OUTPUT );
    ADD_SLOT( ChemicalAngles    , chemical_angles    , INPUT_OUTPUT );
    ADD_SLOT( ChemicalTorsions  , chemical_torsions  , INPUT_OUTPUT );
    ADD_SLOT( ChemicalImpropers , chemical_impropers , INPUT_OUTPUT );

//    REGISTER_OPERATOR(set_connectivity);
  public:
    inline void execute () override final
    {
      
      extract_connectivity( domain->grid_dimension(), *grid, *chemical_bonds, *chemical_angles, *chemical_torsions, *chemical_impropers, *id_map, *id_map_ghosts );

    }
  };


  template<class GridT> using MolExtractConnectivityTmpl = MolExtractConnectivity<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "mol_extract_connectivity", make_grid_variant_operator< MolExtractConnectivityTmpl > );
  }

}
