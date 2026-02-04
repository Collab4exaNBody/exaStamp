/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <exaStamp/molecule/extract_connectivity.h>
#include <exaStamp/molecule/get_r.h>
#include <exaStamp/molecule/id_map.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_id, field::_type >
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
  ONIKA_AUTORUN_INIT(mol_extract_connectivity)
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "mol_extract_connectivity", make_grid_variant_operator< MolExtractConnectivityTmpl > );
  }

}
