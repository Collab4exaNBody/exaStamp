/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>

#include "force_from_grid.h"

namespace exanb
{

  template<class GridT, class = AssertGridHasFields<GridT,field::_rx,field::_ry,field::_rz,field::_fx,field::_fy,field::_fz> >
  class ForceFromGradientGrid : public OperatorNode
  {
    ADD_SLOT( GridT          , grid             , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( GridCellValues , grid_cell_values , INPUT        , REQUIRED );
    ADD_SLOT( double         , force_factor     , INPUT        , 1.0, DocString{"Prefactor applied to -igar_dEdx, -igar_dEdy and -igar_dEdz before adding them to particle forces."} );

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(
Adds a force to particles from the grid-cell-value gradient fields igar_dEdx,
igar_dEdy and igar_dEdz. For each particle, the operator locates the sub-cell
containing the particle and adds:

  f += -force_factor * (igar_dEdx, igar_dEdy, igar_dEdz)

The default force_factor is 1, corresponding to f = -grad(E).
)EOF";
    }

    inline void execute() override final
    {
      ParticleCellProjectionTools::get_particle_force_from_gradient_grid(
          ldbg
        , *grid
        , *grid_cell_values
        , *force_factor );
    }
  };

  template<class GridT> using ForceFromGradientGridTmpl = ForceFromGradientGrid<GridT>;

  ONIKA_AUTORUN_INIT(force_from_gradient_grid)
  {
    OperatorNodeFactory::instance()->register_factory( "force_from_gradient_grid", make_grid_variant_operator<ForceFromGradientGridTmpl> );
  }

}
