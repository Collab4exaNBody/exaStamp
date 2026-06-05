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

#include "igar_force_from_gradient.h"

namespace exaStamp
{

  using namespace exanb;

  template<class GridT, class = AssertGridHasFields<GridT,field::_rx,field::_ry,field::_rz,field::_fx,field::_fy,field::_fz> >
  class IgarForceFromGradient : public OperatorNode
  {
    ADD_SLOT( GridT          , grid             , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( GridCellValues , grid_cell_values , INPUT        , REQUIRED );
    ADD_SLOT( double         , energy_factor    , INPUT        , 1.0, DocString{"Prefactor applied to -igar_dEdx, -igar_dEdy and -igar_dEdz before adding them to particle forces."} );

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(
Applies precomputed IGAR gradient fields to particles, accumulating forces
and potential energy.

Requires igar_compute_gradient to have been called beforehand to populate
igar_dEdx, igar_dEdy and igar_dEdz in grid_cell_values.

For each particle, the operator locates its subcell and applies:
  f += -energy_factor * (igar_dEdx, igar_dEdy, igar_dEdz)
using nearest-subcell assignment for forces. The potential energy ep is
accumulated using trilinear interpolation of igar_ep between the 8
surrounding subcell centers.

Pair with igar_compute_gradient placed in setup_system for best performance
when igar_ep is a static field.

Example (LJ_igar.msp — gradient precomputed once at setup):

  setup_system:
    - read_cell_values:
        field_name: "igar_ep"
        field_dim: 1
        grid_subdiv: 64
        file_name: "igar_test.vtk"
    - igar_compute_gradient

  compute_force:
    - lj_compute_force
    - igar_force_from_gradient: { energy_factor: 20000.0 }
)EOF";
    }

    inline void execute() override final
    {
      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args ...); };
      ParticleCellProjectionTools::get_particle_force_from_gradient_grid(
          ldbg
        , *grid
        , pecfunc
        , *grid_cell_values
        , *energy_factor );
    }
  };

  template<class GridT> using IgarForceFromGradientTmpl = IgarForceFromGradient<GridT>;

  ONIKA_AUTORUN_INIT(igar_force_from_gradient)
  {
    OperatorNodeFactory::instance()->register_factory( "igar_force_from_gradient", make_grid_variant_operator<IgarForceFromGradientTmpl> );
  }

}
