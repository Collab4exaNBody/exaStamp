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
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>

#include <exaStamp/compute/physics_functors.h>
#include <exanb/core/grid_particle_field_accessor.h>
#include <exanb/compute/field_combiners.h>
#include <exaStamp/compute/field_combiners.h>

#include "igar_force_interp.h"

namespace exaStamp
{

  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class IgarForceInterp : public OperatorNode
  {
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( GridCellValues            , grid_cell_values    , INPUT );
    ADD_SLOT( double                    , energy_factor       , INPUT        , REQUIRED );

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(
Computes per-particle potential energy and forces from the IGAR scalar energy
field igar_ep stored in grid_cell_values, using triquadratic Lagrange
interpolation.

For each particle, a 3x3x3 stencil of neighboring subcell energy values is
read and scaled by energy_factor. Triquadratic Lagrange basis functions
(nodes at -1, 0, +1 in normalized subcell coordinates) are used to evaluate
both the interpolated energy and its analytical gradient at the particle's
exact position. Forces are accumulated as f += -grad(E).

This operator is self-contained: no gradient precomputation is required.
Energy and forces are C1 continuous across subcell and cell boundaries.

Example (Lennard-Jones + IGAR, LJ_igar.msp):

  setup_system:
    - read_cell_values:
        field_name: "igar_ep"
        field_dim: 1
        grid_subdiv: 64
        file_name: "igar_test.vtk"

  compute_force:
    - lj_compute_force
    - igar_force_interp: { energy_factor: 20000.0 }

Example (REBO + IGAR, CH_rebo_read_xyz.msp):

  compute_force:
    - rebo_force
    - igar_force_interp: { energy_factor: 50000.0 }
)EOF";
    }
    inline void execute() override final
    {
      using namespace ParticleCellProjectionTools;
      if( grid->number_of_cells() == 0 ) return;
      ldbg << "Igar force from grid_cell_values" << std::endl;
      get_particle_force_from_grid( ldbg, *grid, *grid_cell_values, *energy_factor );
    }

  };

  template<class GridT> using IgarForceInterpTmpl = IgarForceInterp<GridT>;

  ONIKA_AUTORUN_INIT(igar_force_interp)
  {
    OperatorNodeFactory::instance()->register_factory("igar_force_interp", make_grid_variant_operator<IgarForceInterpTmpl>);
  }

}
