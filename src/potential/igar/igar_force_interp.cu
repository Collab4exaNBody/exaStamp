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
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <exaStamp/compute/field_combiners.h>
#include <exaStamp/compute/physics_functors.h>
#include <exanb/compute/field_combiners.h>
#include <exanb/core/grid_particle_field_accessor.h>

#include "igar_force_interp.h"

namespace exaStamp
{

  using namespace exanb;

  template <class GridT, class = AssertGridHasFields<GridT, field::_rx, field::_ry, field::_rz, field::_ep, field::_fx, field::_fy, field::_fz>> class IgarForceInterp : public OperatorNode
  {
    ADD_SLOT(GridT, grid, INPUT_OUTPUT);
    ADD_SLOT(GridCellValues, grid_cell_values, INPUT);
    ADD_SLOT(double, energy_factor, INPUT, REQUIRED, DocString{"Energy scale of the igar_ep field (an energy: internal units unless a unit is given, e.g. 2.0 eV). Energy = energy_factor * igar_ep, force = -energy_factor * grad(igar_ep)."});
    ADD_SLOT(std::string, interpolation, INPUT, std::string("catmull_rom"), DocString{"Interpolation of the igar_ep subcell values: catmull_rom (tricubic, 4x4x4 subcells, interpolating, energy and force continuous), bspline (quadratic B-spline, 3x3x3 subcells, continuous, smooths the field) or lagrange (triquadratic, 3x3x3 subcells, energy discontinuous at subcell faces: no energy conservation)."});

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Computes per-particle potential energy and forces from the IGAR scalar energy
field igar_ep stored in grid_cell_values (one value per subcell, located at
the subcell center).

The energy at the particle position is interpolated from the neighboring
subcell values and multiplied by energy_factor; the force is minus its
analytical gradient, so forces are the exact derivative of the energy.

interpolation selects the scheme:
  catmull_rom (default)  tricubic Catmull-Rom, 4x4x4 subcells. Interpolating,
                         energy and force continuous: energy is conserved.
  bspline                quadratic B-spline, 3x3x3 subcells. Energy and force
                         continuous, cheaper, but smooths the field (energy
                         error ~ subcell_size^2 * curvature / 8).
  lagrange               triquadratic Lagrange, 3x3x3 subcells centered on the
                         particle's subcell. The energy jumps when a particle
                         crosses a subcell face: energy is not conserved.

The stencil reads igar_ep up to 2 subcells (catmull_rom) or 1 subcell away
from the particle's subcell, in neighboring cells and ghost cells included:
the ghost layer must cover that reach (values outside the grid count as 0).

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
      auto pecfunc = [self = this](auto... args) { return self->parallel_execution_context(args...); };

      IgarInterpolation scheme = IgarInterpolation::CatmullRom;
      if (*interpolation == "lagrange") scheme = IgarInterpolation::Lagrange;
      else if (*interpolation == "bspline") scheme = IgarInterpolation::BSpline;
      else if (*interpolation == "catmull_rom") scheme = IgarInterpolation::CatmullRom;
      else fatal_error() << "igar_force_interp: unknown interpolation '" << *interpolation << "' (lagrange, bspline or catmull_rom)" << std::endl;

      ParticleCellProjectionTools::get_particle_force_from_grid(ldbg, *grid, pecfunc, *grid_cell_values, *energy_factor, scheme);
    }
  };

  template <class GridT> using IgarForceInterpTmpl = IgarForceInterp<GridT>;

  ONIKA_AUTORUN_INIT(igar_force_interp) { OperatorNodeFactory::instance()->register_factory("igar_force_interp", make_grid_variant_operator<IgarForceInterpTmpl>); }

} // namespace exaStamp
