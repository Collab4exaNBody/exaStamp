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
#include <exanb/core/domain.h>
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
    ADD_SLOT(Domain, domain, INPUT, REQUIRED);
    ADD_SLOT(IgarField, igar_field, INPUT, REQUIRED, DocString{"IGAR field, read by igar_read_field."});
    ADD_SLOT(double, energy_factor, INPUT, REQUIRED, DocString{"Energy scale of the IGAR field (an energy: internal units unless a unit is given, e.g. 2.0 eV). Energy = energy_factor * field, force = -energy_factor * grad(field)."});
    ADD_SLOT(std::string, interpolation, INPUT, std::string("catmull_rom"), DocString{"Interpolation of the field values: catmull_rom (tricubic, 4x4x4 voxels, interpolating, energy and force continuous), bspline (quadratic B-spline, 3x3x3 voxels, continuous, smooths the field) or lagrange (triquadratic, 3x3x3 voxels, energy discontinuous at voxel faces: no energy conservation)."});

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Computes per-particle potential energy and forces from the IGAR scalar field
read by igar_read_field (one value per voxel, located at the voxel center, the
voxels splitting the simulation domain evenly).

The field value at the particle position is interpolated from the neighboring
voxel values and multiplied by energy_factor; the force is minus its
analytical gradient, so forces are the exact derivative of the energy.

interpolation selects the scheme:
  catmull_rom (default)  tricubic Catmull-Rom, 4x4x4 voxels. Interpolating,
                         energy and force continuous: energy is conserved.
  bspline                quadratic B-spline, 3x3x3 voxels. Energy and force
                         continuous, cheaper, but smooths the field (energy
                         error ~ voxel_size^2 * curvature / 8).
  lagrange               triquadratic Lagrange, 3x3x3 voxels centered on the
                         particle's voxel. The energy jumps when a particle
                         crosses a voxel face: energy is not conserved.

The field is attached to the simulation box: its voxels split the domain in
grid coordinates, so a deformed box (xform) deforms the field with it. Forces
are real-space derivatives of the energy (gradient times transpose(inv_xform)).

The stencil wraps around periodic boundaries; along a non-periodic axis, voxels
outside the domain count as 0. The field is held by every MPI process: no
ghost exchange of the field is needed, whatever the particle update scheme.

Example (REBO + IGAR, 8-bit image):

  setup_system:
    - igar_read_field: { file_name: "image_1024.vtk", storage: uint8 }

  compute_force:
    - rebo_compute_force
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

      if (igar_field->empty()) fatal_error() << "igar_force_interp: empty IGAR field (missing igar_read_field?)" << std::endl;
      igar_field_force(*grid, *domain, pecfunc, *igar_field, *energy_factor, scheme);
    }
  };

  template <class GridT> using IgarForceInterpTmpl = IgarForceInterp<GridT>;

  ONIKA_AUTORUN_INIT(igar_force_interp) { OperatorNodeFactory::instance()->register_factory("igar_force_interp", make_grid_variant_operator<IgarForceInterpTmpl>); }

} // namespace exaStamp
