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
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/core/domain.h>

#include <onika/cuda/cuda.h>

#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/compute_pair_optional_args.h>

#include <exaStamp/unit_system.h>

#include <cmath>

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  template <class XFormT>
  struct CylinderWallComputeFunc
  {
    const Vec3d origin;
    const Vec3d axis; // unit vector
    const double radius;
    const double R;
    const long exposant = 12;
    const double epsilon = onika::physics::make_quantity(1.0e-19, "J").convert();
    XFormT xform;

    ONIKA_HOST_DEVICE_FUNC inline void operator()(double rx, double ry, double rz, double &fx, double &fy, double &fz, double &ep) const
    {
      Vec3d r{rx, ry, rz};
      r = xform.transformCoord(r);

      Vec3d ra = r - origin;
      Vec3d r_perp = ra - dot(ra, axis) * axis;
      const double dist = norm(r_perp);

      double d_sign = dist - radius;
      double d = abs(d_sign);

      if (d <= R && dist > 0.)
      {
        const Vec3d N = r_perp / dist;

        const double ratio = 1.0 - R / d;
        const double ratio_puis_exposant = pow(ratio, exposant);

        double f_contrib = -epsilon * exposant * (R / (d * d)) * ratio_puis_exposant / ratio;
        Vec3d F = N * f_contrib;
        if (d_sign < 0.)
          F = -F;

        // Energy
        ep += epsilon * ratio_puis_exposant;

        // Forces
        fx += F.x;
        fy += F.y;
        fz += F.z;
      }
    }
  };

}

namespace exanb
{
  template <class XFormT>
  struct ComputeCellParticlesTraits<exaStamp::CylinderWallComputeFunc<XFormT>>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template <typename GridT, class = AssertGridHasFields<GridT, field::_fx, field::_fy, field::_fz, field::_ep>>
  class CylinderWall : public OperatorNode
  {
    static inline constexpr double default_epsilon = ONIKA_CONST_QUANTITY(1.0e-19 * J).convert(exaStamp::UNIT_SYSTEM);

    ADD_SLOT(GridT, grid, INPUT_OUTPUT);
    ADD_SLOT(Vec3d, origin, INPUT, Vec3d{0.0, 0.0, 0.0});
    ADD_SLOT(Vec3d, axis, INPUT, Vec3d{0.0, 0.0, 1.0});
    ADD_SLOT(double, radius, INPUT, REQUIRED);
    ADD_SLOT(double, cutoff, INPUT, REQUIRED);
    ADD_SLOT(Domain, domain, INPUT, REQUIRED);
    ADD_SLOT(long, exponent, INPUT, 12);
    ADD_SLOT(double, epsilon, INPUT, default_epsilon);

    static constexpr FieldSet<field::_rx, field::_ry, field::_rz, field::_fx, field::_fy, field::_fz, field::_ep> compute_field_set{};

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Applies a repulsive potential wall to every particle: a cylindrical barrier defined by
`origin` (a point on the axis), `axis` (unit vector, cylinder direction) and `radius`,
repelling particles within `cutoff` of the cylinder's surface along a power-law of
exponent `exponent`, scaled by `epsilon`. `axis` must be a unit vector.

For a particle at radial vector r_perp = (r - origin) projected perpendicular to axis,
signed distance d = |r_perp| - radius (|d| <= cutoff), pushed along N = r_perp / |r_perp|:
ep += epsilon * (1 - cutoff/d)^exponent
f  = -epsilon * exponent * (cutoff/d^2) * (1 - cutoff/d)^(exponent - 1) * N

Particles just inside the cylinder are pushed further inward, particles just outside
are pushed further outward (a thin cylindrical shell/membrane barrier), same convention
as the planar `wall` operator. The cylinder is infinite along axis (unbounded ends).

YAML example:

myoperator:
  - cylinder_wall:
      origin: [0.0, 0.0, 0.0]
      axis: [0.0, 0.0, 1.0]
      radius: 20.0 ang
      cutoff: 5.0 ang
      epsilon: 1.0e-19 J
      exponent: 12
)EOF";
    }

    inline void execute() override final
    {
      if (!domain->xform_is_identity())
      {
        CylinderWallComputeFunc<LinearXForm> func{*origin, *axis, *radius, *cutoff, *exponent, *epsilon, LinearXForm{domain->xform()}};
        compute_cell_particles(*grid, false, func, compute_field_set, parallel_execution_context());
      }
      else
      {
        CylinderWallComputeFunc<NullXForm> func{*origin, *axis, *radius, *cutoff, *exponent, *epsilon, NullXForm{}};
        compute_cell_particles(*grid, false, func, compute_field_set, parallel_execution_context());
      }
    }
  };

  template <class GridT>
  using CylinderWallTmpl = CylinderWall<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(cylinder_wall)
  {
    OperatorNodeFactory::instance()->register_factory("cylinder_wall", make_grid_variant_operator<CylinderWallTmpl>);
  }

}
