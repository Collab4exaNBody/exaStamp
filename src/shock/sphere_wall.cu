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
  struct SphereWallComputeFunc
  {
    const Vec3d center;
    const double radius;
    const double R;
    const long exposant = 12;
    const double epsilon = onika::physics::make_quantity(1.0e-19, "J").convert();
    XFormT xform;

    ONIKA_HOST_DEVICE_FUNC inline void operator()(double rx, double ry, double rz, double &fx, double &fy, double &fz, double &ep) const
    {
      Vec3d r{rx, ry, rz};
      r = xform.transformCoord(r);

      Vec3d rc = r - center;
      const double dist = norm(rc);

      double d_sign = dist - radius;
      double d = abs(d_sign);

      if (d <= R && dist > 0.)
      {
        const Vec3d N = rc / dist;

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
  struct ComputeCellParticlesTraits<exaStamp::SphereWallComputeFunc<XFormT>>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template <typename GridT, class = AssertGridHasFields<GridT, field::_fx, field::_fy, field::_fz, field::_ep>>
  class SphereWall : public OperatorNode
  {
    static inline constexpr double default_epsilon = ONIKA_CONST_QUANTITY(1.0e-19 * J).convert(exaStamp::UNIT_SYSTEM);

    ADD_SLOT(GridT, grid, INPUT_OUTPUT);
    ADD_SLOT(Vec3d, center, INPUT, Vec3d{0.0, 0.0, 0.0});
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
Applies a repulsive potential wall to every particle: a spherical barrier defined by
`center` and `radius`, repelling particles within `cutoff` of the sphere's surface
along a power-law of exponent `exponent`, scaled by `epsilon`.

For a particle at signed distance d = |r - center| - radius (|d| <= cutoff), pushed
along the radial direction N = (r - center) / |r - center|:
ep += epsilon * (1 - cutoff/d)^exponent
f  = -epsilon * exponent * (cutoff/d^2) * (1 - cutoff/d)^(exponent - 1) * N

Particles just inside the sphere are pushed further inward, particles just outside
are pushed further outward (a thin spherical shell/membrane barrier), same convention
as the planar `wall` operator.

YAML example:

myoperator:
  - sphere_wall:
      center: [0.0, 0.0, 0.0]
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
        SphereWallComputeFunc<LinearXForm> func{*center, *radius, *cutoff, *exponent, *epsilon, LinearXForm{domain->xform()}};
        compute_cell_particles(*grid, false, func, compute_field_set, parallel_execution_context());
      }
      else
      {
        SphereWallComputeFunc<NullXForm> func{*center, *radius, *cutoff, *exponent, *epsilon, NullXForm{}};
        compute_cell_particles(*grid, false, func, compute_field_set, parallel_execution_context());
      }
    }
  };

  template <class GridT>
  using SphereWallTmpl = SphereWall<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(sphere_wall)
  {
    OperatorNodeFactory::instance()->register_factory("sphere_wall", make_grid_variant_operator<SphereWallTmpl>);
  }

}
