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
#pragma once

#include <exanb/compute/compute_cell_particles.h>
#include <exanb/core/grid_particle_field_accessor.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types_operators.h>
#include <onika/cuda/cuda.h>
#include <onika/math/quaternion_operators.h>

#include "igar_field.h"

namespace exaStamp
{

  using namespace exanb;

  // Interpolation schemes of the IGAR field (values sit at the voxel centers).
  //  - Lagrange   : triquadratic Lagrange on the 3x3x3 stencil centered on the particle's voxel. Interpolating, but the
  //                 stencil moves when the particle crosses a voxel face: energy discontinuous there.
  //  - BSpline    : quadratic B-spline on the same 3x3x3 stencil. Energy and force continuous (C1), force = exact
  //                 derivative of the energy. Not interpolating: smooths the field (error ~ h^2/8 * curvature).
  //  - CatmullRom : tricubic Catmull-Rom on a 4x4x4 stencil. Interpolating and C1, force = exact derivative.
  enum class IgarInterpolation : int
  {
    Lagrange = 0,
    BSpline = 1,
    CatmullRom = 2
  };

  template <IgarInterpolation Scheme> struct IgarInterpolationWeights;

  // t in [-1/2, 1/2]: offset from the particle's voxel center, in voxel units. Nodes at -1, 0, +1.
  template <> struct IgarInterpolationWeights<IgarInterpolation::Lagrange>
  {
    static inline constexpr int N = 3;
    ONIKA_HOST_DEVICE_FUNC static inline int first_node(double) { return -1; }
    ONIKA_HOST_DEVICE_FUNC static inline void weights(double t, double w[N], double dw[N])
    {
      w[0] = t * (t - 1.0) * 0.5;
      w[1] = 1.0 - t * t;
      w[2] = t * (t + 1.0) * 0.5;
      dw[0] = t - 0.5;
      dw[1] = -2.0 * t;
      dw[2] = t + 0.5;
    }
  };

  template <> struct IgarInterpolationWeights<IgarInterpolation::BSpline>
  {
    static inline constexpr int N = 3;
    ONIKA_HOST_DEVICE_FUNC static inline int first_node(double) { return -1; }
    ONIKA_HOST_DEVICE_FUNC static inline void weights(double t, double w[N], double dw[N])
    {
      w[0] = 0.5 * (0.5 - t) * (0.5 - t);
      w[1] = 0.75 - t * t;
      w[2] = 0.5 * (0.5 + t) * (0.5 + t);
      dw[0] = t - 0.5;
      dw[1] = -2.0 * t;
      dw[2] = t + 0.5;
    }
  };

  // The particle lies between nodes 0 and 1 of the 4 nodes -1, 0, 1, 2 (s in [0,1] from node 0): node 0 is the
  // particle's voxel (t >= 0) or the one below it (t < 0).
  template <> struct IgarInterpolationWeights<IgarInterpolation::CatmullRom>
  {
    static inline constexpr int N = 4;
    ONIKA_HOST_DEVICE_FUNC static inline int first_node(double t) { return (t >= 0.0) ? -1 : -2; }
    ONIKA_HOST_DEVICE_FUNC static inline void weights(double t, double w[N], double dw[N])
    {
      const double s = (t >= 0.0) ? t : (1.0 + t);
      const double s2 = s * s, s3 = s2 * s;
      w[0] = 0.5 * (-s3 + 2.0 * s2 - s);
      w[1] = 0.5 * (3.0 * s3 - 5.0 * s2 + 2.0);
      w[2] = 0.5 * (-3.0 * s3 + 4.0 * s2 + s);
      w[3] = 0.5 * (s3 - s2);
      dw[0] = 0.5 * (-3.0 * s2 + 4.0 * s - 1.0);
      dw[1] = 0.5 * (9.0 * s2 - 10.0 * s);
      dw[2] = 0.5 * (-9.0 * s2 + 8.0 * s + 1.0);
      dw[3] = 0.5 * (3.0 * s2 - 2.0 * s);
    }
  };

  // Wraps a stencil node index along one axis: periodic axes wrap around, nodes outside a non-periodic axis are -1
  ONIKA_HOST_DEVICE_FUNC inline ssize_t igar_wrap_index(ssize_t n, ssize_t dim, bool periodic)
  {
    if (periodic) return ((n % dim) + dim) % dim;
    return (n >= 0 && n < dim) ? n : -1;
  }

  template <IgarInterpolation Scheme, class ValueT> struct IgarFieldToParticleForceFunctor
  {
    using Weights = IgarInterpolationWeights<Scheme>;
    static inline constexpr int N = Weights::N;

    // the field, the voxels and the particle positions live in grid space (real position = xform * grid position)
    Vec3d origin = {0., 0., 0.};    // domain origin, where voxel (0,0,0) starts
    Vec3d inv_h = {0., 0., 0.};     // 1 / voxel size, per axis
    Mat3d inv_xform_t = {1., 0., 0., 0., 1., 0., 0., 0., 1.}; // transpose(inv_xform): grid-space gradient -> real-space gradient
    IJK dims = {0, 0, 0};           // field dimensions
    bool periodic_x = true, periodic_y = true, periodic_z = true;
    double energy_factor = 0.0;     // includes the field value_scale
    const ValueT *__restrict__ field = nullptr;

    ONIKA_HOST_DEVICE_FUNC inline void operator()(double rx, double ry, double rz, double &fx, double &fy, double &fz, double &en) const
    {
      // position in voxel units; the particle's voxel is c, u is the offset from its center (in [-1/2, 1/2[)
      const double gx = (rx - origin.x) * inv_h.x;
      const double gy = (ry - origin.y) * inv_h.y;
      const double gz = (rz - origin.z) * inv_h.z;
      const double cx = floor(gx), cy = floor(gy), cz = floor(gz);
      const double u = gx - cx - 0.5;
      const double v = gy - cy - 0.5;
      const double w = gz - cz - 0.5;

      double wu[N], wv[N], ww[N], dwu[N], dwv[N], dww[N];
      Weights::weights(u, wu, dwu);
      Weights::weights(v, wv, dwv);
      Weights::weights(w, ww, dww);

      ssize_t ix[N], iy[N], iz[N];
      for (int n = 0; n < N; n++)
      {
        ix[n] = igar_wrap_index(ssize_t(cx) + Weights::first_node(u) + n, dims.i, periodic_x);
        iy[n] = igar_wrap_index(ssize_t(cy) + Weights::first_node(v) + n, dims.j, periodic_y);
        iz[n] = igar_wrap_index(ssize_t(cz) + Weights::first_node(w) + n, dims.k, periodic_z);
      }

      // nodes outside a non-periodic domain count as 0
      double ep_interp = 0.0;
      double dEdx = 0.0, dEdy = 0.0, dEdz = 0.0;
      for (int k = 0; k < N; k++)
        for (int j = 0; j < N; j++)
          for (int i = 0; i < N; i++)
          {
            if (ix[i] >= 0 && iy[j] >= 0 && iz[k] >= 0)
            {
              const double e = static_cast<double>(field[ix[i] + dims.i * (iy[j] + dims.j * iz[k])]);
              ep_interp += wu[i] * wv[j] * ww[k] * e;
              dEdx += dwu[i] * wv[j] * ww[k] * e;
              dEdy += wu[i] * dwv[j] * ww[k] * e;
              dEdz += wu[i] * wv[j] * dww[k] * e;
            }
          }

      // force = -dE/dR with R = xform * r : -transpose(inv_xform) * dE/dr
      const Vec3d f = inv_xform_t * Vec3d{dEdx * inv_h.x, dEdy * inv_h.y, dEdz * inv_h.z};
      en += energy_factor * ep_interp;
      fx += -energy_factor * f.x;
      fy += -energy_factor * f.y;
      fz += -energy_factor * f.z;
    }
  };

} // namespace exaStamp

namespace exanb
{
  template <exaStamp::IgarInterpolation Scheme, class ValueT> struct ComputeCellParticlesTraits<exaStamp::IgarFieldToParticleForceFunctor<Scheme, ValueT>>
  {
    static inline constexpr bool CudaCompatible = true;
  };
} // namespace exanb

namespace exaStamp
{
  // Per-particle interpolation of the IGAR field (scheme chosen at run time, see IgarInterpolation): energy and its
  // analytical gradient at the particle position.
  template <class GridT, class PECFuncT>
  inline void igar_field_force(GridT &grid, const Domain &domain, const PECFuncT &parallel_execution_context, const IgarField &igar_field, double energy_factor, IgarInterpolation scheme)
  {
    const Vec3d size = bounds_size(domain.bounds());
    const IJK dims = igar_field.m_dims;
    const Vec3d inv_h = {dims.i / size.x, dims.j / size.y, dims.k / size.z};
    const double scaled_energy_factor = energy_factor * igar_field.m_value_scale;
    const Mat3d inv_xform_t = transpose(domain.inv_xform());

    auto cp_fields = onika::make_flat_tuple(grid.field_accessor(field::rx), grid.field_accessor(field::ry), grid.field_accessor(field::rz), grid.field_accessor(field::fx), grid.field_accessor(field::fy), grid.field_accessor(field::fz), grid.field_accessor(field::ep));

    auto run = [&](auto scheme_tag, const auto *values)
    {
      constexpr IgarInterpolation S = decltype(scheme_tag)::value;
      using ValueT = std::remove_cv_t<std::remove_pointer_t<decltype(values)>>;
      IgarFieldToParticleForceFunctor<S, ValueT> func = {domain.origin(), inv_h, inv_xform_t, dims, domain.periodic_boundary_x(), domain.periodic_boundary_y(), domain.periodic_boundary_z(), scaled_energy_factor, values};
      compute_cell_particles(grid, false, func, cp_fields, parallel_execution_context());
    };
    auto run_scheme = [&](const auto *values)
    {
      switch (scheme)
      {
        case IgarInterpolation::Lagrange: run(std::integral_constant<IgarInterpolation, IgarInterpolation::Lagrange>{}, values); break;
        case IgarInterpolation::BSpline: run(std::integral_constant<IgarInterpolation, IgarInterpolation::BSpline>{}, values); break;
        case IgarInterpolation::CatmullRom: run(std::integral_constant<IgarInterpolation, IgarInterpolation::CatmullRom>{}, values); break;
      }
    };

    if (igar_field.is_uint8()) run_scheme(igar_field.m_u8.data());
    else run_scheme(igar_field.m_f64.data());
  }

} // namespace exaStamp
