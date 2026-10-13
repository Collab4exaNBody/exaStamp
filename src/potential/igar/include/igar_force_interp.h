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
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/grid_cell_particles/grid_cell_values_utils.h>
#include <onika/cuda/cuda.h>
#include <onika/math/quaternion_operators.h>

namespace exaStamp
{

  using namespace exanb;

  // Interpolation schemes of the igar_ep subcell values (values sit at the subcell centers).
  //  - Lagrange   : triquadratic Lagrange on the 3x3x3 stencil centered on the particle's subcell. Interpolating, but the
  //                 stencil moves when the particle crosses a subcell face: energy discontinuous there.
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

  // t in [-1/2, 1/2]: offset from the particle's subcell center, in subcell units. Nodes at -1, 0, +1.
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
  // particle's subcell (t >= 0) or the one below it (t < 0).
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

  template <IgarInterpolation Scheme> struct CellScalarEnergyToParticleForceFunctor
  {
    using Weights = IgarInterpolationWeights<Scheme>;
    static inline constexpr int N = Weights::N;

    Vec3d grid_origin = {0., 0., 0.};
    IJK grid_offset = {0, 0, 0};
    IJK grid_dims = {0, 0, 0};
    ssize_t subdiv = 0;

    double cell_size = 0.0;
    double subcell_size = 0.0;
    double inv_subcell_size = 0.0;

    double energy_factor = 0.0;

    const double *__restrict__ e_ptr = nullptr;

    size_t e_stride = 0;

    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t cell_i, size_t p_i, double rx, double ry, double rz, double &fx, double &fy, double &fz, double &en) const
    {
      using GridCellValuesUtils::localize_subcell;
      using GridCellValuesUtils::subcell_neighbor;
      const IJK cell_loc = grid_index_to_ijk(grid_dims, cell_i);
      const Vec3d cell_origin = grid_origin + ((grid_offset + cell_loc) * cell_size);
      const Vec3d rco = {rx - cell_origin.x, ry - cell_origin.y, rz - cell_origin.z};

      // between two neighbor list updates, a particle may lie outside its cell: localize_subcell then returns
      // the cell it actually lies in, relative to cell_loc
      IJK rel_cell_loc, center_subcell_loc;
      localize_subcell(rco, cell_size, subcell_size, subdiv, rel_cell_loc, center_subcell_loc);
      const IJK center_cell_loc = rel_cell_loc + cell_loc;

      if (grid_contains(grid_dims, center_cell_loc))
      {
        // offset from the center of the particle's subcell, in subcell units (in [-1/2, 1/2])
        const Vec3d rsc = rco - (rel_cell_loc * cell_size);
        const double u = rsc.x * inv_subcell_size - (center_subcell_loc.i + 0.5);
        const double v = rsc.y * inv_subcell_size - (center_subcell_loc.j + 0.5);
        const double w = rsc.z * inv_subcell_size - (center_subcell_loc.k + 0.5);

        double wu[N], wv[N], ww[N], dwu[N], dwv[N], dww[N];
        Weights::weights(u, wu, dwu);
        Weights::weights(v, wv, dwv);
        Weights::weights(w, ww, dww);
        const IJK first = {Weights::first_node(u), Weights::first_node(v), Weights::first_node(w)};

        // stencil values read through subcell_neighbor (crosses cell boundaries), out-of-grid values count as 0
        double ep_interp = 0.0;
        double dEdx = 0.0, dEdy = 0.0, dEdz = 0.0;
        for (int k = 0; k < N; k++)
          for (int j = 0; j < N; j++)
            for (int i = 0; i < N; i++)
            {
              IJK nbh_cell_loc, nbh_subcell_loc;
              subcell_neighbor(center_cell_loc, center_subcell_loc, subdiv, first + IJK{i, j, k}, nbh_cell_loc, nbh_subcell_loc);
              if (grid_contains(grid_dims, nbh_cell_loc))
              {
                const ssize_t nbh_cell_i = grid_ijk_to_index(grid_dims, nbh_cell_loc);
                const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                const double e = e_ptr[nbh_cell_i * e_stride + nbh_subcell_i];
                ep_interp += wu[i] * wv[j] * ww[k] * e;
                dEdx += dwu[i] * wv[j] * ww[k] * e;
                dEdy += wu[i] * dwv[j] * ww[k] * e;
                dEdz += wu[i] * wv[j] * dww[k] * e;
              }
            }

        en += energy_factor * ep_interp;
        fx += -energy_factor * dEdx * inv_subcell_size;
        fy += -energy_factor * dEdy * inv_subcell_size;
        fz += -energy_factor * dEdz * inv_subcell_size;
      }
    }
  };

} // namespace exaStamp

namespace exanb
{
  template <exaStamp::IgarInterpolation Scheme> struct ComputeCellParticlesTraits<exaStamp::CellScalarEnergyToParticleForceFunctor<Scheme>>
  {
    static inline constexpr bool CudaCompatible = true;
  };
} // namespace exanb

namespace exaStamp
{
  namespace ParticleCellProjectionTools
  {
    using namespace GridCellValuesUtils;

    // Per-particle interpolation of igar_ep (scheme chosen at run time, see IgarInterpolation): energy and its
    // analytical gradient at the particle position.
    template <class LDBGT, class GridT, class PECFuncT>
    inline void get_particle_force_from_grid(LDBGT &ldbg, GridT &grid, const PECFuncT &parallel_execution_context, const GridCellValues &grid_cell_values, double energy_factor, IgarInterpolation scheme)
    {
      static constexpr const char *e_name = "igar_ep";

      if (!grid_cell_values.has_field(e_name))
      {
        ldbg << "get_particle_force_from_grid: no '" << e_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }

      const GridCellField &e_field = grid_cell_values.field(e_name);
      const ssize_t subdiv = static_cast<ssize_t>(e_field.m_subdiv);

      const double cell_size = grid.cell_size();
      const double subcell_size = cell_size / subdiv;
      const double inv_subcell_size = 1.0 / subcell_size;

      auto e_accessor = grid_cell_values.field_data(e_name);
      const double *__restrict__ e_ptr = e_accessor.m_data_ptr;
      const size_t e_stride = e_accessor.m_stride;

      auto cp_fields = onika::make_flat_tuple(grid.field_accessor(field::rx), grid.field_accessor(field::ry), grid.field_accessor(field::rz), grid.field_accessor(field::fx), grid.field_accessor(field::fy), grid.field_accessor(field::fz), grid.field_accessor(field::ep));

      auto run = [&](auto scheme_tag)
      {
        constexpr IgarInterpolation S = decltype(scheme_tag)::value;
        CellScalarEnergyToParticleForceFunctor<S> func = {grid.origin(), grid.offset(), grid.dimension(), subdiv, cell_size, subcell_size, inv_subcell_size, energy_factor, e_ptr, e_stride};
        compute_cell_particles(grid, false, func, cp_fields, parallel_execution_context());
      };

      switch (scheme)
      {
        case IgarInterpolation::Lagrange: run(std::integral_constant<IgarInterpolation, IgarInterpolation::Lagrange>{}); break;
        case IgarInterpolation::BSpline: run(std::integral_constant<IgarInterpolation, IgarInterpolation::BSpline>{}); break;
        case IgarInterpolation::CatmullRom: run(std::integral_constant<IgarInterpolation, IgarInterpolation::CatmullRom>{}); break;
      }
    }

  } // namespace ParticleCellProjectionTools

} // namespace exaStamp
