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

  struct CellScalarEnergyToParticleForceFunctor
  {
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

      IJK center_cell_loc, center_subcell_loc;
      localize_subcell(rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc);
      center_cell_loc = center_cell_loc + cell_loc;

      if (grid_contains(grid_dims, center_cell_loc))
      {
        const ssize_t center_cell_i = grid_ijk_to_index(grid_dims, center_cell_loc);
        const ssize_t center_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, center_subcell_loc);

        // Normalized offset from particle's subcell center: u,v,w ∈ [-0.5, 0.5].
        // Lagrange nodes (neighboring subcell centers) are at -1, 0, +1 in these units.
        const double u = rco.x * inv_subcell_size - (center_subcell_loc.i + 0.5);
        const double v = rco.y * inv_subcell_size - (center_subcell_loc.j + 0.5);
        const double w = rco.z * inv_subcell_size - (center_subcell_loc.k + 0.5);
        // Build 3x3x3 stencil of subcell energy values centered on particle's subcell.
        // subcell_neighbor handles cross-cell-boundary indexing transparently.
        // stencil[ck+1][cj+1][ci+1], center at stencil[1][1][1]. Out-of-domain = 0.
        double stencil[3][3][3] = {};
        for (int ck = -1; ck <= 1; ck++)
          for (int cj = -1; cj <= 1; cj++)
            for (int ci = -1; ci <= 1; ci++)
            {
              IJK nbh_cell_loc, nbh_subcell_loc;
              subcell_neighbor(center_cell_loc, center_subcell_loc, subdiv, IJK{ci, cj, ck}, nbh_cell_loc, nbh_subcell_loc);
              if (grid_contains(grid_dims, nbh_cell_loc))
              {
                const ssize_t nbh_cell_i = grid_ijk_to_index(grid_dims, nbh_cell_loc);
                const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                stencil[ck + 1][cj + 1][ci + 1] = e_ptr[nbh_cell_i * e_stride + nbh_subcell_i];
              }
            }

        // Lagrange basis:  L0(t)=t(t-1)/2  L1(t)=1-t^2  L2(t)=t(t+1)/2
        // Derivatives:    dL0=t-1/2        dL1=-2t       dL2=t+1/2
        const double Lu[3] = {u * (u - 1.0) * 0.5, 1.0 - u * u, u * (u + 1.0) * 0.5};
        const double Lv[3] = {v * (v - 1.0) * 0.5, 1.0 - v * v, v * (v + 1.0) * 0.5};
        const double Lw[3] = {w * (w - 1.0) * 0.5, 1.0 - w * w, w * (w + 1.0) * 0.5};
        const double dLu[3] = {u - 0.5, -2.0 * u, u + 0.5};
        const double dLv[3] = {v - 0.5, -2.0 * v, v + 0.5};
        const double dLw[3] = {w - 0.5, -2.0 * w, w + 0.5};

        double ep_interp = 0.0;
        double dEdx = 0.0, dEdy = 0.0, dEdz = 0.0;
        for (int k = 0; k < 3; k++)
          for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
            {
              const double e = stencil[k][j][i];
              ep_interp += Lu[i] * Lv[j] * Lw[k] * e;
              dEdx += dLu[i] * Lv[j] * Lw[k] * e;
              dEdy += Lu[i] * dLv[j] * Lw[k] * e;
              dEdz += Lu[i] * Lv[j] * dLw[k] * e;
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
  template <> struct ComputeCellParticlesTraits<exaStamp::CellScalarEnergyToParticleForceFunctor>
  {
    static inline constexpr bool CudaCompatible = true;
  };
} // namespace exanb

namespace exaStamp
{
  namespace ParticleCellProjectionTools
  {
    using namespace GridCellValuesUtils;

    // Per-particle triquadratic Lagrange interpolation of igar_ep.
    // Builds a 3x3x3 subcell stencil around each particle and evaluates
    // both energy and its gradient analytically (C1 continuous).
    template <class LDBGT, class GridT, class PECFuncT> inline void get_particle_force_from_grid(LDBGT &ldbg, GridT &grid, const PECFuncT &parallel_execution_context, const GridCellValues &grid_cell_values, double energy_factor)
    {
      static constexpr const char *e_name = "igar_ep";

      if (!grid_cell_values.has_field(e_name))
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << e_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }

      const GridCellField &e_field = grid_cell_values.field(e_name);
      const ssize_t subdiv = static_cast<ssize_t>(e_field.m_subdiv);
      const ssize_t subcell_count = subdiv * subdiv * subdiv;

      const IJK dims = grid.dimension();
      const ssize_t gl = grid.ghost_layers();
      const IJK dims_no_gl = dims - 2 * gl;
      const double cell_size = grid.cell_size();
      const double subcell_size = cell_size / subdiv;
      const double inv_subcell_size = 1.0 / subcell_size;

      auto e_accessor = grid_cell_values.field_data(e_name);
      const double *__restrict__ e_ptr = e_accessor.m_data_ptr;
      const size_t e_stride = e_accessor.m_stride;

      auto cells = grid.cells();

      CellScalarEnergyToParticleForceFunctor func = {grid.origin(), grid.offset(), grid.dimension(), subdiv, cell_size, subcell_size, inv_subcell_size, energy_factor, e_ptr, e_stride};

      auto cp_fields = onika::make_flat_tuple(grid.field_accessor(field::rx), grid.field_accessor(field::ry), grid.field_accessor(field::rz), grid.field_accessor(field::fx), grid.field_accessor(field::fy), grid.field_accessor(field::fz), grid.field_accessor(field::ep));
      compute_cell_particles(grid, false, func, cp_fields, parallel_execution_context());

      /*
      #pragma omp parallel
            {
              GRID_OMP_FOR_BEGIN(dims_no_gl, _cell_flat_i, cell_loc_no_gl, schedule(static))
              {
                const IJK cell_loc = cell_loc_no_gl + gl;
                const size_t cell_i = grid_ijk_to_index(dims, cell_loc);
                const Vec3d cell_origin = grid.cell_position(cell_loc);
                const unsigned int np = cells[cell_i].size();
                const auto *__restrict__ rx = cells[cell_i][field::rx];
                const auto *__restrict__ ry = cells[cell_i][field::ry];
                const auto *__restrict__ rz = cells[cell_i][field::rz];
                auto *__restrict__ fx = cells[cell_i][field::fx];
                auto *__restrict__ fy = cells[cell_i][field::fy];
                auto *__restrict__ fz = cells[cell_i][field::fz];
                auto *__restrict__ ep = cells[cell_i][field::ep];

                for (unsigned int p = 0; p < np; p++)
                {
                  // Position relative to cell corner
                  const Vec3d rco = {rx[p] - cell_origin.x, ry[p] - cell_origin.y, rz[p] - cell_origin.z};

                  // Locate the subcell the particle falls in
                  IJK center_cell_loc, center_subcell_loc;
                  localize_subcell(rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc);
                  center_cell_loc = center_cell_loc + cell_loc; // relative → absolute

                  // Normalized offset from particle's subcell center: u,v,w ∈ [-0.5, 0.5].
                  // Lagrange nodes (neighboring subcell centers) are at -1, 0, +1 in these units.
                  const double u = rco.x * inv_subcell_size - (center_subcell_loc.i + 0.5);
                  const double v = rco.y * inv_subcell_size - (center_subcell_loc.j + 0.5);
                  const double w = rco.z * inv_subcell_size - (center_subcell_loc.k + 0.5);

                  // Build 3x3x3 stencil of subcell energy values centered on particle's subcell.
                  // subcell_neighbor handles cross-cell-boundary indexing transparently.
                  // stencil[ck+1][cj+1][ci+1], center at stencil[1][1][1]. Out-of-domain = 0.
                  double stencil[3][3][3] = {};
                  for (int ck = -1; ck <= 1; ck++)
                    for (int cj = -1; cj <= 1; cj++)
                      for (int ci = -1; ci <= 1; ci++)
                      {
                        IJK nbh_cell_loc, nbh_subcell_loc;
                        subcell_neighbor(center_cell_loc, center_subcell_loc, subdiv, IJK{ci, cj, ck}, nbh_cell_loc, nbh_subcell_loc);
                        if (grid.contains(nbh_cell_loc))
                        {
                          const ssize_t nbh_cell_i = grid_ijk_to_index(dims, nbh_cell_loc);
                          const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                          stencil[ck + 1][cj + 1][ci + 1] = energy_ptr[nbh_cell_i * energy_stride + nbh_subcell_i];
                        }
                      }

                  // Lagrange basis:  L0(t)=t(t-1)/2  L1(t)=1-t^2  L2(t)=t(t+1)/2
                  // Derivatives:    dL0=t-1/2        dL1=-2t       dL2=t+1/2
                  const double Lu[3] = {u * (u - 1.0) * 0.5, 1.0 - u * u, u * (u + 1.0) * 0.5};
                  const double Lv[3] = {v * (v - 1.0) * 0.5, 1.0 - v * v, v * (v + 1.0) * 0.5};
                  const double Lw[3] = {w * (w - 1.0) * 0.5, 1.0 - w * w, w * (w + 1.0) * 0.5};
                  const double dLu[3] = {u - 0.5, -2.0 * u, u + 0.5};
                  const double dLv[3] = {v - 0.5, -2.0 * v, v + 0.5};
                  const double dLw[3] = {w - 0.5, -2.0 * w, w + 0.5};

                  double ep_interp = 0.0;
                  double dEdx = 0.0, dEdy = 0.0, dEdz = 0.0;
                  for (int k = 0; k < 3; k++)
                    for (int j = 0; j < 3; j++)
                      for (int i = 0; i < 3; i++)
                      {
                        const double e = stencil[k][j][i];
                        ep_interp += Lu[i] * Lv[j] * Lw[k] * e;
                        dEdx += dLu[i] * Lv[j] * Lw[k] * e;
                        dEdy += Lu[i] * dLv[j] * Lw[k] * e;
                        dEdz += Lu[i] * Lv[j] * dLw[k] * e;
                      }

                  ep[p] += energy_factor * ep_interp;
                  fx[p] += -energy_factor * dEdx * inv_subcell_size;
                  fy[p] += -energy_factor * dEdy * inv_subcell_size;
                  fz[p] += -energy_factor * dEdz * inv_subcell_size;
                }
              }
              GRID_OMP_FOR_END;
            }

  */
    }

  } // namespace ParticleCellProjectionTools

} // namespace exaStamp
