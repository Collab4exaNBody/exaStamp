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

#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/grid_cell_particles/grid_cell_values_utils.h>
#include <exanb/core/grid_particle_field_accessor.h>
#include <onika/math/quaternion_operators.h>
#include <onika/cuda/cuda.h>

namespace exanb
{

  namespace ParticleCellProjectionTools
  {
    using namespace GridCellValuesUtils;

    // Per-particle triquadratic Lagrange interpolation of igar_ep.
    // Builds a 3x3x3 subcell stencil around each particle and evaluates
    // both energy and its gradient analytically (C1 continuous).
    template<class LDBGT, class GridT>
    inline void get_particle_force_from_grid(
        LDBGT& ldbg
      , GridT& grid
      , GridCellValues& grid_cell_values
      , double energy_factor )
    {
      static constexpr const char* energy_field_name = "igar_ep";

      if( !grid_cell_values.has_field(energy_field_name) )
      {
        ldbg << "get_particle_force_from_grid: no '" << energy_field_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }

      ldbg << "get_particle_force_from_grid: subcell-level triquadratic Lagrange interpolation" << std::endl;

      const IJK dims = grid.dimension();
      const ssize_t gl = grid.ghost_layers();
      const IJK dims_no_gl = dims - 2*gl;
      const double cell_size = grid.cell_size();

      auto accessor = grid_cell_values.field_data(energy_field_name);
      const double* __restrict__ energy_ptr = accessor.m_data_ptr;
      const size_t energy_stride = accessor.m_stride;
      const ssize_t subdiv = static_cast<ssize_t>( grid_cell_values.field(energy_field_name).m_subdiv );
      const double subcell_size = cell_size / subdiv;
      const double inv_subcell_size = 1.0 / subcell_size;

      auto cells = grid.cells();

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_gl, _cell_flat_i, cell_loc_no_gl, schedule(static))
        {
          const IJK cell_loc = cell_loc_no_gl + gl;
          const size_t cell_i = grid_ijk_to_index(dims, cell_loc);
          const Vec3d cell_origin = grid.cell_position(cell_loc);
          const unsigned int np = cells[cell_i].size();
          const auto* __restrict__ rx = cells[cell_i][field::rx];
          const auto* __restrict__ ry = cells[cell_i][field::ry];
          const auto* __restrict__ rz = cells[cell_i][field::rz];
          auto* __restrict__ fx = cells[cell_i][field::fx];
          auto* __restrict__ fy = cells[cell_i][field::fy];
          auto* __restrict__ fz = cells[cell_i][field::fz];
          auto* __restrict__ ep = cells[cell_i][field::ep];

          for(unsigned int p = 0; p < np; p++)
          {
            // Position relative to cell corner
            const Vec3d rco = { rx[p] - cell_origin.x, ry[p] - cell_origin.y, rz[p] - cell_origin.z };

            // Locate the subcell the particle falls in
            IJK center_cell_loc, center_subcell_loc;
            localize_subcell(rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc);
            center_cell_loc = center_cell_loc + cell_loc;  // relative → absolute

            // Normalized offset from particle's subcell center: u,v,w ∈ [-0.5, 0.5].
            // Lagrange nodes (neighboring subcell centers) are at -1, 0, +1 in these units.
            const double u = rco.x * inv_subcell_size - (center_subcell_loc.i + 0.5);
            const double v = rco.y * inv_subcell_size - (center_subcell_loc.j + 0.5);
            const double w = rco.z * inv_subcell_size - (center_subcell_loc.k + 0.5);

            // Build 3x3x3 stencil of subcell energy values centered on particle's subcell.
            // subcell_neighbor handles cross-cell-boundary indexing transparently.
            // stencil[ck+1][cj+1][ci+1], center at stencil[1][1][1]. Out-of-domain = 0.
            double stencil[3][3][3] = {};
            for(int ck=-1; ck<=1; ck++)
            for(int cj=-1; cj<=1; cj++)
            for(int ci=-1; ci<=1; ci++)
            {
              IJK nbh_cell_loc, nbh_subcell_loc;
              subcell_neighbor(center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc);
              if( grid.contains(nbh_cell_loc) )
              {
                const ssize_t nbh_cell_i    = grid_ijk_to_index(dims, nbh_cell_loc);
                const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv,subdiv,subdiv}, nbh_subcell_loc);
                stencil[ck+1][cj+1][ci+1]  = energy_ptr[ nbh_cell_i * energy_stride + nbh_subcell_i ];
              }
            }

            // Lagrange basis:  L0(t)=t(t-1)/2  L1(t)=1-t^2  L2(t)=t(t+1)/2
            // Derivatives:    dL0=t-1/2        dL1=-2t       dL2=t+1/2
            const double Lu[3]  = { u*(u-1.0)*0.5,  1.0-u*u,  u*(u+1.0)*0.5 };
            const double Lv[3]  = { v*(v-1.0)*0.5,  1.0-v*v,  v*(v+1.0)*0.5 };
            const double Lw[3]  = { w*(w-1.0)*0.5,  1.0-w*w,  w*(w+1.0)*0.5 };
            const double dLu[3] = { u-0.5,  -2.0*u,  u+0.5 };
            const double dLv[3] = { v-0.5,  -2.0*v,  v+0.5 };
            const double dLw[3] = { w-0.5,  -2.0*w,  w+0.5 };

            double ep_interp = 0.0;
            double dEdx = 0.0, dEdy = 0.0, dEdz = 0.0;
            for(int k=0; k<3; k++)
            for(int j=0; j<3; j++)
            for(int i=0; i<3; i++)
            {
              const double e = stencil[k][j][i];
              ep_interp += Lu[i]  * Lv[j]  * Lw[k] * e;
              dEdx      += dLu[i] * Lv[j]  * Lw[k] * e;
              dEdy      += Lu[i]  * dLv[j] * Lw[k] * e;
              dEdz      += Lu[i]  * Lv[j]  * dLw[k] * e;
            }

            ep[p] += energy_factor * ep_interp;
            fx[p] += -energy_factor * dEdx * inv_subcell_size;
            fy[p] += -energy_factor * dEdy * inv_subcell_size;
            fz[p] += -energy_factor * dEdz * inv_subcell_size;
          }
        }
        GRID_OMP_FOR_END;
      }
    }

  } // ParticleCellProjectionTools

} // exanb
