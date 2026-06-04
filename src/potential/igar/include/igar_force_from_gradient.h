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

    // Nearest-subcell lookup of precomputed igar_dEdx/y/z gradient fields.
    // Requires igar_compute_gradient to have been called first.
    template<class LDBGT, class GridT>
    inline void get_particle_force_from_gradient_grid(
        LDBGT& ldbg
      , GridT& grid
      , const GridCellValues& grid_cell_values
      , double energy_factor )
    {
      static constexpr const char* dEdx_name = "igar_dEdx";
      static constexpr const char* dEdy_name = "igar_dEdy";
      static constexpr const char* dEdz_name = "igar_dEdz";
      static constexpr const char* e_name = "igar_ep";

      if( !grid_cell_values.has_field(dEdx_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << dEdx_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      if( !grid_cell_values.has_field(dEdy_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << dEdy_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      if( !grid_cell_values.has_field(dEdz_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << dEdz_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      if( !grid_cell_values.has_field(e_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << e_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      const GridCellField& dEdx_field = grid_cell_values.field(dEdx_name);
      const ssize_t subdiv = static_cast<ssize_t>( dEdx_field.m_subdiv );
      assert( subdiv > 0 );
#     ifndef NDEBUG
      {
        const GridCellField& e_field    = grid_cell_values.field(e_name);
        const GridCellField& dEdy_field = grid_cell_values.field(dEdy_name);
        const GridCellField& dEdz_field = grid_cell_values.field(dEdz_name);
        const ssize_t subcell_count = subdiv * subdiv * subdiv;
        assert( e_field.m_subdiv    == static_cast<size_t>(subdiv) );
        assert( dEdy_field.m_subdiv == static_cast<size_t>(subdiv) );
        assert( dEdz_field.m_subdiv == static_cast<size_t>(subdiv) );
        assert( e_field.m_components    >= static_cast<size_t>(subcell_count) );
        assert( dEdx_field.m_components >= static_cast<size_t>(subcell_count) );
        assert( dEdy_field.m_components >= static_cast<size_t>(subcell_count) );
        assert( dEdz_field.m_components >= static_cast<size_t>(subcell_count) );
        assert( e_field.m_components    % static_cast<size_t>(subcell_count) == 0 );
        assert( dEdx_field.m_components % static_cast<size_t>(subcell_count) == 0 );
        assert( dEdy_field.m_components % static_cast<size_t>(subcell_count) == 0 );
        assert( dEdz_field.m_components % static_cast<size_t>(subcell_count) == 0 );
        assert( e_field.m_components    / static_cast<size_t>(subcell_count) == 1 );
        assert( dEdx_field.m_components / static_cast<size_t>(subcell_count) == 1 );
        assert( dEdy_field.m_components / static_cast<size_t>(subcell_count) == 1 );
        assert( dEdz_field.m_components / static_cast<size_t>(subcell_count) == 1 );
      }
#     endif

      const IJK dims = grid.dimension();
      const ssize_t gl = grid.ghost_layers();
      const IJK dims_no_gl = dims - 2*gl;
      const double cell_size = grid.cell_size();
      const double subcell_size = cell_size / subdiv;
      const double inv_subcell_size = 1.0 / subcell_size;

      auto e_accessor = grid_cell_values.field_data(e_name);
      auto dEdx_accessor = grid_cell_values.field_data(dEdx_name);
      auto dEdy_accessor = grid_cell_values.field_data(dEdy_name);
      auto dEdz_accessor = grid_cell_values.field_data(dEdz_name);
      const double* __restrict__ e_ptr = e_accessor.m_data_ptr;
      const double* __restrict__ dEdx_ptr = dEdx_accessor.m_data_ptr;
      const double* __restrict__ dEdy_ptr = dEdy_accessor.m_data_ptr;
      const double* __restrict__ dEdz_ptr = dEdz_accessor.m_data_ptr;
      const size_t e_stride = e_accessor.m_stride;
      const size_t dEdx_stride = dEdx_accessor.m_stride;
      const size_t dEdy_stride = dEdy_accessor.m_stride;
      const size_t dEdz_stride = dEdz_accessor.m_stride;

      auto cells = grid.cells();

      ldbg << "get_particle_force_from_gradient_grid: apply force from igar_dEdx/igar_dEdy/igar_dEdz"
           << " subdiv=" << subdiv << " energy_factor=" << energy_factor << std::endl;

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
          auto* __restrict__ e = cells[cell_i][field::ep];
          auto* __restrict__ fx = cells[cell_i][field::fx];
          auto* __restrict__ fy = cells[cell_i][field::fy];
          auto* __restrict__ fz = cells[cell_i][field::fz];

          for(unsigned int p = 0; p < np; p++)
          {
            const Vec3d rco = { rx[p] - cell_origin.x, ry[p] - cell_origin.y, rz[p] - cell_origin.z };

            IJK center_cell_loc, center_subcell_loc;
            localize_subcell(rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc);
            center_cell_loc = center_cell_loc + cell_loc;

            if( grid.contains(center_cell_loc) )
            {
              const ssize_t center_cell_i = grid_ijk_to_index(dims, center_cell_loc);
              const ssize_t center_subcell_i = grid_ijk_to_index(IJK{subdiv,subdiv,subdiv}, center_subcell_loc);

              fx[p] += -energy_factor * dEdx_ptr[center_cell_i * dEdx_stride + center_subcell_i];
              fy[p] += -energy_factor * dEdy_ptr[center_cell_i * dEdy_stride + center_subcell_i];
              fz[p] += -energy_factor * dEdz_ptr[center_cell_i * dEdz_stride + center_subcell_i];

              // Trilinear interpolation of ep between the 8 surrounding subcell centers.
              // u,v,w ∈ [-0.5, 0.5]: normalized offset from the particle's subcell center.
              const double u = rco.x * inv_subcell_size - (center_subcell_loc.i + 0.5);
              const double v = rco.y * inv_subcell_size - (center_subcell_loc.j + 0.5);
              const double w = rco.z * inv_subcell_size - (center_subcell_loc.k + 0.5);
              const int di = (u >= 0.0) ? 1 : -1;
              const int dj = (v >= 0.0) ? 1 : -1;
              const int dk = (w >= 0.0) ? 1 : -1;
              const double tu = (u >= 0.0) ? u : (1.0 + u);
              const double tv = (v >= 0.0) ? v : (1.0 + v);
              const double tw = (w >= 0.0) ? w : (1.0 + w);
              double ep_interp = 0.0;
              for(int bk=0; bk<=1; bk++)
              for(int bj=0; bj<=1; bj++)
              for(int bi=0; bi<=1; bi++)
              {
                IJK nbh_cell_loc, nbh_subcell_loc;
                subcell_neighbor(center_cell_loc, center_subcell_loc, subdiv,
                                 IJK{bi*di, bj*dj, bk*dk}, nbh_cell_loc, nbh_subcell_loc);
                if( grid.contains(nbh_cell_loc) )
                {
                  const ssize_t nbh_cell_i    = grid_ijk_to_index(dims, nbh_cell_loc);
                  const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv,subdiv,subdiv}, nbh_subcell_loc);
                  const double wx = (bi == 0) ? (1.0 - tu) : tu;
                  const double wy = (bj == 0) ? (1.0 - tv) : tv;
                  const double wz = (bk == 0) ? (1.0 - tw) : tw;
                  ep_interp += wx * wy * wz * e_ptr[nbh_cell_i * e_stride + nbh_subcell_i];
                }
              }
              e[p] += ep_interp;
            }
          }
        }
        GRID_OMP_FOR_END;
      }
    }

  } // ParticleCellProjectionTools

} // exanb
