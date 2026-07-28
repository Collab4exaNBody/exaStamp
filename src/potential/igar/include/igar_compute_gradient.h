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

#include <exanb/core/grid_particle_field_accessor.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/grid_cell_particles/grid_cell_values_utils.h>
#include <onika/cuda/cuda.h>
#include <onika/math/quaternion_operators.h>

namespace exanb
{

  namespace ParticleCellProjectionTools
  {

    using namespace GridCellValuesUtils;

    template <class LDBGT, class GridT> inline void compute_energy_gradient(LDBGT &ldbg, GridT &grid, GridCellValues &grid_cell_values)
    {
      static constexpr const char *energy_field_name = "igar_ep";
      static constexpr const char *dEdx_name = "igar_dEdx";
      static constexpr const char *dEdy_name = "igar_dEdy";
      static constexpr const char *dEdz_name = "igar_dEdz";

      if (!grid_cell_values.has_field(energy_field_name))
      {
        ldbg << "compute_energy_gradient: no '" << energy_field_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }

      const IJK dims = grid.dimension();
      const ssize_t gl = grid.ghost_layers();
      const IJK dims_no_gl = dims - 2 * gl;
      const double cell_size = grid.cell_size();

      const GridCellField &energy_field = grid_cell_values.field(energy_field_name);
      const ssize_t subdiv = static_cast<ssize_t>(energy_field.m_subdiv);
      const ssize_t subcell_count = subdiv * subdiv * subdiv;
      assert(subdiv > 0);
      assert(energy_field.m_components >= static_cast<size_t>(subcell_count));
      assert(energy_field.m_components % static_cast<size_t>(subcell_count) == 0);
      const size_t energy_comps = energy_field.m_components / static_cast<size_t>(subcell_count);

      if (!grid_cell_values.has_field(dEdx_name))
      {
        grid_cell_values.add_field(dEdx_name, subdiv, 1);
      }
      if (!grid_cell_values.has_field(dEdy_name))
      {
        grid_cell_values.add_field(dEdy_name, subdiv, 1);
      }
      if (!grid_cell_values.has_field(dEdz_name))
      {
        grid_cell_values.add_field(dEdz_name, subdiv, 1);
      }

      const GridCellField &dEdx_field = grid_cell_values.field(dEdx_name);
      const GridCellField &dEdy_field = grid_cell_values.field(dEdy_name);
      const GridCellField &dEdz_field = grid_cell_values.field(dEdz_name);
      assert(dEdx_field.m_subdiv == static_cast<size_t>(subdiv));
      assert(dEdy_field.m_subdiv == static_cast<size_t>(subdiv));
      assert(dEdz_field.m_subdiv == static_cast<size_t>(subdiv));
      assert(dEdx_field.m_components >= static_cast<size_t>(subcell_count));
      assert(dEdy_field.m_components >= static_cast<size_t>(subcell_count));
      assert(dEdz_field.m_components >= static_cast<size_t>(subcell_count));
      assert(dEdx_field.m_components % static_cast<size_t>(subcell_count) == 0);
      assert(dEdy_field.m_components % static_cast<size_t>(subcell_count) == 0);
      assert(dEdz_field.m_components % static_cast<size_t>(subcell_count) == 0);
      const size_t dEdx_comps = dEdx_field.m_components / static_cast<size_t>(subcell_count);
      const size_t dEdy_comps = dEdy_field.m_components / static_cast<size_t>(subcell_count);
      const size_t dEdz_comps = dEdz_field.m_components / static_cast<size_t>(subcell_count);

      auto energy_accessor = grid_cell_values.field_data(energy_field_name);
      auto dEdx_accessor = grid_cell_values.field_data(dEdx_name);
      auto dEdy_accessor = grid_cell_values.field_data(dEdy_name);
      auto dEdz_accessor = grid_cell_values.field_data(dEdz_name);
      const double *__restrict__ energy_ptr = energy_accessor.m_data_ptr;
      double *__restrict__ dEdx_ptr = dEdx_accessor.m_data_ptr;
      double *__restrict__ dEdy_ptr = dEdy_accessor.m_data_ptr;
      double *__restrict__ dEdz_ptr = dEdz_accessor.m_data_ptr;
      const size_t energy_stride = energy_accessor.m_stride;
      const size_t dEdx_stride = dEdx_accessor.m_stride;
      const size_t dEdy_stride = dEdy_accessor.m_stride;
      const size_t dEdz_stride = dEdz_accessor.m_stride;
      const double subcell_size = cell_size / subdiv;
      const double inv_2_subcell_size = 0.5 / subcell_size;

      ldbg << "compute_energy_gradient: subcell-centered finite differences on '" << energy_field_name << "' subdiv=" << subdiv << std::endl;

#pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN (dims_no_gl, _cell_flat_i, cell_loc_no_gl, schedule(static))
        {
          const IJK cell_loc = cell_loc_no_gl + gl;
          const size_t cell_i = grid_ijk_to_index(dims, cell_loc);

          for (ssize_t ck = 0; ck < subdiv; ck++)
            for (ssize_t cj = 0; cj < subdiv; cj++)
              for (ssize_t ci = 0; ci < subdiv; ci++)
              {
                const IJK subcell_loc{ci, cj, ck};
                const ssize_t subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, subcell_loc);

                double e_xm = 0.0, e_xp = 0.0;
                double e_ym = 0.0, e_yp = 0.0;
                double e_zm = 0.0, e_zp = 0.0;

                IJK nbh_cell_loc, nbh_subcell_loc;
                subcell_neighbor(cell_loc, subcell_loc, subdiv, IJK{-1, 0, 0}, nbh_cell_loc, nbh_subcell_loc);
                if (grid.contains(nbh_cell_loc))
                {
                  const ssize_t nbh_cell_i = grid_ijk_to_index(dims, nbh_cell_loc);
                  const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                  e_xm = energy_ptr[nbh_cell_i * energy_stride + nbh_subcell_i * energy_comps];
                }
                subcell_neighbor(cell_loc, subcell_loc, subdiv, IJK{1, 0, 0}, nbh_cell_loc, nbh_subcell_loc);
                if (grid.contains(nbh_cell_loc))
                {
                  const ssize_t nbh_cell_i = grid_ijk_to_index(dims, nbh_cell_loc);
                  const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                  e_xp = energy_ptr[nbh_cell_i * energy_stride + nbh_subcell_i * energy_comps];
                }

                subcell_neighbor(cell_loc, subcell_loc, subdiv, IJK{0, -1, 0}, nbh_cell_loc, nbh_subcell_loc);
                if (grid.contains(nbh_cell_loc))
                {
                  const ssize_t nbh_cell_i = grid_ijk_to_index(dims, nbh_cell_loc);
                  const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                  e_ym = energy_ptr[nbh_cell_i * energy_stride + nbh_subcell_i * energy_comps];
                }
                subcell_neighbor(cell_loc, subcell_loc, subdiv, IJK{0, 1, 0}, nbh_cell_loc, nbh_subcell_loc);
                if (grid.contains(nbh_cell_loc))
                {
                  const ssize_t nbh_cell_i = grid_ijk_to_index(dims, nbh_cell_loc);
                  const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                  e_yp = energy_ptr[nbh_cell_i * energy_stride + nbh_subcell_i * energy_comps];
                }

                subcell_neighbor(cell_loc, subcell_loc, subdiv, IJK{0, 0, -1}, nbh_cell_loc, nbh_subcell_loc);
                if (grid.contains(nbh_cell_loc))
                {
                  const ssize_t nbh_cell_i = grid_ijk_to_index(dims, nbh_cell_loc);
                  const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                  e_zm = energy_ptr[nbh_cell_i * energy_stride + nbh_subcell_i * energy_comps];
                }
                subcell_neighbor(cell_loc, subcell_loc, subdiv, IJK{0, 0, 1}, nbh_cell_loc, nbh_subcell_loc);
                if (grid.contains(nbh_cell_loc))
                {
                  const ssize_t nbh_cell_i = grid_ijk_to_index(dims, nbh_cell_loc);
                  const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                  e_zp = energy_ptr[nbh_cell_i * energy_stride + nbh_subcell_i * energy_comps];
                }

                dEdx_ptr[cell_i * dEdx_stride + subcell_i * dEdx_comps] = (e_xp - e_xm) * inv_2_subcell_size;
                dEdy_ptr[cell_i * dEdy_stride + subcell_i * dEdy_comps] = (e_yp - e_ym) * inv_2_subcell_size;
                dEdz_ptr[cell_i * dEdz_stride + subcell_i * dEdz_comps] = (e_zp - e_zm) * inv_2_subcell_size;
              }
        }
        GRID_OMP_FOR_END;
      }
    }

  } // namespace ParticleCellProjectionTools

} // namespace exanb
