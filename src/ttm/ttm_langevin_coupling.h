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
#pragma once

#include <exanb/core/grid_algorithm.h>
#include <exanb/core/grid_fields.h>
#include <exanb/grid_cell_particles/grid_cell_values_utils.h>
#include <exanb/compute/compute_cell_particles.h>
#include <onika/cuda/cuda.h>

#include "ttm_langevin_rng.h"

namespace exaStamp
{
  using namespace exanb;

  // Pass B of ionic_electronic_heat_transfer: additive Langevin electron-phonon coupling.
  // For each particle: gathers a splat-weighted local Te from the 27-neighbor subcell stencil
  // (same kernel as the Ti deposit), applies an additive Langevin force using Stage 5's stateless
  // RNG, and atomically deposits the work done back into the per-subcell cell_energy_transfer sink.
  // GPU-portable: HasTypeField/HasIdField are resolved at compile time (mirrors get_mass's old
  // tag-dispatch, but via field_pointer_or_null read directly off the captured cell accessor,
  // rather than routing optional fields through compute_cell_particles' FieldSet mechanism).
  template<class CellsT, bool HasTypeField, bool HasIdField>
  struct TtmLangevinCouplingFunctor
  {
    CellsT m_cells;

    Vec3d grid_origin = { 0., 0., 0. };
    IJK grid_offset = { 0, 0, 0 };
    IJK grid_dims = { 0, 0, 0 }; // local dims, including ghost layers
    ssize_t gl = 0;              // ghost layer count, for the is_local_cell (ghost) test
    ssize_t subdiv = 0;

    double cell_size = 0.0;
    double subcell_size = 0.0;
    double sp_size = 0.0;

    double gamma_p = 0.0;
    double gamma_s = 0.0;
    double v_0_sq = 0.0;
    double kB = 0.0;
    double delta_t = 0.0;
    double noise_variance_factor = 2.0;
    bool use_lammps_noise = false;
    uint64_t rng_md_step = 0;

    const double * __restrict__ masses = nullptr;   // persistent GPU-visible species-mass table
    const double * __restrict__ te_ptr = nullptr;    // grid_cell_values "te" field
    size_t te_stride = 0;
    double * __restrict__ cell_energy_transfer_ptr = nullptr; // size n_cells * subdiv^3

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t cell_i, size_t p_i, double rx, double ry, double rz, double vx, double vy, double vz, double& fx, double& fy, double& fz ) const
    {
      using namespace GridCellValuesUtils;
      const ssize_t n_subcells = subdiv * subdiv * subdiv;

      double mass = masses[0];
      if constexpr ( HasTypeField )
      {
        const auto* __restrict__ atom_type = m_cells[cell_i].field_pointer_or_null( field::type );
        mass = masses[ atom_type[p_i] ];
      }
      uint64_t p_id = ( uint64_t(cell_i) << 32 ) | uint64_t(p_i);
      if constexpr ( HasIdField )
      {
        const auto* __restrict__ ids = m_cells[cell_i].field_pointer_or_null( field::id );
        p_id = ids[p_i];
      }

      const IJK cell_loc = grid_index_to_ijk( grid_dims, cell_i );
      const bool is_local_cell = ! inside_grid_shell( grid_dims, 0, gl, cell_loc );
      const Vec3d cell_origin = grid_origin + ( (grid_offset + cell_loc) * cell_size );
      const Vec3d r { rx, ry, rz };
      const Vec3d v { vx, vy, vz };
      const Vec3d rco = r - cell_origin;

      IJK center_cell_loc, center_subcell_loc;
      localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
      center_cell_loc = center_cell_loc + cell_loc;

      // gather local Te (weighted average, same splat kernel as the Ti deposit) and remember
      // per-neighbor weights to re-deposit the coupling work below.
      double nbh_w[27];
      size_t nbh_idx[27];
      int nbh_count = 0;
      double Te_local = 0.0;
      for(int ck=-1;ck<=1;ck++)
      for(int cj=-1;cj<=1;cj++)
      for(int ci=-1;ci<=1;ci++)
      {
        IJK nbh_cell_loc, nbh_subcell_loc;
        subcell_neighbor( center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc );
        if( grid_contains( grid_dims, nbh_cell_loc ) )
        {
          const ssize_t nbh_cell_i = grid_ijk_to_index( grid_dims, nbh_cell_loc );
          const ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
          const Vec3d nbh_cell_origin = grid_origin + ( (grid_offset + nbh_cell_loc) * cell_size );
          const AABB subcell_box = { nbh_cell_origin + nbh_subcell_loc*subcell_size , nbh_cell_origin + (nbh_subcell_loc+1)*subcell_size };
          const double w = particle_weight( r, sp_size, subcell_box );
          const size_t scindex = nbh_cell_i * n_subcells + nbh_subcell_i;
          const double te_val = te_ptr[ nbh_cell_i * te_stride + nbh_subcell_i ];
          Te_local += w * te_val;
          nbh_w[nbh_count] = w;
          nbh_idx[nbh_count] = scindex;
          ++nbh_count;
        }
      }

      if( Te_local > 0.0 )
      {
        // electronic stopping: friction boosted above the v_0 velocity threshold (LAMMPS fix-ttm
        // convention) -- only the friction term is boosted, not the noise term.
        double friction = gamma_p;
        if( gamma_s != 0.0 && norm2(v) > v_0_sq ) { friction += gamma_s; }

        const double noise_x = use_lammps_noise ? ttm_langevin_uniform_rand(p_id,rng_md_step,0) : ttm_langevin_gauss_rand(p_id,rng_md_step,0);
        const double noise_y = use_lammps_noise ? ttm_langevin_uniform_rand(p_id,rng_md_step,1) : ttm_langevin_gauss_rand(p_id,rng_md_step,1);
        const double noise_z = use_lammps_noise ? ttm_langevin_uniform_rand(p_id,rng_md_step,2) : ttm_langevin_gauss_rand(p_id,rng_md_step,2);

        const double noise_amplitude = std::sqrt( noise_variance_factor * kB * gamma_p * Te_local / delta_t );
        const Vec3d f_langevin {
          -friction * v.x + noise_x * noise_amplitude ,
          -friction * v.y + noise_y * noise_amplitude ,
          -friction * v.z + noise_z * noise_amplitude
        };

        // energy actually injected into this particle over the whole MD step by holding f_langevin
        // constant over delta_t -- see ionic_electronic_heat_transfer.cu's Pass B history for why
        // this needs the quadratic (F.F) term, not just the linear (F.v) one.
        const double work = delta_t * dot(f_langevin,v) + 0.5 * delta_t * delta_t * dot(f_langevin,f_langevin) / mass;

        if( is_local_cell )
        {
          fx += f_langevin.x;
          fy += f_langevin.y;
          fz += f_langevin.z;
        }

        for(int k=0;k<nbh_count;k++)
        {
          ONIKA_CU_ATOMIC_ADD( cell_energy_transfer_ptr[ nbh_idx[k] ] , nbh_w[k] * work );
        }
      }
    }
  };

}

namespace exanb
{
  template<class CellsT, bool HasTypeField, bool HasIdField>
  struct ComputeCellParticlesTraits<exaStamp::TtmLangevinCouplingFunctor<CellsT,HasTypeField,HasIdField>>
  {
    static inline constexpr bool CudaCompatible = true;
  };
}
