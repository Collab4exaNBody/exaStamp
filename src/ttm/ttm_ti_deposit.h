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
#include <exanb/grid_cell_particles/grid_cell_values_utils.h>
#include <exanb/compute/compute_cell_particles.h>
#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // Pass A of ionic_electronic_heat_transfer: deposits each particle's kinetic energy
  // (mass*v^2), splat-weighted over the 27-neighbor subcell stencil, into the per-subcell
  // Ti scratch buffer. GPU-portable: only reads field values + this struct's by-value grid
  // geometry (no host-side Grid member calls -- uses the grid_contains/grid_ijk_to_index free
  // functions instead, same convention as igar_force_from_gradient.h's particle-to-grid functor).
  struct TtmTiDepositFunctor
  {
    Vec3d grid_origin = { 0., 0., 0. };
    IJK grid_offset = { 0, 0, 0 };
    IJK grid_dims = { 0, 0, 0 }; // local dims, including ghost layers -- same space Ti is sized in
    ssize_t subdiv = 0;

    double cell_size = 0.0;
    double subcell_size = 0.0;
    double subcell_volume = 0.0;
    double sp_size = 0.0;

    // Points into a persistent GPU-visible CudaMMVector<double> (see ionic_electronic_heat_transfer.cu),
    // not captured by value: embedding a MAX_PARTICLE_TYPES-sized array by value blows past onika's
    // GPU kernel-functor size cap (confirmed: 2280 bytes vs. a 1016-byte limit).
    const double * __restrict__ masses = nullptr;

    double * __restrict__ ti_ptr = nullptr; // size n_cells * subdiv^3

    ONIKA_HOST_DEVICE_FUNC inline void deposit( size_t cell_i, double rx, double ry, double rz, double vx, double vy, double vz, double mass ) const
    {
      using namespace GridCellValuesUtils;
      const ssize_t n_subcells = subdiv * subdiv * subdiv;
      const IJK cell_loc = grid_index_to_ijk( grid_dims, cell_i );
      const Vec3d cell_origin = grid_origin + ( (grid_offset + cell_loc) * cell_size );
      const Vec3d r { rx, ry, rz };
      const double v2 = vx*vx + vy*vy + vz*vz;
      const Vec3d rco = r - cell_origin;

      IJK center_cell_loc, center_subcell_loc;
      localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
      center_cell_loc = center_cell_loc + cell_loc;

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
          const double contrib = ( mass * w * v2 ) / subcell_volume;
          ONIKA_CU_ATOMIC_ADD( ti_ptr[scindex] , contrib );
        }
      }
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t cell_i, double rx, double ry, double rz, double vx, double vy, double vz, uint8_t type ) const
    {
      deposit( cell_i, rx, ry, rz, vx, vy, vz, masses[type] );
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t cell_i, double rx, double ry, double rz, double vx, double vy, double vz ) const
    {
      deposit( cell_i, rx, ry, rz, vx, vy, vz, masses[0] );
    }
  };

}

namespace exanb
{
  template<>
  struct ComputeCellParticlesTraits<exaStamp::TtmTiDepositFunctor>
  {
    static inline constexpr bool CudaCompatible = true;
  };
}
