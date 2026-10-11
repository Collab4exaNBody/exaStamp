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
#include <onika/cuda/cuda.h>
#include <onika/parallel/parallel_for.h>

namespace exaStamp
{
  using namespace exanb;

  // Pass D of ionic_electronic_heat_transfer: 27-point discrete Laplacian of Te, one flattened
  // (local cell, subcell) index per thread, dispatched over [0, n_local_cells*subdiv^3) via
  // onika::parallel::parallel_for (thread-per-index). NOT block_parallel_for: that one runs the
  // functor once per BLOCK, with every thread of the block executing the SAME index -- harmless for a
  // pure overwrite like this one, but it silently multi-applies any accumulating kernel (see
  // ttm_te_update.h). GPU opt-in trait is therefore onika::parallel::ParallelForFunctorTraits.
  // GPU-portable: only reads Te
  // + this struct's by-value grid geometry -- uses the grid_contains/grid_index_to_ijk free
  // functions and GridCellValuesUtils::subcell_neighbor (all ONIKA_HOST_DEVICE_FUNC) instead of the
  // host-only Grid::contains/gcv_subcell_neighbor the old CPU-only loop used.
  //
  // te_stride vs subdiv^3: Te is a grid_cell_values FIELD, whose per-cell stride
  // (cell_te_data.m_stride) is the total component count across every field ever added to that
  // GridCellValues container, not just "te"'s own subdiv^3 -- they only coincide when "te" is the
  // sole field. The old host loop conflated the two (read Te via a bare subdiv^3-based index),
  // a real bug whenever grid_cell_values holds any other field alongside "te"; fixed here by taking
  // te_stride as its own parameter, matching the convention every other Te access in this plugin
  // already uses (Pass B/E, init_ttm.cpp). lap_te_ptr has no such ambiguity: it's a dedicated
  // scratch buffer (ttm_scratch_lap_te), always tightly packed to exactly n_cells*subdiv^3.
  struct TtmLaplacianFunctor
  {
    IJK grid_dims = { 0, 0, 0 }; // local dims, including ghost layers
    ssize_t ghost_layers = 0;    // layers skipped on each side: only local (non-ghost) cells are computed, ghost Te comes from a ghost exchange
    ssize_t subdiv = 0;

    // 1/subcell_size^2, precomputed on the host: an FP64 division costs ~10 DFMA-equivalents per
    // thread, and this kernel is FP64-pipe bound on GPUs with a weak FP64 rate (ncu: 89% FP64 pipe).
    double inv_subcell_size_sq = 0.0;

    const double * __restrict__ te_ptr = nullptr; // grid_cell_values "te" field, te_stride-strided
    size_t te_stride = 0;

    double * __restrict__ lap_te_ptr = nullptr; // size n_cells * subdiv^3, tightly packed

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t idx ) const
    {
      using namespace GridCellValuesUtils;
      const ssize_t n_subcells = subdiv * subdiv * subdiv;
      // idx enumerates local cells x subcells; cell_i is the index in the full grid (ghosts included),
      // and lap_te_ptr keeps the full-grid layout (cell_i*n_subcells + subcell)
      const IJK cell_loc = grid_index_to_ijk( grid_dims - 2*ghost_layers , ssize_t(idx) / n_subcells ) + ghost_layers;
      const ssize_t cell_i = grid_ijk_to_index( grid_dims, cell_loc );
      const ssize_t sc_i = ssize_t(idx) % n_subcells;
      const IJK sc = grid_index_to_ijk( IJK{subdiv,subdiv,subdiv} , sc_i );

      // inspired from https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
      static constexpr double Lap27Norm = 26.0;
      static constexpr double Lap27_compact[4] = { -88/Lap27Norm , 6/Lap27Norm, 3/Lap27Norm, 2/Lap27Norm };

      double L_Te = 0.0;
      for(int nk=-1;nk<=1;nk++)
      for(int nj=-1;nj<=1;nj++)
      for(int ni=-1;ni<=1;ni++)
      {
        IJK nbh_cell_loc, nbh_subcell_loc;
        subcell_neighbor( cell_loc, sc, subdiv, IJK{ni,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
        if( grid_contains( grid_dims, nbh_cell_loc ) )
        {
          const ssize_t nbh_cell_i = grid_ijk_to_index( grid_dims, nbh_cell_loc );
          const ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
          const size_t nbh_j = size_t(nbh_cell_i) * te_stride + size_t(nbh_subcell_i);
          // manual abs: nvcc silently miscompiles std::abs called from device code the same way
          // std::min/std::max do (constexpr-host-only stdlib trap) -- see feedback memory.
          const int lap_compact_index = (ni<0?-ni:ni) + (nj<0?-nj:nj) + (nk<0?-nk:nk);
          L_Te += te_ptr[nbh_j] * Lap27_compact[lap_compact_index];
        }
      }
      // normalize by h^2: raw stencil sum = h^2 * laplacian(Te) + O(h^4)
      lap_te_ptr[ cell_i*n_subcells + sc_i ] = L_Te * inv_subcell_size_sq;
    }
  };

}

namespace onika
{
  namespace parallel
  {
    template<>
    struct ParallelForFunctorTraits<exaStamp::TtmLaplacianFunctor>
    {
      static inline constexpr bool CudaCompatible = true;
    };
  }
}
