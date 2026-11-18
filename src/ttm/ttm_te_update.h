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
#include <onika/cuda/cuda.h>
#include <onika/parallel/parallel_for.h>
#include <exanb/grid_cell_particles/reduce_grid_cell_values.h>

#include "ttm_piecewise_poly.h"

namespace exaStamp
{
  using namespace exanb;

  // Pass E of ionic_electronic_heat_transfer: explicit Te update
  //   dTe = ( Ke*Laplacian(Te) - Te/Ti coupling sink + Se ) / (Ce*rho_e),   Te += dTe * inner_dt
  // or, with a Te-dependent Ce(Te), in energy form: E(Te) += power * inner_dt, E(Te) = integral_0^Te Ce dT
  // one flattened (local cell, subcell) index per thread, dispatched over [0, n_local_cells*subdiv^3) with
  // onika::parallel::parallel_for (thread-per-index), same shape as TtmLaplacianFunctor (ttm_laplacian.h).
  // Must NOT go through block_parallel_for: that runs the functor once per block with ALL threads of the
  // block on the same index, so this accumulating "Te +=" got applied 1..blockDim times per subcell.
  // Every subcell is independent (dTe reads the Laplacian computed by the previous kernel, never a
  // neighbouring Te), so this is race-free; only local cells are updated (ghost Te is refreshed by a ghost exchange between substeps).
  //
  // power()/dTe() are separate from operator() and do not read Te, so the host-side diagnostic reduction
  // that runs after the kernel (sum_dTe) can call them again and get exactly the value the kernel used.
  //
  // Se/Si-style source terms are NOT evaluated here: se_ptr is the array precomputed once per outer
  // MD step (ttm_scratch_se), so no virtual ScalarSourceTerm call reaches device code.
  //
  // Indexing: lap/transfer/Se buffers are tightly packed (idx = cell_i*subdiv^3 + sc), Te is a
  // grid_cell_values field with its own te_stride (see TtmLaplacianFunctor).
  struct TtmTeUpdateFunctor
  {
    IJK grid_dims = { 0, 0, 0 }; // local dims, including ghost layers
    ssize_t ghost_layers = 0;    // only local (non-ghost) cells are updated: ghost Te comes from a ghost exchange
    ssize_t subdiv = 0;

    double * __restrict__ te_ptr = nullptr;                       // "te" field, te_stride-strided (read+write)
    size_t te_stride = 0;
    const double * __restrict__ lap_te_ptr = nullptr;             // Laplacian of Te, from TtmLaplacianFunctor
    const double * __restrict__ energy_transfer_ptr = nullptr;    // total Te<->Ti energy per subcell over the whole MD step
    const double * __restrict__ se_ptr = nullptr;                 // electronic source term, precomputed per MD step

    double Te_cond = 0.0;       // Ke
    // Reciprocals precomputed on the host: this kernel is FP64-pipe bound on GPUs with a weak FP64 rate
    // (ncu: 89% FP64 pipe, 38 DFMA/thread, almost all from 3 FP64 divisions), so no division here.
    double inv_Ce_rho_e = 1.0;  // 1 / (Ce*rho_e)
    double sink_scale = 1.0;    // 1 / (subcell_volume * delta_t), delta_t = full MD step (the sink is spread over it)
    double inner_dt = 0.0;      // diffusion substep
    // Te-dependent heat capacity Ce(Te) (disabled: constant inv_Ce_rho_e above). When enabled, the update is
    // done in energy form: the energy density E = integral_0^Te Ce dT is advanced by power()*inner_dt and Te
    // is recovered by inverting E(Te), so the electronic energy changes by exactly the applied power (an
    // explicit Te += power/Ce(Te) step gains 1/2 Ce'(Te) dTe^2 per step and subcell).
    TtmPiecewisePoly ce;

    // power density of the Te equation, Ke*Laplacian(Te) - Te/Ti coupling sink + Se (does not read Te)
    ONIKA_HOST_DEVICE_FUNC inline double power ( size_t idx ) const
    {
      // energy_transfer is a FIXED total for the whole outer MD step: scale by 1/(volume*delta_t) (not
      // inner_dt) to get the power density sink applied at every inner step.
      const double coupling_sink = energy_transfer_ptr[idx] * sink_scale;
      return Te_cond*lap_te_ptr[idx] - coupling_sink + se_ptr[idx];
    }

    // constant Ce: Te rate (does not read Te)
    ONIKA_HOST_DEVICE_FUNC inline double dTe ( size_t idx ) const
    {
      return power( idx ) * inv_Ce_rho_e;
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t idx ) const
    {
      // idx enumerates local cells x subcells; the scratch buffers (lap/transfer/Se) and dTe() keep the
      // full-grid layout (cell_i*n_subcells + subcell), ghosts included
      const size_t n_subcells = size_t(subdiv) * size_t(subdiv) * size_t(subdiv);
      const IJK cell_loc = grid_index_to_ijk( grid_dims - 2*ghost_layers , ssize_t(idx / n_subcells) ) + ghost_layers;
      const size_t cell_i = size_t( grid_ijk_to_index( grid_dims, cell_loc ) );
      const size_t scindex = idx % n_subcells;
      const size_t j = cell_i*n_subcells + scindex;
      double& te = te_ptr[ cell_i*te_stride + scindex ];
      if( ce.enabled() )
      {
        const double dE = power( j ) * inner_dt;
        const double c0 = ce.value( te );
        const double guess = te + ( c0 > 0.0 ? dE / c0 : 0.0 ); // explicit Te step: Newton converges in a few iterations
        te = ce.inverse_integral( ce.integral( te ) + dE , guess );
      }
      else
      {
        te += dTe( j ) * inner_dt;
      }
    }
  };

  // Per-step Te statistics for the Te-dependent parameters, one reduction over local cells
  // (exanb::reduce_grid_cell_values): energy density integral of Ce, Te range (validity warning), min Ce (stability).
  struct TtmTeStats
  {
    double sum_te = 0.0;
    double sum_energy = 0.0;   // sum of integral_0^Te Ce dT (0 if Ce is constant)
    double min_te = 1.e300;
    double max_te = -1.e300;
    double min_ce = 1.e300;    // stays 1e300 if Ce is constant
  };

  struct TtmTeStatsFunc
  {
    TtmPiecewisePoly ce;       // disabled: constant Ce, only the Te range is reduced
    ONIKA_HOST_DEVICE_FUNC inline void accumulate( TtmTeStats& acc , const double* v ) const
    {
      const double te = v[0];
      acc.sum_te += te;
      acc.min_te = te < acc.min_te ? te : acc.min_te; // no std::min in device code
      acc.max_te = te > acc.max_te ? te : acc.max_te;
      if( ce.enabled() )
      {
        acc.sum_energy += ce.integral( te );
        const double c = ce.value( te );
        acc.min_ce = c < acc.min_ce ? c : acc.min_ce;
      }
    }
    ONIKA_HOST_DEVICE_FUNC inline void combine( TtmTeStats& acc , const TtmTeStats& p ) const
    {
      acc.sum_te += p.sum_te;
      acc.sum_energy += p.sum_energy;
      acc.min_te = p.min_te < acc.min_te ? p.min_te : acc.min_te;
      acc.max_te = p.max_te > acc.max_te ? p.max_te : acc.max_te;
      acc.min_ce = p.min_ce < acc.min_ce ? p.min_ce : acc.min_ce;
    }
  };

  // Device-side fill of a managed scratch buffer: a host-side std::vector::assign() would touch the
  // pages on the CPU and force a UVM migration + kernel stall on the next GPU kernel that reads them.
  struct TtmFillFunctor
  {
    double * __restrict__ ptr = nullptr;
    double value = 0.0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t idx ) const { ptr[idx] = value; }
  };
}

namespace onika
{
  namespace parallel
  {
    template<>
    struct ParallelForFunctorTraits<exaStamp::TtmTeUpdateFunctor>
    {
      static inline constexpr bool CudaCompatible = true;
    };

    template<>
    struct ParallelForFunctorTraits<exaStamp::TtmFillFunctor>
    {
      static inline constexpr bool CudaCompatible = true;
    };
  }
}

namespace exanb
{
  template<> struct ReduceGridCellValuesTraits<exaStamp::TtmTeStatsFunc> { static inline constexpr bool CudaCompatible = true; };
}
