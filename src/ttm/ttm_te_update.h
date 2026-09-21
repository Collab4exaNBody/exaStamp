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

#include <onika/cuda/cuda.h>
#include <onika/parallel/parallel_for.h>

namespace exaStamp
{
  // Pass E of ionic_electronic_heat_transfer: explicit Te update
  //   dTe = ( Ke*Laplacian(Te) - Te/Ti coupling sink + Se ) / (Ce*rho_e),   Te += dTe * inner_dt
  // one flattened (cell_i, subcell_i) index per thread, dispatched over [0, n_cells*subdiv^3) with
  // onika::parallel::parallel_for (thread-per-index), same shape as TtmLaplacianFunctor (ttm_laplacian.h).
  // Must NOT go through block_parallel_for: that runs the functor once per block with ALL threads of the
  // block on the same index, so this accumulating "Te +=" got applied 1..blockDim times per subcell.
  // Every subcell is independent (dTe reads the Laplacian computed by the previous kernel, never a
  // neighbouring Te), so this is race-free; ghost cells are updated too, same as the old host loop.
  //
  // dTe() is separate from operator() and does not read Te, so the host-side diagnostic reduction
  // that runs after the kernel (sum_dTe) can call it again and get exactly the value the kernel used.
  //
  // Se/Si-style source terms are NOT evaluated here: se_ptr is the array precomputed once per outer
  // MD step (ttm_scratch_se), so no virtual ScalarSourceTerm call reaches device code.
  //
  // Indexing: lap/transfer/Se buffers are tightly packed (idx = cell_i*subdiv^3 + sc), Te is a
  // grid_cell_values field with its own te_stride (see TtmLaplacianFunctor).
  struct TtmTeUpdateFunctor
  {
    ssize_t subdiv = 0;

    double * __restrict__ te_ptr = nullptr;                       // "te" field, te_stride-strided (read+write)
    size_t te_stride = 0;
    const double * __restrict__ lap_te_ptr = nullptr;             // Laplacian of Te, from TtmLaplacianFunctor
    const double * __restrict__ energy_transfer_ptr = nullptr;    // total Te<->Ti energy per subcell over the whole MD step
    const double * __restrict__ se_ptr = nullptr;                 // electronic source term, precomputed per MD step

    double Te_cond = 0.0;       // Ke
    double Ce_rho_e = 1.0;
    double subcell_volume = 1.0;
    double delta_t = 0.0;       // full MD step (the sink is spread over it)
    double inner_dt = 0.0;      // diffusion substep

    ONIKA_HOST_DEVICE_FUNC inline double dTe ( size_t idx ) const
    {
      // energy_transfer is a FIXED total for the whole outer MD step: divide by delta_t (not
      // inner_dt) to get the power density sink applied at every inner step.
      const double coupling_sink = energy_transfer_ptr[idx] / subcell_volume / delta_t;
      return ( Te_cond*lap_te_ptr[idx] - coupling_sink + se_ptr[idx] ) / Ce_rho_e;
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t idx ) const
    {
      const size_t n_subcells = size_t(subdiv) * size_t(subdiv) * size_t(subdiv);
      const size_t cell_i = idx / n_subcells;
      const size_t scindex = idx % n_subcells;
      te_ptr[ cell_i*te_stride + scindex ] += dTe(idx) * inner_dt;
    }
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
  }
}
