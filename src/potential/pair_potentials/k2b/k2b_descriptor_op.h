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

#include <cmath>

#include <onika/memory/allocator.h>
#include <exanb/core/concurent_add_contributions.h>

#include "potential.h"

namespace exaStamp
{
  using namespace exanb;

  // Descriptor-only counterpart to k2b's own SymetricForceOp (pair_potential_singlemat_symetric.cpp):
  // computes the coefficient-free per-atom 2-body kernel descriptor D_i[k] = sum_j g_k(r_ij), where
  // g_k is the k-th Gaussian RBF from k2b_potential_compute_force (potential.h) with w/Delta
  // stripped out, plus its compact per-atom derivative aggregate -- same central+=/neighbor-=
  // scatter idiom as POD's PodDescriptorOp (mlip-pod/include/pod_descriptor_op.h), but computed
  // directly here since k2b has no external descriptor library (no EAPOD equivalent).
  //
  // Sign convention: buf.drx/dry/drz point from the central particle to the neighbor (same
  // convention SymetricForceOp relies on: F_central = +de*dr, F_neighbor = -de*dr for
  // F = -dE/dx). For a plain (non-force) derivative d(g_k)/dx there is no extra "-dE/dx" sign
  // flip, so the central/neighbor contributions come out negated relative to that force pattern:
  // d(g_k)/dx_central = -dgdr*dr/r, d(g_k)/dx_neighbor = +dgdr*dr/r.
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) K2bDescriptorOp
  {
    K2bPotentialParametersRO             m_params;
    const size_t * const __restrict__    m_cell_particle_offset = nullptr;
    double * const __restrict__          m_descriptors = nullptr;
    // compute_derivative only: [k*3+xyz] -> per-particle aggregate array, or nullptr to skip.
    double * const * const __restrict__  m_deriv_agg_ptrs = nullptr;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () (int jnum, ComputeBufferT& buf, CellParticlesT cells) const
    {
      const int    K        = m_params.n_rbf;
      const double r_min    = m_params.r_min;
      const double r_cut    = m_params.r_cut;
      const double sigma    = m_params.sigma;
      const double inv_2sig2 = 0.5 / (sigma * sigma);
      const double inv_sig2  = 1.0 / (sigma * sigma);
      const double h = (K > 1) ? (r_cut - r_min) / static_cast<double>(K - 1) : 0.0;

      const size_t p = m_cell_particle_offset[buf.cell] + buf.part;
      double * const __restrict__ out = m_descriptors + static_cast<size_t>(K) * p;

      for (int jj = 0; jj < jnum; jj++)
      {
        const double drx = buf.drx[jj], dry = buf.dry[jj], drz = buf.drz[jj];
        const double r = std::sqrt(buf.d2[jj]);
        if (r <= 0.0) continue;

        size_t nbh_cell=0, nbh_part=0;
        buf.nbh.get(jj, nbh_cell, nbh_part);
        const size_t nbh_p = m_cell_particle_offset[nbh_cell] + nbh_part;

        for (int k = 0; k < K; k++)
        {
          const double s_k = r_min + k * h;
          const double d   = r - s_k;               // displacement from grid centre
          const double arg = d * d * inv_2sig2;
          if (arg > 20.0) continue;

          const double g = std::exp(-arg);
          out[k] += g;

          if (m_deriv_agg_ptrs != nullptr)
          {
            const double dgdr = -d * inv_sig2 * g;   // d(g_k)/dr
            const double c = -dgdr / r;              // see sign-convention note above
            const double vx = c * drx, vy = c * dry, vz = c * drz;
            atomic_add_contribution(m_deriv_agg_ptrs[k*3+0][p],  vx);
            atomic_add_contribution(m_deriv_agg_ptrs[k*3+1][p],  vy);
            atomic_add_contribution(m_deriv_agg_ptrs[k*3+2][p],  vz);
            atomic_add_contribution(m_deriv_agg_ptrs[k*3+0][nbh_p], -vx);
            atomic_add_contribution(m_deriv_agg_ptrs[k*3+1][nbh_p], -vy);
            atomic_add_contribution(m_deriv_agg_ptrs[k*3+2][nbh_p], -vz);
          }
        }
      }
    }
  };
}
