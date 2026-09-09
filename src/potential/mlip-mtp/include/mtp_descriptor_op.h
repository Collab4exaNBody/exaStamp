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

#include <onika/memory/allocator.h>
#include <exanb/core/concurent_add_contributions.h>

#include "../emtp.h"

namespace exaStamp
{

  using namespace exanb;

  // Descriptor-only counterpart to MtpForceOp (mtp_force_op.h): computes the coefficient-free
  // MTP basis-function vector B_k for the central particle and writes it to its own output slot.
  // Unlike PodDescriptorOp, the central atom's type IS needed here (MTP's radial basis, and hence
  // the descriptor itself, is per-(central,neighbor)-species-pair) -- see emtp.h.
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) MtpDescriptorOp
  {
    std::vector<std::shared_ptr<EMTP>>& m_mtp_interfaces;
    const std::vector<int>&             m_type_map;   // exaStamp 0-indexed -> MTP 0-indexed (positional)
    const size_t * const __restrict__   m_cell_particle_offset = nullptr;
    double * const __restrict__         m_descriptors = nullptr;
    // compute_derivative only: [k*3+xyz] -> per-particle aggregate array, or nullptr to skip.
    double * const * const __restrict__ m_deriv_agg_ptrs = nullptr;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () (int jnum, ComputeBufferT& buf, int type, CellParticlesT cells) const
    {
      EMTP& mtp = *m_mtp_interfaces[omp_get_thread_num()];

      mtp.peratom_descriptors_soa(buf.drx, buf.dry, buf.drz, type, buf.ext.type, jnum, m_type_map.data());

      const int K = mtp.alpha_scalar_moments;
      const size_t p = m_cell_particle_offset[buf.cell] + buf.part;
      double * const __restrict__ out = m_descriptors + static_cast<size_t>(K) * p;
      for (int k = 0; k < K; k++) out[k] = mtp.bd[k];

      // Raw per-neighbor-pair Jacobian mtp.bdd[3*jj + 3*jnum*k] = d(B_k)/d(rij), rij = r_neighbor -
      // r_central. Reduce it into a compact per-atom aggregate using the same central+=/neighbor-=
      // scatter idiom as POD/SNAP/k2b.
      if (m_deriv_agg_ptrs != nullptr)
      {
        for (int jj = 0; jj < jnum; jj++)
        {
          size_t nbh_cell=0, nbh_part=0;
          buf.nbh.get(jj, nbh_cell, nbh_part);
          const size_t nbh_p = m_cell_particle_offset[nbh_cell] + nbh_part;
          for (int k = 0; k < K; k++)
          {
            const size_t base = 3*static_cast<size_t>(jj) + 3*static_cast<size_t>(jnum)*k;
            const double vx = mtp.bdd[base+0];
            const double vy = mtp.bdd[base+1];
            const double vz = mtp.bdd[base+2];

            const int k3 = k*3;
            atomic_add_contribution(m_deriv_agg_ptrs[k3+0][p],  vx);
            atomic_add_contribution(m_deriv_agg_ptrs[k3+1][p],  vy);
            atomic_add_contribution(m_deriv_agg_ptrs[k3+2][p],  vz);
            atomic_add_contribution(m_deriv_agg_ptrs[k3+0][nbh_p], -vx);
            atomic_add_contribution(m_deriv_agg_ptrs[k3+1][nbh_p], -vy);
            atomic_add_contribution(m_deriv_agg_ptrs[k3+2][nbh_p], -vz);
          }
        }
      }
    }
  };

}
