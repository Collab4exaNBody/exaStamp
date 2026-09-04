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

#include "eapod.h"

namespace exaStamp
{

  using namespace exanb;

  // Descriptor-only counterpart to PodForceOp (pod_force_op.h): computes the coefficient-free
  // POD base descriptor vector for the central particle and writes it to its own output slot --
  // no neighbor scatter, no locks needed (matches why compute_descriptor_snap's non-derivative
  // path also needs none).
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PodDescriptorOp
  {
    std::vector<std::shared_ptr<EAPOD>>& m_pod_interfaces;
    const std::vector<int>&              m_type_map;        // exaStamp 0-indexed -> POD 1-indexed
    const size_t * const __restrict__    m_cell_particle_offset = nullptr;
    double * const __restrict__          m_descriptors = nullptr;
    // compute_derivative only: [(m+Mdesc*k)*3+xyz] -> per-particle aggregate array, or nullptr to
    // skip the derivative scatter entirely.
    double * const * const __restrict__  m_deriv_agg_ptrs = nullptr;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () (int jnum, ComputeBufferT& buf, int type, CellParticlesT cells) const
    {
      EAPOD& pod = *m_pod_interfaces[omp_get_thread_num()];

      pod.peratombase_descriptors_soa(buf.drx, buf.dry, buf.drz, buf.ext.type, jnum, m_type_map.data());

      const int Mdesc = pod.Mdesc;
      const int nClusters = pod.nClusters;
      const size_t p = m_cell_particle_offset[buf.cell] + buf.part;
      double * const __restrict__ out = m_descriptors + static_cast<size_t>(Mdesc) * nClusters * p;

      if (nClusters > 1)
      {
        const int ti0 = m_type_map[type] - 1;
        pod.peratomenvironment_descriptors(pod.pd, pod.pdd, pod.bd, pod.bdd, pod.tmpmem, ti0, jnum);
        for (int k = 0; k < nClusters; k++)
          for (int m = 0; m < Mdesc; m++)
            out[m + Mdesc*k] = pod.pd[k] * pod.bd[m];
      }
      else
      {
        for (int m = 0; m < Mdesc; m++) out[m] = pod.bd[m];
      }

      // Raw per-neighbor-pair Jacobian bdd[xyz+3*jj+3*jnum*m] = d(bd[m])/d(rij[xyz]), rij =
      // r_neighbor - r_central (see EAPOD::peratombase_descriptors / peratomenergyforce's DGEMV
      // consumption of bdd). Reduce it into a compact per-atom aggregate using the exact same
      // central+=/neighbor-=  scatter LAMMPS's own compute_podd_atom.cpp uses per local-atom row,
      // just summed straight into each atom's own slot instead of kept as a dense global matrix.
      // nClusters>1: apply the same product rule as PodGlobalOp (out[m,k]=pd[k]*bd[m], so its
      // derivative is bdd[m]*pd[k] + bd[m]*pdd[k]) -- pod.pd/pod.pdd were already populated above
      // by peratomenvironment_descriptors.
      if (m_deriv_agg_ptrs != nullptr)
      {
        for (int jj = 0; jj < jnum; jj++)
        {
          size_t nbh_cell=0, nbh_part=0;
          buf.nbh.get(jj, nbh_cell, nbh_part);
          const size_t nbh_p = m_cell_particle_offset[nbh_cell] + nbh_part;
          for (int m = 0; m < Mdesc; m++)
          {
            const int base = 3*jj + 3*jnum*m;
            const int nk_loop = (nClusters > 1) ? nClusters : 1;
            for (int k = 0; k < nk_loop; k++)
            {
              double vx, vy, vz;
              if (nClusters > 1)
              {
                const int nk = 3*jj + 3*jnum*k;
                vx = pod.bdd[0+base]*pod.pd[k] + pod.bd[m]*pod.pdd[0+nk];
                vy = pod.bdd[1+base]*pod.pd[k] + pod.bd[m]*pod.pdd[1+nk];
                vz = pod.bdd[2+base]*pod.pd[k] + pod.bd[m]*pod.pdd[2+nk];
              }
              else { vx = pod.bdd[0+base]; vy = pod.bdd[1+base]; vz = pod.bdd[2+base]; }

              const int mk3 = (m + Mdesc*k)*3;
              atomic_add_contribution(m_deriv_agg_ptrs[mk3+0][p],  vx);
              atomic_add_contribution(m_deriv_agg_ptrs[mk3+1][p],  vy);
              atomic_add_contribution(m_deriv_agg_ptrs[mk3+2][p],  vz);
              atomic_add_contribution(m_deriv_agg_ptrs[mk3+0][nbh_p], -vx);
              atomic_add_contribution(m_deriv_agg_ptrs[mk3+1][nbh_p], -vy);
              atomic_add_contribution(m_deriv_agg_ptrs[mk3+2][nbh_p], -vz);
            }
          }
        }
      }
    }
  };

}
