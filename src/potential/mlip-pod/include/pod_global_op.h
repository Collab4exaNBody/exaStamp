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

  // Global-array analogue of LAMMPS's compute pod/global (ML-POD/compute_pod_global.cpp): row 0 is
  // the system-wide per-element summed descriptor vector (incl. the nl1 one-body atom-count term),
  // rows 1..3*natoms are the gradient of that same row w.r.t. every atom's x/y/z -- indexed by atom
  // id-1 (not internal cell/particle order), which both matches LAMMPS's own row layout for direct
  // comparison and makes ghost/real folding automatic: a ghost neighbor carries the same field::id
  // as its real counterpart, so its contribution lands on the correct row with no separate
  // update_opt_from_ghost step. Single MPI rank only (see compute_descriptor_pod_global.cu) -- the
  // real compute_pod_global.cpp itself refuses nprocs>1, so this is an honest match, not a
  // limitation introduced here.
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PodGlobalOp
  {
    std::vector<std::shared_ptr<EAPOD>>& m_pod_interfaces;
    const std::vector<int>&              m_type_map;        // exaStamp 0-indexed -> POD 1-indexed
    const long                           m_ncols = 0;        // nCoeffAll
    double * const __restrict__          m_array = nullptr;  // row-major, (1+3*natoms) x m_ncols

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () (int jnum, ComputeBufferT& buf, int type, CellParticlesT cells) const
    {
      EAPOD& pod = *m_pod_interfaces[omp_get_thread_num()];
      const int ti0 = m_type_map[type] - 1;
      const int Mdesc = pod.Mdesc, nClusters = pod.nClusters;
      const int nCoeffPerElement = pod.nCoeffPerElement, nl1 = pod.nl1;
      const long ncols = m_ncols;
      // exaStamp's field::id is 0-indexed (unlike LAMMPS's 1-indexed atom->tag, hence no "-1" here) --
      // assumes ids are contiguous 0..natoms-1, matching how compute_descriptor_pod_global.cu sizes
      // the array (rows = 1+3*natoms).
      const uint64_t central_id = cells[buf.cell][field::id][buf.part];
      const long central_row = 1 + 3*static_cast<long>(central_id);

      if (nl1 > 0) atomic_add_contribution(m_array[nCoeffPerElement*ti0], 1.0);
      if (jnum == 0) return;

      pod.peratombase_descriptors_soa(buf.drx, buf.dry, buf.drz, buf.ext.type, jnum, m_type_map.data());
      if (nClusters > 1) pod.peratomenvironment_descriptors(pod.pd, pod.pdd, pod.bd, pod.bdd, pod.tmpmem, ti0, jnum);

      for (int m = 0; m < Mdesc; m++)
      {
        const int nj_loop = (nClusters > 1) ? nClusters : 1;
        for (int j = 0; j < nj_loop; j++)
        {
          const long k = (nClusters > 1)
                        ? nCoeffPerElement*ti0 + nl1 + m + j*Mdesc
                        : nCoeffPerElement*ti0 + nl1 + m;
          atomic_add_contribution(m_array[k], (nClusters > 1) ? pod.pd[j]*pod.bd[m] : pod.bd[m]);

          for (int jj = 0; jj < jnum; jj++)
          {
            size_t nbh_cell=0, nbh_part=0;
            buf.nbh.get(jj, nbh_cell, nbh_part);
            const uint64_t nbh_id = cells[nbh_cell][field::id][nbh_part];
            const long nbh_row = 1 + 3*static_cast<long>(nbh_id);
            const int nm = 3*jj + 3*jnum*m;

            double vx, vy, vz;
            if (nClusters > 1)
            {
              const int nk = 3*jj + 3*jnum*j;
              vx = pod.bdd[0+nm]*pod.pd[j] + pod.bd[m]*pod.pdd[0+nk];
              vy = pod.bdd[1+nm]*pod.pd[j] + pod.bd[m]*pod.pdd[1+nk];
              vz = pod.bdd[2+nm]*pod.pd[j] + pod.bd[m]*pod.pdd[2+nk];
            }
            else { vx = pod.bdd[0+nm]; vy = pod.bdd[1+nm]; vz = pod.bdd[2+nm]; }

            atomic_add_contribution(m_array[(central_row+0)*ncols+k],  vx);
            atomic_add_contribution(m_array[(central_row+1)*ncols+k],  vy);
            atomic_add_contribution(m_array[(central_row+2)*ncols+k],  vz);
            atomic_add_contribution(m_array[(nbh_row+0)*ncols+k],     -vx);
            atomic_add_contribution(m_array[(nbh_row+1)*ncols+k],     -vy);
            atomic_add_contribution(m_array[(nbh_row+2)*ncols+k],     -vz);
          }
        }
      }
    }
  };
}
