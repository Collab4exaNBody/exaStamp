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

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>
#include <exanb/core/concurent_add_contributions.h>

#include "eapod.h"
#include "pod_config.h"

namespace exaStamp
{

  using namespace exanb;

  // Scratch arrays (rij, fij, ti, tj) are owned by each EAPOD instance and
  // allocated in allocate_temp_memory — the compute buffer stays small.
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PodComputeBuffer
  {
    alignas(onika::memory::DEFAULT_ALIGNMENT) int type[exanb::MAX_PARTICLE_NEIGHBORS];
  };

  struct CopyParticleType
  {
    template<class ComputeBufferT, class FieldArraysT>
    inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, FieldArraysT cells, size_t cell_b, size_t p_b, double weight) const noexcept
    {
      assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
      tab.ext.type[tab.count] = cells[cell_b][field::type][p_b];
      exanb::DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, weight );
    }
  };

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PodForceOp
  {
    std::vector<std::shared_ptr<EAPOD>>& m_pod_interfaces;
    const std::vector<int>&              m_type_map;        // exaStamp 0-indexed → POD 1-indexed
    const bool conv_energy_units = true;
    static constexpr double conv_energy_factor = EXASTAMP_CONST_QUANTITY( 1. * eV );
    const bool eflag;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (size_t n, ComputeBufferT& buf, double& en, double& fx, double& fy, double& fz,
       int type, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (size_t n, ComputeBufferT& buf, double& en, double& fx, double& fy, double& fz,
       int type, Mat3d& virial, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (size_t n, ComputeBufferT& buf, double& en, double& fx, double& fy, double& fz,
       int type, CellParticlesT cells, GridCellLocksT locks, ParticleLockT& lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (int jnum, ComputeBufferT& buf, double& en, double& fx, double& fy, double& fz,
       int type, Mat3dT& virial, CellParticlesT cells, GridCellLocksT locks, ParticleLockT& lock_a) const
    {
      static constexpr bool vflag = std::is_same_v<Mat3dT, Mat3d>;

      Mat3dT _vir;
      double _en = 0.;
      double _fx = 0.;
      double _fy = 0.;
      double _fz = 0.;

      size_t tid = omp_get_thread_num();
      EAPOD& pod = *m_pod_interfaces[tid];

      // peratomenergyforce2_soa packs SoA→AoS using EAPOD-owned scratch,
      // applies type_map for both central and neighbour types, returns per-atom energy.
      const double e = pod.peratomenergyforce2_soa(
          buf.drx, buf.dry, buf.drz,
          type, buf.ext.type, jnum, m_type_map.data());

      for (int jj = 0; jj < jnum; ++jj) {
        const double fx_j = pod.soa_fij[jj*3+0] * conv_energy_factor;
        const double fy_j = pod.soa_fij[jj*3+1] * conv_energy_factor;
        const double fz_j = pod.soa_fij[jj*3+2] * conv_energy_factor;

        _fx += fx_j;
        _fy += fy_j;
        _fz += fz_j;

        size_t cell_b = 0, p_b = 0;
        buf.nbh.get(jj, cell_b, p_b);
        atomic_add_contribution(cells[cell_b][field::fx][p_b], -fx_j);
        atomic_add_contribution(cells[cell_b][field::fy][p_b], -fy_j);
        atomic_add_contribution(cells[cell_b][field::fz][p_b], -fz_j);

        //
        if constexpr (vflag)
        {
          const Mat3d virial_contrib = tensor(Vec3d{fx_j, fy_j, fz_j}, Vec3d{buf.drx[jj], buf.dry[jj], buf.drz[jj]});
          _vir -= virial_contrib;
        }

      }

      if constexpr (vflag)
        atomic_add_contribution(virial, _vir);
      if (eflag)
        atomic_add_contribution(en, e * conv_energy_factor);
      atomic_add_contribution(fx, _fx);
      atomic_add_contribution(fy, _fy);
      atomic_add_contribution(fz, _fz);
    }
  };

}
