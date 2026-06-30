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

namespace exaStamp
{

  using namespace exanb;

  // This is the additional information needed on central particle's neighbors
  // Here: we need the neighbors type 
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PaceComputeBuffer
  {
    alignas(onika::memory::DEFAULT_ALIGNMENT) int type[exanb::MAX_PARTICLE_NEIGHBORS];
  };

  // This is the functor that populates the extra array (type) from the grid type field
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

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PaceForceOp
  {
    std::vector<std::shared_ptr<ACEImpl>>& m_pace_interfaces;
    const bool conv_energy_units = true;
    static constexpr double conv_energy_factor = EXASTAMP_CONST_QUANTITY( 1. * eV );
    const bool eflag;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells
      ) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3d& virial,
      CellParticlesT cells
      ) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {
      FakeMat3d virial;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (
      int jnum,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3dT& virial,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {

      static constexpr bool vflag = std::is_same_v<Mat3dT, Mat3d>;
    
      Mat3dT _vir;
      double _en = 0.;
      double _fx = 0.;
      double _fy = 0.;
      double _fz = 0.;

      size_t tid = omp_get_thread_num();
      ACEImpl& acecalc = *m_pace_interfaces[tid];
      acecalc.ace->resize_neighbours_cache(jnum);

      acecalc.ace->compute_atom_drxyz(type, buf.drx, buf.dry, buf.drz, buf.ext.type, jnum);

      double fij[3];
      for (int jj = 0; jj < jnum; ++jj) {
        fij[0] = acecalc.ace->neighbours_forces(jj,0) * conv_energy_factor;
        fij[1] = acecalc.ace->neighbours_forces(jj,1) * conv_energy_factor;
        fij[2] = acecalc.ace->neighbours_forces(jj,2) * conv_energy_factor;

        // Always collect forces contribution
        _fx += fij[0];
        _fy += fij[1];
        _fz += fij[2];

        // Always send reciprocal force contribution to neighbors
        size_t cell_b=0, p_b=0;
        buf.nbh.get(jj, cell_b, p_b);        
        atomic_add_contribution( cells[cell_b][field::fx][p_b], -fij[0]);
        atomic_add_contribution( cells[cell_b][field::fy][p_b], -fij[1]);
        atomic_add_contribution( cells[cell_b][field::fz][p_b], -fij[2]);

        // Collect the virial cobtribution to the central atom only if needed
        if constexpr ( vflag ) 
        { 
          const Mat3d virial_contrib = tensor( Vec3d{ fij[0], fij[1], fij[2] }, Vec3d{ buf.drx[jj], buf.dry[jj], buf.drz[jj] } );
          _vir -= virial_contrib;
        }

      }

      // Collect the virial cobtribution to the central atom only if needed
      if constexpr ( vflag ) atomic_add_contribution( virial, _vir );
      // Collect the energy cobtribution to the central atom only if needed
      if ( eflag ) atomic_add_contribution( en, acecalc.ace->e_atom * conv_energy_factor);
      // Always collect forces contributions
      atomic_add_contribution( fx, _fx);
      atomic_add_contribution( fy, _fy);
      atomic_add_contribution( fz, _fz);
    }
  };

}
