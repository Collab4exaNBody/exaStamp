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

// Stillinger-Weber force functor, one thread per central atom i, full neighbor list.
//  - two-body : i takes the full pair force and half the pair energy (j does the same on its side),
//               nothing is written to j.
//  - three-body : every j<k pair of i's short list; forces on j and k are atomically scattered
//               (ghost j/k included, they are folded back by update_force_energy_from_ghost),
//               the triplet energy goes to the central atom i.
// Per-atom virial follows the exaStamp convention sum_a tensor(F_a, r_a) split per central atom.

#include <cmath>
#include <type_traits>
#include <exanb/compute/compute_pair_buffer.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <exanb/compute/compute_pair_traits.h>
#include <exanb/core/concurent_add_contributions.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/cuda/cuda.h>

#include "sw_params.h"

#ifndef XSTAMP_SW_MAX_NEIGHBORS
#define XSTAMP_SW_MAX_NEIGHBORS 64
#endif

namespace exaStamp
{
  using namespace exanb;

  static constexpr size_t SW_MAX_NEIGHBORS = XSTAMP_SW_MAX_NEIGHBORS;

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) SWComputeBuffer
  {
    alignas(onika::memory::DEFAULT_ALIGNMENT) int type[SW_MAX_NEIGHBORS];
  };

  // Stores neighbor type alongside the default buffer fields. Never writes past capacity :
  // on overflow count keeps growing so that the functor can detect it and abort.
  struct SWCopyParticleType
  {
    template <class ComputeBufferT, class FieldArraysT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator()(ComputeBufferT &tab, const Vec3d &dr, double d2, FieldArraysT cells, size_t cell_b, size_t p_b, double weight) const noexcept
    {
      if( size_t(tab.count) < ComputeBufferT::MaxNeighbors )
      {
        tab.ext.type[tab.count] = cells[cell_b][field::type][p_b];
        exanb::DefaultComputePairBufferAppendFunc{}(tab, tab.count, dr, d2, cells, cell_b, p_b, weight);
      }
      ++ tab.count;
    }
  };

  struct SWForceOp
  {
    const SWTriplet * m_triplets = nullptr;
    int m_ntypes = 0;
    bool m_threebody = true;

    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE const SWTriplet & param(int i, int j, int k) const
    {
      return m_triplets[ ( i * m_ntypes + j ) * m_ntypes + k ];
    }

    // pair term : fforce such that f_i = (r_i - r_j) * fforce, eng = phi2(r)
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static void twobody(const SWTriplet & p, double rsq, double & fforce, double & eng)
    {
      const double r = std::sqrt(rsq);
      const double rinvsq = 1.0 / rsq;
      const double rp = ( p.powerp == 4.0 ) ? rinvsq * rinvsq : std::pow( r , -p.powerp );
      const double rq = ( p.powerq == 0.0 ) ? 1.0 : std::pow( r , -p.powerq );
      const double rainv = 1.0 / ( r - p.cut );
      const double rainvsq = rainv * rainv * r;
      const double expsrainv = std::exp( p.sigma * rainv );
      fforce = ( p.c1 * rp - p.c2 * rq + ( p.c3 * rp - p.c4 * rq ) * rainvsq ) * expsrainv * rinvsq;
      eng = ( p.c5 * rp - p.c6 * rq ) * expsrainv;
    }

    // triplet term j-i-k, delr1 = r_j - r_i, delr2 = r_k - r_i ; returns forces on j and k, f_i = -(fj+fk)
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static void threebody(const SWTriplet & pij, const SWTriplet & pik, const SWTriplet & pijk,
                                                                    double rsq1, double rsq2, const Vec3d & delr1, const Vec3d & delr2,
                                                                    Vec3d & fj, Vec3d & fk, double & eng)
    {
      const double r1 = std::sqrt(rsq1);
      const double rinvsq1 = 1.0 / rsq1;
      const double rainv1 = 1.0 / ( r1 - pij.cut );
      const double gsrainv1 = pij.sigma_gamma * rainv1;
      const double gsrainvsq1 = gsrainv1 * rainv1 / r1;
      const double expgsrainv1 = std::exp(gsrainv1);

      const double r2 = std::sqrt(rsq2);
      const double rinvsq2 = 1.0 / rsq2;
      const double rainv2 = 1.0 / ( r2 - pik.cut );
      const double gsrainv2 = pik.sigma_gamma * rainv2;
      const double gsrainvsq2 = gsrainv2 * rainv2 / r2;
      const double expgsrainv2 = std::exp(gsrainv2);

      const double rinv12 = 1.0 / ( r1 * r2 );
      const double cs = dot(delr1,delr2) * rinv12;
      const double delcs = cs - pijk.costheta;
      const double delcssq = delcs * delcs;

      const double facexp = expgsrainv1 * expgsrainv2;
      const double facrad = pijk.lambda_epsilon * facexp * delcssq;
      const double frad1 = facrad * gsrainvsq1;
      const double frad2 = facrad * gsrainvsq2;
      const double facang = pijk.lambda_epsilon2 * facexp * delcs;
      const double facang12 = rinv12 * facang;
      const double csfacang = cs * facang;
      const double csfac1 = rinvsq1 * csfacang;
      const double csfac2 = rinvsq2 * csfacang;

      fj = delr1 * ( frad1 + csfac1 ) - delr2 * facang12;
      fk = delr2 * ( frad2 + csfac2 ) - delr1 * facang12;
      eng = facrad;
    }

    template <class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3d &virial, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    // All force updates are true atomics (host and device), locks are not needed.
    template <class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3dT &virial, CellParticlesT cells, GridCellLocksT, ParticleLockT &) const
    {
      static constexpr bool compute_virial = std::is_same_v<Mat3dT, Mat3d>;

      if( n > ComputeBufferT::MaxNeighbors )
      {
        printf("sw_force: %d neighbors within SW cutoff exceed buffer capacity (%d), rebuild with larger XSTAMP_SW_MAX_NEIGHBORS\n", int(n), int(ComputeBufferT::MaxNeighbors));
        ONIKA_CU_ABORT();
      }

      const int itype = type;

      // short list : keep neighbors inside their pair cutoff (i,j,j)
      int ns = 0;
      for(size_t jj = 0; jj < n; ++jj)
      {
        const int jtype = buf.ext.type[jj];
        if( buf.d2[jj] < param(itype,jtype,jtype).cutsq )
        {
          if( size_t(ns) != jj ) { buf.copy(jj, ns); buf.ext.type[ns] = jtype; }
          ++ ns;
        }
      }

      double _en = 0.0;
      Vec3d fi = { 0.0, 0.0, 0.0 };
      Mat3d vir = { 0.,0.,0., 0.,0.,0., 0.,0.,0. };

      // two-body, full list : whole force and half energy on i only
      for(int jj = 0; jj < ns; ++jj)
      {
        const SWTriplet & pij = param(itype, buf.ext.type[jj], buf.ext.type[jj]);
        const Vec3d dr = { buf.drx[jj], buf.dry[jj], buf.drz[jj] }; // r_j - r_i
        double fpair = 0.0, eng = 0.0;
        twobody(pij, buf.d2[jj], fpair, eng);
        fi = fi - dr * fpair;
        _en += 0.5 * eng;
        if constexpr ( compute_virial ) { vir = vir + tensor(dr, dr) * ( 0.5 * fpair ); }
      }

      // three-body, j<k pairs of the short list
      if( m_threebody )
      {
        for(int jj = 0; jj < ns - 1; ++jj)
        {
          const int jtype = buf.ext.type[jj];
          const SWTriplet & pij = param(itype, jtype, jtype);
          const Vec3d delr1 = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
          const double rsq1 = buf.d2[jj];
          Vec3d fjsum = { 0.0, 0.0, 0.0 };

          for(int kk = jj + 1; kk < ns; ++kk)
          {
            const int ktype = buf.ext.type[kk];
            const Vec3d delr2 = { buf.drx[kk], buf.dry[kk], buf.drz[kk] };
            Vec3d fj, fk;
            double eng = 0.0;
            threebody(pij, param(itype,ktype,ktype), param(itype,jtype,ktype), rsq1, buf.d2[kk], delr1, delr2, fj, fk, eng);

            fi = fi - fj - fk;
            fjsum = fjsum + fj;
            _en += eng;
            if constexpr ( compute_virial ) { vir = vir + tensor(fj, delr1) + tensor(fk, delr2); }

            size_t cell_k = 0, p_k = 0;
            buf.nbh.get(kk, cell_k, p_k);
            atomic_add_contribution(cells[cell_k][field::fx][p_k], fk.x);
            atomic_add_contribution(cells[cell_k][field::fy][p_k], fk.y);
            atomic_add_contribution(cells[cell_k][field::fz][p_k], fk.z);
          }

          size_t cell_j = 0, p_j = 0;
          buf.nbh.get(jj, cell_j, p_j);
          atomic_add_contribution(cells[cell_j][field::fx][p_j], fjsum.x);
          atomic_add_contribution(cells[cell_j][field::fy][p_j], fjsum.y);
          atomic_add_contribution(cells[cell_j][field::fz][p_j], fjsum.z);
        }
      }

      // central atom force may be concurrently updated by other atoms' triplets
      atomic_add_contribution(fx, fi.x);
      atomic_add_contribution(fy, fi.y);
      atomic_add_contribution(fz, fi.z);
      en += _en;
      if constexpr ( compute_virial ) { virial += vir; }
    }
  };

}

namespace exanb
{
  template <> struct ComputePairTraits<exaStamp::SWForceOp>
  {
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = true;
    static inline constexpr bool HasParticleContextStart = false;
    static inline constexpr bool HasParticleContext = false;
    static inline constexpr bool HasParticleContextStop = false;
    static inline constexpr bool RequiresNbhOptionalData = true;
  };
}
