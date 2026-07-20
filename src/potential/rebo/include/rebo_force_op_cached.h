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

// Optimized REBO variant. It keeps the LAMMPS-like active set for force
// derivatives, while avoiding repeated central-neighbour switch evaluations and
// filtering secondary pairs by squared cutoff before taking square roots.

#include <cmath>
#include <exaStamp/unit_system.h>
#include <exanb/core/concurent_add_contributions.h>
#include <onika/physics/constants.h>
#include <onika/physics/units.h>

#include "rebo_params.h"

namespace exaStamp
{
  using namespace exanb;

  static constexpr size_t REBO_MAX_NEIGHBORS = 256;
  static constexpr size_t REBO_MAX_REBO_NEIGHBORS = 32;

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) ReboComputeBuffer
  {
    alignas(onika::memory::DEFAULT_ALIGNMENT) int type[REBO_MAX_NEIGHBORS];
  };

  struct CopyParticleType
  {
    template <class ComputeBufferT, class FieldArraysT> ONIKA_HOST_DEVICE_FUNC inline void operator()(ComputeBufferT &tab, const Vec3d &dr, double d2, FieldArraysT cells, size_t cell_b, size_t p_b, double weight) const noexcept
    {
      assert(ssize_t(tab.count) < ssize_t(tab.MaxNeighbors));
      tab.ext.type[tab.count] = cells[cell_b][field::type][p_b];
      exanb::DefaultComputePairBufferAppendFunc{}(tab, dr, d2, cells, cell_b, p_b, weight);
    }
  };

  // ===========================================================================
  // REBO helper functions
  // ===========================================================================

  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_Sp(double Xij, double Xmin, double Xmax, double &dX)
  {
    double cutoff;

    const double t = (Xij - Xmin) / (Xmax - Xmin);
    if (t <= 0.0)
    {
      cutoff = 1.0;
      dX = 0.0;
    }
    else if (t >= 1.0)
    {
      cutoff = 0.0;
      dX = 0.0;
    }
    else
    {
      cutoff = 0.5 * (1.0 + cos(t * M_PI));
      dX = (-0.5 * M_PI * sin(t * M_PI)) / (Xmax - Xmin);
    }
    return cutoff;
  };

  // Cubic switching function — used for Tij planarity (thmin=-1, thmax=-0.995)
  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_Sp2(double Xij, double Xmin, double Xmax, double &dX)
  {
    double cutoff;

    const double t = (Xij - Xmin) / (Xmax - Xmin);
    if (t <= 0.0)
    {
      cutoff = 1.0;
      dX = 0.0;
    }
    else if (t >= 1.0)
    {
      cutoff = 0.0;
      dX = 0.0;
    }
    else
    {
      cutoff = (1.0 - (t * t * (3.0 - 2.0 * t)));
      dX = 6.0 * (t * t - t) / (Xmax - Xmin);
    }
    return cutoff;
  };

  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_Sp5th(double x, const double *c, double &df)
  {
    double f = c[5] * x;
    double d = 5.0 * c[5] * x;
    f += c[4];
    d += 4.0 * c[4];
    f *= x;
    d *= x;
    f += c[3];
    d += 3.0 * c[3];
    f *= x;
    d *= x;
    f += c[2];
    d += 2.0 * c[2];
    f *= x;
    d *= x;
    f += c[1];
    d += c[1];
    f *= x;
    f += c[0];
    df = d;
    return f;
  }

  // 2D bicubic: f(x,y) = sum_{i,j=0}^{3} c[i*4+j] * x^i * y^j
  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_Spbicubic(double x, double y, const double *coeffs, double df[2])
  {
    double f, xn, yn, xn1, yn1, c;
    int i, j;

    f = 0.0;
    df[0] = 0.0;
    df[1] = 0.0;

    xn = 1.0;
    for (i = 0; i < 4; i++)
    {
      yn = 1.0;
      for (j = 0; j < 4; j++)
      {
        c = coeffs[i * 4 + j];

        f += c * xn * yn;
        if (i > 0)
          df[0] += c * ((double)i) * xn1 * yn;
        if (j > 0)
          df[1] += c * ((double)j) * xn * yn1;

        yn1 = yn;
        yn *= y;
      }
      xn1 = xn;
      xn *= x;
    }

    return f;
  }

  // 3D tricubic: f(x,y,z) = sum_{i,j,k=0}^{3} c[16i+4j+k] * x^i * y^j * z^k
  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_Sptricubic(double x, double y, double z, const double *coeffs, double df[3])
  {
    double f, ir, jr, kr, xn, yn, zn, xn1, yn1, zn1, c;
    int i, j, k;

    f = 0.0;
    df[0] = 0.0;
    df[1] = 0.0;
    df[2] = 0.0;

    xn = 1.0;
    for (i = 0; i < 4; i++)
    {
      ir = (double)i;
      yn = 1.0;
      for (j = 0; j < 4; j++)
      {
        jr = (double)j;
        zn = 1.0;
        for (k = 0; k < 4; k++)
        {
          kr = (double)k;
          c = coeffs[16 * i + 4 * j + k];
          f += c * xn * yn * zn;
          if (i > 0)
            df[0] += c * ir * xn1 * yn * zn;
          if (j > 0)
            df[1] += c * jr * xn * yn1 * zn;
          if (k > 0)
            df[2] += c * kr * xn * yn * zn1;
          zn1 = zn;
          zn *= z;
        }
        yn1 = yn;
        yn *= y;
      }
      xn1 = xn;
      xn *= x;
    }

    return f;
  }

  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_gSpline(double costh, double N, int typei, const ReboParamsRO &p, double &dgdc, double &dgdN)
  {
    dgdc = dgdN = 0.0;
    if (typei == 0)
    {
      if (costh < p.gCdom[0])
        costh = p.gCdom[0];
      if (costh > p.gCdom[4])
        costh = p.gCdom[4];
      int iv = 3;
      for (int i = 0; i < 4; ++i)
        if (costh >= p.gCdom[i] && costh <= p.gCdom[i + 1])
        {
          iv = i;
          break;
        }
      double dg1, dg2;
      const double g1 = rebo_Sp5th(costh, p.gC1 + iv * 6, dg1);
      const double g2 = rebo_Sp5th(costh, p.gC2 + iv * 6, dg2);
      if (N >= p.NCmax)
      {
        dgdc = dg2;
        return g2;
      }
      else if (N <= p.NCmin)
      {
        dgdc = dg1;
        return g1;
      }
      else
      {
        double dcut;
        const double cut = rebo_Sp(N, p.NCmin, p.NCmax, dcut);
        dgdc = dg2 + cut * (dg1 - dg2);
        dgdN = dcut * (g1 - g2);
        return g2 + cut * (g1 - g2);
      }
    }
    else
    {
      if (costh < p.gHdom[0])
        costh = p.gHdom[0];
      if (costh > p.gHdom[3])
        costh = p.gHdom[3];
      int iv = 2;
      for (int i = 0; i < 3; ++i)
        if (costh >= p.gHdom[i] && costh <= p.gHdom[i + 1])
        {
          iv = i;
          break;
        }
      return rebo_Sp5th(costh, p.gH + iv * 6, dgdc);
    }
  }

  // P_ij bicubic correction (zero for H central atom)
  // dom layout (flat): [lo_nC, hi_nC, lo_nH, hi_nH]
  // coeff layout: [nC_patch][nH_patch][16], flat = ((x*ny)+y)*16
  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_PijSpline(double NijC, double NijH, int typei, int typej, const ReboParamsRO &p, double dN2[2])
  {
    dN2[0] = dN2[1] = 0.0;
    if (typei == 1)
      return 0.0;
    const double *dom = (typej == 0) ? p.pCCdom : p.pCHdom;
    const double *tbl = (typej == 0) ? p.pCC : p.pCH;
    if (NijC < dom[0])
      NijC = dom[0];
    if (NijC > dom[1])
      NijC = dom[1];
    if (NijH < dom[2])
      NijH = dom[2];
    if (NijH > dom[3])
      NijH = dom[3];
    int x = (int)std::floor(NijC);
    int y = (int)std::floor(NijH);
    if (NijC == dom[1])
      --x;
    if (NijH == dom[3])
      --y;
    const int ny = (int)dom[3];
    return rebo_Spbicubic(NijC, NijH, tbl + (x * ny + y) * 16, dN2);
  }

  // piRC tricubic correction
  // dom layout (flat): [lo0,hi0, lo1,hi1, lo2,hi2]
  // coeff layout: [x][y][z][64], flat = ((x*ny+y)*nz+z)*64
  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_piRCSpline(double Nij, double Nji, double Nijconj, int typei, int typej, const ReboParamsRO &p, double dN3[3])
  {
    dN3[0] = dN3[1] = dN3[2] = 0.0;
    const double *dom;
    const double *tbl;
    if (typei == 0 && typej == 0)
    {
      dom = p.piCCdom;
      tbl = p.piCC;
    }
    else if (typei == 1 && typej == 1)
    {
      dom = p.piHHdom;
      tbl = p.piHH;
    }
    else
    {
      dom = p.piCHdom;
      tbl = p.piCH;
    }
    if (Nij < dom[0])
      Nij = dom[0];
    if (Nij > dom[1])
      Nij = dom[1];
    if (Nji < dom[2])
      Nji = dom[2];
    if (Nji > dom[3])
      Nji = dom[3];
    if (Nijconj < dom[4])
      Nijconj = dom[4];
    if (Nijconj > dom[5])
      Nijconj = dom[5];
    int x = (int)std::floor(Nij), y = (int)std::floor(Nji), z = (int)std::floor(Nijconj);
    if (Nij == dom[1])
      --x;
    if (Nji == dom[3])
      --y;
    if (Nijconj == dom[5])
      --z;
    const int ny = (int)dom[3], nz = (int)dom[5];
    return rebo_Sptricubic(Nij, Nji, Nijconj, tbl + ((x * ny + y) * nz + z) * 64, dN3);
  }

  ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE static double rebo_TijSpline(double Nij, double Nji, double Nijconj, const ReboParamsRO &p, double dN3[3])
  {
    dN3[0] = dN3[1] = dN3[2] = 0.0;
    const double *dom = p.Tijdom;
    if (Nij < dom[0])
      Nij = dom[0];
    if (Nij > dom[1])
      Nij = dom[1];
    if (Nji < dom[2])
      Nji = dom[2];
    if (Nji > dom[3])
      Nji = dom[3];
    if (Nijconj < dom[4])
      Nijconj = dom[4];
    if (Nijconj > dom[5])
      Nijconj = dom[5];
    int x = (int)std::floor(Nij), y = (int)std::floor(Nji), z = (int)std::floor(Nijconj);
    if (Nij == dom[1])
      --x;
    if (Nji == dom[3])
      --y;
    if (Nijconj == dom[5])
      --z;
    const int ny = (int)dom[3], nz = (int)dom[5];
    return rebo_Sptricubic(Nij, Nji, Nijconj, p.Tijc + ((x * ny + y) * nz + z) * 64, dN3);
  }

  // ===========================================================================
  // REBOPairwiseForceOp
  // Pure two-body term: the Brenner repulsive potential VR depends only on
  // r_ij and the two atom types — no coordination numbers, no three/four-body
  // sums. Split out from REBOManyBodyForceOp (below) so the cheap i-j-only
  // part and the expensive many-body part run as two separate kernels and can
  // be profiled independently.
  //
  // del = r_i - r_j  (LAMMPS convention, used internally)
  // ===========================================================================

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) REBOPairwiseForceOp
  {
    const ReboParamsRO *m_params = nullptr;

    template <class ComputeBufferT, class CellParticlesT> ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT> ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3d &virial, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT> ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT> ONIKA_HOST_DEVICE_FUNC inline void operator()(int jnum, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3dT &virial, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      static constexpr double TOL = 1.0e-9;
      static constexpr double conv = EXASTAMP_CONST_QUANTITY(1. * eV);

      double _en = 0.0, _fx = 0.0, _fy = 0.0, _fz = 0.0;
      const int itype = type;
      const ReboParamsRO &p = *m_params;

      for (int jj = 0; jj < jnum; ++jj)
      {
        const int jtype = buf.ext.type[jj];
        const double rsq = buf.d2[jj];
        if (rsq >= p.rcmaxsq[itype][jtype])
          continue;
        const double rij = std::sqrt(rsq);
        if (rij <= TOL)
          continue;

        double dwij;
        const double wij = rebo_Sp(rij, p.rcmin[itype][jtype], p.rcmax[itype][jtype], dwij);
        if (wij <= TOL)
          continue;

        // del = r_i - r_j  (LAMMPS convention)
        const double delx = -buf.drx[jj];
        const double dely = -buf.dry[jj];
        const double delz = -buf.drz[jj];

        const double Qij = p.Q[itype][jtype];
        const double Aij = p.A[itype][jtype];
        const double alph = p.alpha[itype][jtype];
        const double ea = std::exp(-alph * rij);
        double VR = wij * (1.0 + Qij / rij) * Aij * ea;
        const double pre = wij * Aij * ea;
        double dVR = pre * (-alph - Qij / rsq - Qij * alph / rij) + VR / wij * dwij;
        VR *= conv;
        dVR *= conv;

        _en += 0.5 * VR;

        // --- Direct pair force (×0.5 for full-list) ---
        const double fpair = -dVR / rij;
        _fx += 0.5 * delx * fpair;
        _fy += 0.5 * dely * fpair;
        _fz += 0.5 * delz * fpair;

        size_t cell_j, p_j;
        buf.nbh.get(jj, cell_j, p_j);
        atomic_add_contribution(cells[cell_j][field::fx][p_j], -0.5 * delx * fpair);
        atomic_add_contribution(cells[cell_j][field::fy][p_j], -0.5 * dely * fpair);
        atomic_add_contribution(cells[cell_j][field::fz][p_j], -0.5 * delz * fpair);
      }

      atomic_add_contribution(en, _en);
      atomic_add_contribution(fx, _fx);
      atomic_add_contribution(fy, _fy);
      atomic_add_contribution(fz, _fz);
    }
  };

  // ===========================================================================
  // REBOManyBodyForceOp
  // Attractive term VA scaled by the many-body bond order bij = 0.5*(pij+pji)
  // + piRC + Tij*Etmp_Tij. This is the expensive part: three-body sigma-pi sums
  // (k/l loops) and the four-body Tij dihedral double loop (k×l), plus their
  // higher-order coordination-number corrections (mm sub-loops). Split out
  // from the pure pairwise VR term (REBOPairwiseForceOp, above) so the two can
  // be profiled and launched as separate kernels.
  //
  // Full-list implementation: every pair (i,j) is visited twice (once as i's
  // neighbour, once as j's neighbour).  ALL energy and force contributions are
  // scaled by 0.5 so that the two visits sum to the correct single-pair value.
  //
  // buf.drx[jj] = r_j - r_i  (exaStamp convention, "from i to j")
  // del = r_i - r_j           (LAMMPS convention, used internally)
  //
  // All of j's REBO neighbours are accessible through i's buffer because the
  // REBO cutoff (~2 Å) << rcut_max (~10 Å), so r_{j,l} < 2 Å and
  // r_{i,l} < r_{i,j} + r_{j,l} < 4 Å < rcut_max.
  // ===========================================================================

  template <class NijcFieldT, class NijhFieldT, class NconjFieldT> struct alignas(onika::memory::DEFAULT_ALIGNMENT) REBOManyBodyForceOp
  {
    const ReboParamsRO *m_params = nullptr;
    NijcFieldT m_nijc_field = {};   // read neighbour's nC
    NijhFieldT m_nijh_field = {};   // read neighbour's nH
    NconjFieldT m_nconj_field = {}; // read neighbour's Nconj

    template <class ComputeBufferT, class CellParticlesT> ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, nC_central, nH_central, Nconj_central, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT> ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3d &virial, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, nC_central, nH_central, Nconj_central, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT> ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, nC_central, nH_central, Nconj_central, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT> ONIKA_HOST_DEVICE_FUNC inline void operator()(int jnum, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3dT &virial, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      // std::min/std::max are constexpr host functions, not device-callable
      // without --expt-relaxed-constexpr; onika::cuda::min/max are the
      // device-safe drop-in equivalents already used elsewhere (e.g. EAM).
      using onika::cuda::max;
      using onika::cuda::min;

      // static constexpr bool   vflag = std::is_same_v<Mat3dT, Mat3d>;
      static constexpr double TOL = 1.0e-9;
      static constexpr double thmin = -1.0;
      static constexpr double thmax = -0.995;
      // CH.Rebo parameters are in LAMMPS metal units (eV, Å); convert to exaStamp internal units
      static constexpr double conv = EXASTAMP_CONST_QUANTITY(1. * eV);

      Mat3dT _vir;
      double _en = 0.0, _fx = 0.0, _fy = 0.0, _fz = 0.0;
      const int itype = type;
      const ReboParamsRO &p = *m_params;

      // -----------------------------------------------------------------------
      // Pass 1 — REBO coordination numbers NiC, NiH for atom i (precomputed fields)
      // -----------------------------------------------------------------------
      const double NiC = nC_central;
      const double NiH = nH_central;

      // -----------------------------------------------------------------------
      // Cache central-neighbour geometry once. The buffer radius is intentionally
      // much larger than the REBO bond radius (3*rcmax), but most i-side loops
      // only need true REBO neighbours of the central atom.
      // -----------------------------------------------------------------------
      int rebo_idx[REBO_MAX_REBO_NEIGHBORS];
      double rebo_r[REBO_MAX_REBO_NEIGHBORS];
      double rebo_w[REBO_MAX_REBO_NEIGHBORS];
      double rebo_dw[REBO_MAX_REBO_NEIGHBORS];
      int rebo_count = 0;

      for (int aa = 0; aa < jnum; ++aa)
      {
        const int atype = buf.ext.type[aa];
        const double rsq = buf.d2[aa];
        const double r = std::sqrt(rsq);
        double dw;
        const double w = rebo_Sp(r, p.rcmin[itype][atype], p.rcmax[itype][atype], dw);
        if (rsq < p.rcmaxsq[itype][atype] && r > TOL)
        {
          if (rebo_count >= int(REBO_MAX_REBO_NEIGHBORS))
          {
            assert(false);
            continue;
          }
          rebo_idx[rebo_count] = aa;
          rebo_r[rebo_count] = r;
          rebo_w[rebo_count] = w;
          rebo_dw[rebo_count] = dw;
          ++rebo_count;
        }
      }

      // -----------------------------------------------------------------------
      // Cache (i,k) data that does NOT depend on which j is currently being
      // visited: k's cell/particle handle, and — for carbon k — its
      // coordination excluding the i-k bond (Nki) and the Sp cutoff evaluated
      // there. The piRC and Tij "third order" coordination corrections used
      // to recompute this once per (j,k) pair (i.e. once per k for every j),
      // even though it only depends on k. Computed once per k instead.
      // -----------------------------------------------------------------------
      size_t rebo_cell_k[REBO_MAX_REBO_NEIGHBORS];
      size_t rebo_p_k[REBO_MAX_REBO_NEIGHBORS];
      double rebo_SpNk[REBO_MAX_REBO_NEIGHBORS];
      double rebo_dNk[REBO_MAX_REBO_NEIGHBORS];

      for (int kpos = 0; kpos < rebo_count; ++kpos)
      {
        const int kk = rebo_idx[kpos];
        buf.nbh.get(kk, rebo_cell_k[kpos], rebo_p_k[kpos]);
        if (buf.ext.type[kk] == 0)
        {
          const double Nk = cells[rebo_cell_k[kpos]][m_nijc_field][rebo_p_k[kpos]] + cells[rebo_cell_k[kpos]][m_nijh_field][rebo_p_k[kpos]] - rebo_w[kpos];
          rebo_SpNk[kpos] = rebo_Sp(Nk, p.Nmin, p.Nmax, rebo_dNk[kpos]);
        }
        else
        {
          rebo_SpNk[kpos] = 0.0;
          rebo_dNk[kpos] = 0.0;
        }
      }

      // -----------------------------------------------------------------------
      // Pass 2 — pair loop over j
      // -----------------------------------------------------------------------
      for (int jpos = 0; jpos < rebo_count; ++jpos)
      {
        const int jj = rebo_idx[jpos];
        const int jtype = buf.ext.type[jj];
        const double rsq = buf.d2[jj];
        const double rij = rebo_r[jpos];

        // del = r_i - r_j  (LAMMPS convention)
        const double delx = -buf.drx[jj];
        const double dely = -buf.dry[jj];
        const double delz = -buf.drz[jj];

        double dwij;
        const double wij = rebo_w[jpos];
        dwij = rebo_dw[jpos];
        if (wij <= TOL)
          continue;

        // --- VA (attractive, 3 terms). VR (repulsive) is handled by REBOPairwiseForceOp. ---
        double VA = 0.0, dVA = 0.0;
        for (int m = 0; m < 3; ++m)
        {
          const double B = p.BIJc[itype][jtype][m];
          const double beta = p.Beta[itype][jtype][m];
          const double term = -wij * B * std::exp(-beta * rij);
          VA += term;
          dVA += -beta * term;
        }
        dVA += VA / wij * dwij;
        // Convert eV → exaStamp internal units; all force/energy terms derive from VA
        VA *= conv;
        dVA *= conv;

        // --- Coordination of i excluding j ---
        const double NijC = NiC - wij * (jtype == 0 ? 1.0 : 0.0);
        const double NijH = NiH - wij * (jtype == 1 ? 1.0 : 0.0);
        const double Nij = NijC + NijH;

        size_t cell_j, p_j;
        buf.nbh.get(jj, cell_j, p_j);
        const double nC_j = cells[cell_j][m_nijc_field][p_j];
        const double nH_j = cells[cell_j][m_nijh_field][p_j];
        const double Nconj_j = cells[cell_j][m_nconj_field][p_j];

        // ===================================================================
        // pij — i's sigma-pi bond order to j
        // NconjtmpI = Nconj[i] minus j's C contribution (O(1) with precomputed fields)
        // ===================================================================
        double unused_dSpNj;
        const double SpNj_I = rebo_Sp(nC_j + nH_j - wij, p.Nmin, p.Nmax, unused_dSpNj);
        const double NconjtmpI = Nconj_central - (jtype == 0 ? wij * SpNj_I : 0.0);

        // ===================================================================
        // pij — i's sigma-pi bond order to j: energy k-loop → compute pij → force k-loop
        // ===================================================================
        // g/dgdc/exp_lam per k, cached here so the force k-loop below doesn't
        // repeat the branchy rebo_gSpline eval and exp() call for the same
        // (j,k) pair (these depend on j via cosjik, so they can't be hoisted
        // above the j-loop like rebo_SpNk/rebo_dNk).
        double kj_g[REBO_MAX_REBO_NEIGHBORS], kj_dgdc[REBO_MAX_REBO_NEIGHBORS], kj_exp_lam[REBO_MAX_REBO_NEIGHBORS];
        double Etmp_pij = 0.0, tmp3_pij = 0.0;

        for (int kpos = 0; kpos < rebo_count; ++kpos)
        {
          const int kk = rebo_idx[kpos];
          if (kk == jj)
            continue;
          const int ktype = buf.ext.type[kk];
          const double rik = rebo_r[kpos];
          const double wik = rebo_w[kpos];
          // lambda (nonzero only when i is H)
          const double lam = (itype == 1) ? 4.0 * ((p.rho[ktype][1] - rik) - (p.rho[jtype][1] - rij)) : 0.0;
          const double exp_lam = std::exp(lam);
          kj_exp_lam[kpos] = exp_lam;

          // cos_jik: dot((ri-rj),(ri-rk))/(rij*rik); signs cancel with buf convention
          double cosjik = (buf.drx[jj] * buf.drx[kk] + buf.dry[jj] * buf.dry[kk] + buf.drz[jj] * buf.drz[kk]) / (rij * rik);
          cosjik = min(1.0, max(-1.0, cosjik));

          double dgdN;
          const double g = rebo_gSpline(cosjik, Nij, itype, p, kj_dgdc[kpos], dgdN);
          kj_g[kpos] = g;
          Etmp_pij += wik * g * exp_lam;
          tmp3_pij += wik * dgdN * exp_lam;
        }

        double dN2_pij[2];
        const double PijS = rebo_PijSpline(NijC, NijH, itype, jtype, p, dN2_pij);
        const double pij = 1.0 / std::sqrt(max(1.0 + Etmp_pij + PijS, TOL));
        const double tmp_pij = -0.5 * pij * pij * pij;

        // pij three-body forces (k-loop, ×0.5 for full-list)
        for (int kpos = 0; kpos < rebo_count; ++kpos)
        {
          const int kk = rebo_idx[kpos];
          if (kk == jj)
            continue;
          const int ktype = buf.ext.type[kk];
          const double rik = rebo_r[kpos];
          double dwik;
          const double wik = rebo_w[kpos];
          dwik = rebo_dw[kpos];
          const double exp_lam = kj_exp_lam[kpos];

          const double rikx = -buf.drx[kk], riky = -buf.dry[kk], rikz = -buf.drz[kk];

          double cosjik = (delx * rikx + dely * riky + delz * rikz) / (rij * rik);
          cosjik = min(1.0, max(-1.0, cosjik));

          const double g = kj_g[kpos];
          const double dgdc = kj_dgdc[kpos];

          const double inv_rij = 1.0 / rij, inv_rik = 1.0 / rik;
          const double inv_rij2 = inv_rij * inv_rij, inv_rik2 = inv_rik * inv_rik;
          const double irr = inv_rij * inv_rik;

          const double dcdrix = (delx + rikx) * irr - cosjik * (delx * inv_rij2 + rikx * inv_rik2);
          const double dcdriy = (dely + riky) * irr - cosjik * (dely * inv_rij2 + riky * inv_rik2);
          const double dcdriz = (delz + rikz) * irr - cosjik * (delz * inv_rij2 + rikz * inv_rik2);
          const double dcdrjx = -rikx * irr + cosjik * delx * inv_rij2;
          const double dcdrjy = -riky * irr + cosjik * dely * inv_rij2;
          const double dcdrjz = -rikz * irr + cosjik * delz * inv_rij2;
          const double dcdrkx = -delx * irr + cosjik * rikx * inv_rik2;
          const double dcdrky = -dely * irr + cosjik * riky * inv_rik2;
          const double dcdrkz = -delz * irr + cosjik * rikz * inv_rik2;

          // Scale = VA * 0.5 * tmp_pij * ... * 0.5(full-list) = VA * 0.25 * tmp_pij * ...
          const double pf_ang = VA * 0.25 * tmp_pij * wik * dgdc * exp_lam;
          double fix = -pf_ang * dcdrix, fiy = -pf_ang * dcdriy, fiz = -pf_ang * dcdriz;
          double fjx = -pf_ang * dcdrjx, fjy = -pf_ang * dcdrjy, fjz = -pf_ang * dcdrjz;
          double fkx = -pf_ang * dcdrkx, fky = -pf_ang * dcdrky, fkz = -pf_ang * dcdrkz;

          if (itype == 1)
          {
            const double pf_lam = VA * 0.25 * tmp_pij * wik * g * exp_lam * 4.0;
            fjx -= pf_lam * (-delx * inv_rij);
            fjy -= pf_lam * (-dely * inv_rij);
            fjz -= pf_lam * (-delz * inv_rij);
            fix -= pf_lam * ((-rikx * inv_rik) + (delx * inv_rij));
            fiy -= pf_lam * ((-riky * inv_rik) + (dely * inv_rij));
            fiz -= pf_lam * ((-rikz * inv_rik) + (delz * inv_rij));
            fkx -= pf_lam * (rikx * inv_rik);
            fky -= pf_lam * (riky * inv_rik);
            fkz -= pf_lam * (rikz * inv_rik);
          }

          const double pf_wik = VA * 0.25 * tmp_pij * dwik * g * exp_lam / rik;
          fix -= pf_wik * rikx;
          fiy -= pf_wik * riky;
          fiz -= pf_wik * rikz;
          fkx += pf_wik * rikx;
          fky += pf_wik * riky;
          fkz += pf_wik * rikz;

          const double pf_dN = VA * 0.25 * tmp_pij * tmp3_pij * dwik / rik;
          fix -= pf_dN * rikx;
          fiy -= pf_dN * riky;
          fiz -= pf_dN * rikz;
          fkx += pf_dN * rikx;
          fky += pf_dN * riky;
          fkz += pf_dN * rikz;

          // PijS forces
          const double pf_pij = VA * 0.25 * tmp_pij * dN2_pij[ktype] * dwik / rik;
          fix -= pf_pij * rikx;
          fiy -= pf_pij * riky;
          fiz -= pf_pij * rikz;
          fkx += pf_pij * rikx;
          fky += pf_pij * riky;
          fkz += pf_pij * rikz;

          _fx += fix;
          _fy += fiy;
          _fz += fiz;
          atomic_add_contribution(cells[cell_j][field::fx][p_j], fjx);
          atomic_add_contribution(cells[cell_j][field::fy][p_j], fjy);
          atomic_add_contribution(cells[cell_j][field::fz][p_j], fjz);
          const size_t cell_k = rebo_cell_k[kpos], p_k = rebo_p_k[kpos];
          atomic_add_contribution(cells[cell_k][field::fx][p_k], fkx);
          atomic_add_contribution(cells[cell_k][field::fy][p_k], fky);
          atomic_add_contribution(cells[cell_k][field::fz][p_k], fkz);
        }

        // ===================================================================
        // NjiC, NjiH — j's coordination (precomputed fields, excluding i)
        // ===================================================================
        const double NjiC = nC_j - (itype == 0 ? wij : 0.0);
        const double NjiH = nH_j - (itype == 1 ? wij : 0.0);
        const double Nji = NjiC + NjiH;

        // ===================================================================
        // pji — j's sigma-pi bond order to i: energy l-loop → compute pji → force l-loop
        // NconjtmpJ = Nconj[j] minus i's C contribution (O(1) with precomputed fields)
        // ===================================================================
        double unused_dSpNi;
        const double SpNi_J = rebo_Sp(nC_central + nH_central - wij, p.Nmin, p.Nmax, unused_dSpNi);
        const double NconjtmpJ = Nconj_j - (itype == 0 ? wij * SpNi_J : 0.0);

        // -----------------------------------------------------------------------
        // Build j's REBO-neighbour list once for this j: reused below by the
        // pji three-body sum (energy + force) and by the piRC / Tij "l-loop"
        // coordination corrections, instead of each of those four blocks
        // independently rescanning the whole (up to jnum-entry) buffer.
        // -----------------------------------------------------------------------
        int l_idx[REBO_MAX_REBO_NEIGHBORS];
        size_t l_cell[REBO_MAX_REBO_NEIGHBORS];
        size_t l_p[REBO_MAX_REBO_NEIGHBORS];
        double l_jlx[REBO_MAX_REBO_NEIGHBORS], l_jly[REBO_MAX_REBO_NEIGHBORS], l_jlz[REBO_MAX_REBO_NEIGHBORS];
        double l_rjl[REBO_MAX_REBO_NEIGHBORS];
        double l_wjl[REBO_MAX_REBO_NEIGHBORS], l_dwjl[REBO_MAX_REBO_NEIGHBORS];
        double l_SpNlj[REBO_MAX_REBO_NEIGHBORS], l_dNlj[REBO_MAX_REBO_NEIGHBORS];
        int l_count = 0;

        for (int ll = 0; ll < jnum; ++ll)
        {
          if (ll == jj)
            continue;
          const int ltype = buf.ext.type[ll];
          const double jlx = buf.drx[ll] - buf.drx[jj];
          const double jly = buf.dry[ll] - buf.dry[jj];
          const double jlz = buf.drz[ll] - buf.drz[jj];
          const double rjlsq = jlx * jlx + jly * jly + jlz * jlz;
          if (rjlsq >= p.rcmaxsq[jtype][ltype])
            continue;
          if (l_count >= int(REBO_MAX_REBO_NEIGHBORS))
          {
            assert(false);
            continue;
          }

          const double rjl = std::sqrt(rjlsq);
          double dwjl;
          const double wjl = rebo_Sp(rjl, p.rcmin[jtype][ltype], p.rcmax[jtype][ltype], dwjl);

          const int lpos = l_count++;
          l_idx[lpos] = ll;
          buf.nbh.get(ll, l_cell[lpos], l_p[lpos]);
          l_jlx[lpos] = jlx;
          l_jly[lpos] = jly;
          l_jlz[lpos] = jlz;
          l_rjl[lpos] = rjl;
          l_wjl[lpos] = wjl;
          l_dwjl[lpos] = dwjl;

          if (ltype == 0)
          {
            const double Nlj = cells[l_cell[lpos]][m_nijc_field][l_p[lpos]] + cells[l_cell[lpos]][m_nijh_field][l_p[lpos]] - wjl;
            l_SpNlj[lpos] = rebo_Sp(Nlj, p.Nmin, p.Nmax, l_dNlj[lpos]);
          }
          else
          {
            l_SpNlj[lpos] = 0.0;
            l_dNlj[lpos] = 0.0;
          }
        }

        // g/dgdc/exp_lam per l, cached here so the force l-loop below doesn't
        // repeat the branchy rebo_gSpline eval and exp() call for the same
        // (j,l) pair.
        double lj_g[REBO_MAX_REBO_NEIGHBORS], lj_dgdc[REBO_MAX_REBO_NEIGHBORS], lj_exp_lam[REBO_MAX_REBO_NEIGHBORS];
        double Etmp_pji = 0.0, tmp3_pji = 0.0;

        for (int lpos = 0; lpos < l_count; ++lpos)
        {
          const int ll = l_idx[lpos];
          const int ltype = buf.ext.type[ll];
          const double jlx = l_jlx[lpos], jly = l_jly[lpos], jlz = l_jlz[lpos];
          const double rjl = l_rjl[lpos];
          const double wjl = l_wjl[lpos];
          // lambda_ijl (nonzero when j is H)
          const double lam = (jtype == 1) ? 4.0 * ((p.rho[ltype][1] - rjl) - (p.rho[itype][1] - rij)) : 0.0;
          const double exp_lam = std::exp(lam);
          lj_exp_lam[lpos] = exp_lam;

          // cos_ijl: angle at j between bonds j-i and j-l
          double cosijl = (delx * jlx + dely * jly + delz * jlz) / (rij * rjl);
          cosijl = min(1.0, max(-1.0, cosijl));

          double dgdN;
          const double g = rebo_gSpline(cosijl, Nji, jtype, p, lj_dgdc[lpos], dgdN);
          lj_g[lpos] = g;
          Etmp_pji += wjl * g * exp_lam;
          tmp3_pji += wjl * dgdN * exp_lam;
        }

        double dN2_pji[2];
        const double PjiS = rebo_PijSpline(NjiC, NjiH, jtype, itype, p, dN2_pji);
        const double pji = 1.0 / std::sqrt(max(1.0 + Etmp_pji + PjiS, TOL));
        const double tmp_pji = -0.5 * pji * pji * pji;

        // pji three-body forces (l-loop, ×0.5 for full-list)
        for (int lpos = 0; lpos < l_count; ++lpos)
        {
          const int ll = l_idx[lpos];
          const int ltype = buf.ext.type[ll];
          const double jlx = l_jlx[lpos], jly = l_jly[lpos], jlz = l_jlz[lpos];
          const double rjl = l_rjl[lpos];
          const double wjl = l_wjl[lpos];
          const double dwjl = l_dwjl[lpos];
          const double exp_lam = lj_exp_lam[lpos];

          double cosijl = (delx * jlx + dely * jly + delz * jlz) / (rij * rjl);
          cosijl = min(1.0, max(-1.0, cosijl));

          const double g = lj_g[lpos];
          const double dgdc = lj_dgdc[lpos];

          const double inv_rij = 1.0 / rij, inv_rjl = 1.0 / rjl;
          const double inv_rij2 = inv_rij * inv_rij, inv_rjl2 = inv_rjl * inv_rjl;
          const double irr = inv_rij * inv_rjl;

          // dcosijl/dr* (see derivation in header comment)
          const double dcix = jlx * irr - cosijl * delx * inv_rij2;
          const double dciy = jly * irr - cosijl * dely * inv_rij2;
          const double dciz = jlz * irr - cosijl * delz * inv_rij2;
          const double dcjx = (-delx - jlx) * irr + cosijl * (delx * inv_rij2 + jlx * inv_rjl2);
          const double dcjy = (-dely - jly) * irr + cosijl * (dely * inv_rij2 + jly * inv_rjl2);
          const double dcjz = (-delz - jlz) * irr + cosijl * (delz * inv_rij2 + jlz * inv_rjl2);
          const double dclx = delx * irr - cosijl * jlx * inv_rjl2;
          const double dcly = dely * irr - cosijl * jly * inv_rjl2;
          const double dclz = delz * irr - cosijl * jlz * inv_rjl2;

          const double pf_ang = VA * 0.25 * tmp_pji * wjl * dgdc * exp_lam;
          double fix = -pf_ang * dcix, fiy = -pf_ang * dciy, fiz = -pf_ang * dciz;
          double fjx = -pf_ang * dcjx, fjy = -pf_ang * dcjy, fjz = -pf_ang * dcjz;
          double flx = -pf_ang * dclx, fly = -pf_ang * dcly, flz = -pf_ang * dclz;

          if (jtype == 1)
          {
            const double pf_lam = VA * 0.25 * tmp_pji * wjl * g * exp_lam * 4.0;
            fix -= pf_lam * (delx * inv_rij);
            fiy -= pf_lam * (dely * inv_rij);
            fiz -= pf_lam * (delz * inv_rij);
            fjx -= pf_lam * (jlx * inv_rjl - delx * inv_rij);
            fjy -= pf_lam * (jly * inv_rjl - dely * inv_rij);
            fjz -= pf_lam * (jlz * inv_rjl - delz * inv_rij);
            flx -= pf_lam * (-jlx * inv_rjl);
            fly -= pf_lam * (-jly * inv_rjl);
            flz -= pf_lam * (-jlz * inv_rjl);
          }

          // rjl = r_j - r_l = -jlx; d(rjl_mag)/d(r_j) = -jlx/rjl → fj += pf*jlx, fl -= pf*jlx
          const double pf_wjl = VA * 0.25 * tmp_pji * dwjl * g * exp_lam / rjl;
          fjx += pf_wjl * jlx;
          fjy += pf_wjl * jly;
          fjz += pf_wjl * jlz;
          flx -= pf_wjl * jlx;
          fly -= pf_wjl * jly;
          flz -= pf_wjl * jlz;

          const double pf_dN = VA * 0.25 * tmp_pji * tmp3_pji * dwjl / rjl;
          fjx += pf_dN * jlx;
          fjy += pf_dN * jly;
          fjz += pf_dN * jlz;
          flx -= pf_dN * jlx;
          fly -= pf_dN * jly;
          flz -= pf_dN * jlz;

          // PjiS forces
          const double pf_pji = VA * 0.25 * tmp_pji * dN2_pji[ltype] * dwjl / rjl;
          fjx += pf_pji * jlx;
          fjy += pf_pji * jly;
          fjz += pf_pji * jlz;
          flx -= pf_pji * jlx;
          fly -= pf_pji * jly;
          flz -= pf_pji * jlz;

          _fx += fix;
          _fy += fiy;
          _fz += fiz;
          atomic_add_contribution(cells[cell_j][field::fx][p_j], fjx);
          atomic_add_contribution(cells[cell_j][field::fy][p_j], fjy);
          atomic_add_contribution(cells[cell_j][field::fz][p_j], fjz);
          const size_t cell_l = l_cell[lpos], p_l = l_p[lpos];
          atomic_add_contribution(cells[cell_l][field::fx][p_l], flx);
          atomic_add_contribution(cells[cell_l][field::fy][p_l], fly);
          atomic_add_contribution(cells[cell_l][field::fz][p_l], flz);
        }

        // ===================================================================
        // piRC: spline → force k-loop + l-loop
        // ===================================================================
        const double Nijconj = 1.0 + NconjtmpI * NconjtmpI + NconjtmpJ * NconjtmpJ;

        double dN3_pi[3];
        const double piRC = rebo_piRCSpline(Nij, Nji, Nijconj, itype, jtype, p, dN3_pi);

        // ===================================================================
        // piRC forces (×0.5 for full-list)
        // dN3[0]: d(piRC)/d(Nij) → k-loop
        // dN3[1]: d(piRC)/d(Nji) → l-loop
        // dN3[2]: d(piRC)/d(Nijconj) → k-loop (NconjtmpI) + l-loop (NconjtmpJ)
        // ===================================================================
        {
          // k-loop: dN3[0] and dN3[2] forces
          for (int kpos = 0; kpos < rebo_count; ++kpos)
          {
            const int kk = rebo_idx[kpos];
            if (kk == jj)
              continue;
            const int ktype = buf.ext.type[kk];
            const double rik = rebo_r[kpos];
            double dwik;
            const double wik = rebo_w[kpos];
            dwik = rebo_dw[kpos];
            const double rikx = -buf.drx[kk], riky = -buf.dry[kk], rikz = -buf.drz[kk];
            const double pf0 = 0.5 * VA * dN3_pi[0] * dwik / rik;
            _fx -= pf0 * rikx;
            _fy -= pf0 * riky;
            _fz -= pf0 * rikz;
            const size_t cell_k = rebo_cell_k[kpos], p_k = rebo_p_k[kpos];
            atomic_add_contribution(cells[cell_k][field::fx][p_k], pf0 * rikx);
            atomic_add_contribution(cells[cell_k][field::fy][p_k], pf0 * riky);
            atomic_add_contribution(cells[cell_k][field::fz][p_k], pf0 * rikz);
            if (ktype == 0)
            {
              const double SpNki = rebo_SpNk[kpos];
              const double dNki = rebo_dNk[kpos];
              const double pf2a = 0.5 * VA * dN3_pi[2] * 2.0 * NconjtmpI * dwik * SpNki / rik;
              _fx -= pf2a * rikx;
              _fy -= pf2a * riky;
              _fz -= pf2a * rikz;
              atomic_add_contribution(cells[cell_k][field::fx][p_k], pf2a * rikx);
              atomic_add_contribution(cells[cell_k][field::fy][p_k], pf2a * riky);
              atomic_add_contribution(cells[cell_k][field::fz][p_k], pf2a * rikz);
              if (std::fabs(dNki) > TOL)
              {
                for (int mm = 0; mm < jnum; ++mm)
                {
                  if (mm == kk)
                    continue;
                  const int mtype = buf.ext.type[mm];
                  const double kmx = buf.drx[mm] - buf.drx[kk], kmy = buf.dry[mm] - buf.dry[kk], kmz = buf.drz[mm] - buf.drz[kk];
                  const double rkmsq = kmx * kmx + kmy * kmy + kmz * kmz;
                  if (rkmsq >= p.rcmaxsq[ktype][mtype])
                    continue;
                  const double rkm = std::sqrt(rkmsq);
                  if (rkm <= TOL)
                    continue;
                  double dwkm;
                  rebo_Sp(rkm, p.rcmin[ktype][mtype], p.rcmax[ktype][mtype], dwkm);
                  const double pf2b = 0.5 * VA * dN3_pi[2] * 2.0 * NconjtmpI * wik * dNki * dwkm / rkm;
                  size_t cell_m, p_m;
                  buf.nbh.get(mm, cell_m, p_m);
                  atomic_add_contribution(cells[cell_k][field::fx][p_k], pf2b * kmx);
                  atomic_add_contribution(cells[cell_k][field::fy][p_k], pf2b * kmy);
                  atomic_add_contribution(cells[cell_k][field::fz][p_k], pf2b * kmz);
                  atomic_add_contribution(cells[cell_m][field::fx][p_m], -pf2b * kmx);
                  atomic_add_contribution(cells[cell_m][field::fy][p_m], -pf2b * kmy);
                  atomic_add_contribution(cells[cell_m][field::fz][p_m], -pf2b * kmz);
                }
              }
            }
          }
          // l-loop: dN3[1] and dN3[2] forces
          for (int lpos = 0; lpos < l_count; ++lpos)
          {
            const int ll = l_idx[lpos];
            const int ltype = buf.ext.type[ll];
            const double jlx = l_jlx[lpos], jly = l_jly[lpos], jlz = l_jlz[lpos];
            const double rjl = l_rjl[lpos];
            const double wjl = l_wjl[lpos];
            const double dwjl = l_dwjl[lpos];
            // jlx = r_l-r_j = -rjl_LAMMPS → fj += +pf*jlx, fl += -pf*jlx
            const double pf1 = 0.5 * VA * dN3_pi[1] * dwjl / rjl;
            atomic_add_contribution(cells[cell_j][field::fx][p_j], pf1 * jlx);
            atomic_add_contribution(cells[cell_j][field::fy][p_j], pf1 * jly);
            atomic_add_contribution(cells[cell_j][field::fz][p_j], pf1 * jlz);
            const size_t cell_l = l_cell[lpos], p_l = l_p[lpos];
            atomic_add_contribution(cells[cell_l][field::fx][p_l], -pf1 * jlx);
            atomic_add_contribution(cells[cell_l][field::fy][p_l], -pf1 * jly);
            atomic_add_contribution(cells[cell_l][field::fz][p_l], -pf1 * jlz);
            if (ltype == 0)
            {
              const double SpNlj = l_SpNlj[lpos];
              const double dNlj = l_dNlj[lpos];
              const double pf2a = 0.5 * VA * dN3_pi[2] * 2.0 * NconjtmpJ * dwjl * SpNlj / rjl;
              atomic_add_contribution(cells[cell_j][field::fx][p_j], pf2a * jlx);
              atomic_add_contribution(cells[cell_j][field::fy][p_j], pf2a * jly);
              atomic_add_contribution(cells[cell_j][field::fz][p_j], pf2a * jlz);
              atomic_add_contribution(cells[cell_l][field::fx][p_l], -pf2a * jlx);
              atomic_add_contribution(cells[cell_l][field::fy][p_l], -pf2a * jly);
              atomic_add_contribution(cells[cell_l][field::fz][p_l], -pf2a * jlz);
              if (std::fabs(dNlj) > TOL)
              {
                // lmx = r_m-r_l = -rlm_LAMMPS → fl += +pf*lmx, fm += -pf*lmx
                for (int mm = 0; mm < jnum; ++mm)
                {
                  if (mm == ll || mm == jj)
                    continue;
                  const int mtype = buf.ext.type[mm];
                  const double lmx = buf.drx[mm] - buf.drx[ll], lmy = buf.dry[mm] - buf.dry[ll], lmz = buf.drz[mm] - buf.drz[ll];
                  const double rlmsq = lmx * lmx + lmy * lmy + lmz * lmz;
                  if (rlmsq >= p.rcmaxsq[ltype][mtype])
                    continue;
                  const double rlm = std::sqrt(rlmsq);
                  if (rlm <= TOL)
                    continue;
                  double dwlm;
                  rebo_Sp(rlm, p.rcmin[ltype][mtype], p.rcmax[ltype][mtype], dwlm);
                  const double pf2b = 0.5 * VA * dN3_pi[2] * 2.0 * NconjtmpJ * wjl * dNlj * dwlm / rlm;
                  size_t cell_m, p_m;
                  buf.nbh.get(mm, cell_m, p_m);
                  atomic_add_contribution(cells[cell_l][field::fx][p_l], pf2b * lmx);
                  atomic_add_contribution(cells[cell_l][field::fy][p_l], pf2b * lmy);
                  atomic_add_contribution(cells[cell_l][field::fz][p_l], pf2b * lmz);
                  atomic_add_contribution(cells[cell_m][field::fx][p_m], -pf2b * lmx);
                  atomic_add_contribution(cells[cell_m][field::fy][p_m], -pf2b * lmy);
                  atomic_add_contribution(cells[cell_m][field::fz][p_m], -pf2b * lmz);
                }
                // m = i: LAMMPS l-loop excludes only j, so i is included; i is absent from
                // own buffer so must be handled separately (rln = r_l - r_i = buf.dr[ll])
                {
                  const double rlim = std::sqrt(buf.d2[ll]);
                  if (rlim > TOL)
                  {
                    double dwlim;
                    rebo_Sp(rlim, p.rcmin[ltype][itype], p.rcmax[ltype][itype], dwlim);
                    const double limx = -buf.drx[ll], limy = -buf.dry[ll], limz = -buf.drz[ll];
                    const double pf2b = 0.5 * VA * dN3_pi[2] * 2.0 * NconjtmpJ * wjl * dNlj * dwlim / rlim;
                    atomic_add_contribution(cells[cell_l][field::fx][p_l], pf2b * limx);
                    atomic_add_contribution(cells[cell_l][field::fy][p_l], pf2b * limy);
                    atomic_add_contribution(cells[cell_l][field::fz][p_l], pf2b * limz);
                    _fx -= pf2b * limx;
                    _fy -= pf2b * limy;
                    _fz -= pf2b * limz;
                  }
                }
              }
            }
          }
        }

        // Tij (CC pairs only)
        double dN3_tij[3] = {0.0, 0.0, 0.0};
        double Tij = 0.0;
        double Etmp_Tij = 0.0;
        if (itype == 0 && jtype == 0)
          Tij = rebo_TijSpline(Nij, Nji, Nijconj, p, dN3_tij);

        // 4-body dihedral sum for Tij energy and forces
        // atom2=i, atom3=j, atom1=k (i's neighbor), atom4=l (j's neighbor)
        // r32 = rj-ri = buf.dr[jj], r23 = ri-rj = del, r23mag = rij
        if (std::fabs(Tij) > TOL)
        {
          const double r32x = buf.drx[jj], r32y = buf.dry[jj], r32z = buf.drz[jj];
          const double r23x = delx, r23y = dely, r23z = delz;
          const double rij2 = rsq;

          for (int kpos = 0; kpos < rebo_count; ++kpos)
          {
            const int kk = rebo_idx[kpos];
            if (kk == jj)
              continue;
            const int ktype = buf.ext.type[kk];
            // r21 = ri - rk  (LAMMPS: atom2 - atom1)
            const double r21x = -buf.drx[kk], r21y = -buf.dry[kk], r21z = -buf.drz[kk];
            const double r21mag = rebo_r[kpos];
            if (r21mag <= TOL)
              continue;

            const double cos321 = -(r21x * r32x + r21y * r32y + r21z * r32z) / (r21mag * rij);
            const double cos321c = min(1.0, max(-1.0, cos321));
            const double sin321 = std::sqrt(max(0.0, 1.0 - cos321c * cos321c));
            if (sin321 <= TOL)
              continue;

            double dcut321;
            const double rr_ik = 0.5 * (rij2 + buf.d2[kk] - (r21x - r23x) * (r21x - r23x) - (r21y - r23y) * (r21y - r23y) - (r21z - r23z) * (r21z - r23z)) / (rij * r21mag);
            const double tspjik = rebo_Sp2(rr_ik, thmin, thmax, dcut321);

            double dw21;
            const double w21 = rebo_Sp(r21mag, p.rcmin[itype][ktype], p.rcmaxp[itype][ktype], dw21);

            for (int ll = 0; ll < jnum; ++ll)
            {
              if (ll == jj || ll == kk)
                continue;
              const int ltype = buf.ext.type[ll];
              // r34 = rj - rl  (LAMMPS: atom3 - atom4)
              const double r34x = buf.drx[jj] - buf.drx[ll];
              const double r34y = buf.dry[jj] - buf.dry[ll];
              const double r34z = buf.drz[jj] - buf.drz[ll];
              const double r34sq = r34x * r34x + r34y * r34y + r34z * r34z;
              if (r34sq >= p.rcmaxsq[jtype][ltype])
                continue;
              const double r34mag = std::sqrt(r34sq);
              if (r34mag <= TOL)
                continue;

              const double cos234 = (r32x * r34x + r32y * r34y + r32z * r34z) / (rij * r34mag);
              const double cos234c = min(1.0, max(-1.0, cos234));
              const double sin234 = std::sqrt(max(0.0, 1.0 - cos234c * cos234c));
              if (sin234 <= TOL)
                continue;

              double dcut234;
              const double rr_jl = 0.5 * (rij2 + r34mag * r34mag - (r23x + r34x) * (r23x + r34x) - (r23y + r34y) * (r23y + r34y) - (r23z + r34z) * (r23z + r34z)) / (rij * r34mag);
              const double tspijl = rebo_Sp2(rr_jl, thmin, thmax, dcut234);

              double dw34;
              const double w34 = rebo_Sp(r34mag, p.rcmin[jtype][ltype], p.rcmaxp[jtype][ltype], dw34);

              // dihedral: cross(r32, r21) · cross(r23, r34) / denom
              const double cx321_0 = r32y * r21z - r32z * r21y;
              const double cx321_1 = r32z * r21x - r32x * r21z;
              const double cx321_2 = r32x * r21y - r32y * r21x;
              const double cx234_0 = r23y * r34z - r23z * r34y;
              const double cx234_1 = r23z * r34x - r23x * r34z;
              const double cx234_2 = r23x * r34y - r23y * r34x;

              const double cwnum = cx321_0 * cx234_0 + cx321_1 * cx234_1 + cx321_2 * cx234_2;
              const double cwnom = r21mag * r34mag * rij * rij * sin321 * sin234;
              const double om1234 = cwnum / cwnom;
              const double sin2om = 1.0 - om1234 * om1234;

              Etmp_Tij += sin2om * w21 * w34 * (1.0 - tspjik) * (1.0 - tspijl);

              // --- Tij forces (dihedral geometry) ---
              const double prefactor = VA * Tij;
              const double rjk_x = r21x - r23x, rjk_y = r21y - r23y, rjk_z = r21z - r23z;
              const double ril_x = r23x + r34x, ril_y = r23y + r34y, ril_z = r23z + r34z;

              const double sink2i = 1.0 / (sin321 * sin321);
              const double sinl2i = 1.0 / (sin234 * sin234);
              const double rik2i = 1.0 / (r21mag * r21mag);
              const double rjl2i = 1.0 / (r34mag * r34mag);

              const double rijrik = 2.0 * rij * r21mag;
              const double rijrjl = 2.0 * rij * r34mag;
              const double rr21 = rij2 - buf.d2[kk];
              const double rr34 = rij2 - r34mag * r34mag;
              const double rjk2 = rjk_x * rjk_x + rjk_y * rjk_y + rjk_z * rjk_z;
              const double ril2 = ril_x * ril_x + ril_y * ril_y + ril_z * ril_z;

              const double dctik = (-rr21 + rjk2) / (rijrik * buf.d2[kk]);
              const double dctij = (rr21 + rjk2) / (rijrik * rij2);
              const double dctjk = -2.0 / rijrik;
              const double dctjl = (-rr34 + ril2) / (rijrjl * r34mag * r34mag);
              const double dctji = (rr34 + ril2) / (rijrjl * rij2);
              const double dctil = -2.0 / rijrjl;

              const double aa = (prefactor * 2.0 * om1234 / cwnom) * w21 * w34 * (1.0 - tspjik) * (1.0 - tspijl);
              const double aaa2 = -prefactor * sin2om * w21 * w34;
              const double at2 = aa * cwnum;

              const double dt1dik = rik2i - dctik * sink2i * cos321c;
              const double dt1djk = -dctjk * sink2i * cos321c;
              const double dt1djl = rjl2i - dctjl * sinl2i * cos234c;
              const double dt1dil = -dctil * sinl2i * cos234c;
              const double dt1dij = 2.0 / rij2 - dctij * sink2i * cos321c - dctji * sinl2i * cos234c;

              const double dt2dik_0 = -r23z * cx234_1 + r23y * cx234_2;
              const double dt2dik_1 = -r23x * cx234_2 + r23z * cx234_0;
              const double dt2dik_2 = -r23y * cx234_0 + r23x * cx234_1;
              const double dt2djl_0 = -r23y * cx321_2 + r23z * cx321_1;
              const double dt2djl_1 = -r23z * cx321_0 + r23x * cx321_2;
              const double dt2djl_2 = -r23x * cx321_1 + r23y * cx321_0;
              const double dt2dij_0 = r21z * cx234_1 - r34z * cx321_1 - r21y * cx234_2 + r34y * cx321_2;
              const double dt2dij_1 = r21x * cx234_2 - r34x * cx321_2 - r21z * cx234_0 + r34z * cx321_0;
              const double dt2dij_2 = r21y * cx234_0 - r34y * cx321_0 - r21x * cx234_1 + r34x * cx321_1;

              const double fcijpc = (-dt1dij * at2) + aaa2 * (-dcut321 * dctij * (1.0 - tspijl) - dcut234 * dctji * (1.0 - tspjik));
              const double fcikpc = (-dt1dik * at2) + aaa2 * (-dcut321 * dctik * (1.0 - tspijl));
              const double fcjlpc = (-dt1djl * at2) + aaa2 * (-dcut234 * dctjl * (1.0 - tspjik));
              const double fcjkpc = (-dt1djk * at2) + aaa2 * (-dcut321 * dctjk * (1.0 - tspijl));
              const double fcilpc = (-dt1dil * at2) + aaa2 * (-dcut234 * dctil * (1.0 - tspjik));

              // F12 (on k=atom1), F23 (on ij bond), F34 (on l=atom4)
              double F12_0 = fcikpc * r21x + aa * dt2dik_0;
              double F12_1 = fcikpc * r21y + aa * dt2dik_1;
              double F12_2 = fcikpc * r21z + aa * dt2dik_2;
              double F23_0 = fcijpc * r23x + aa * dt2dij_0;
              double F23_1 = fcijpc * r23y + aa * dt2dij_1;
              double F23_2 = fcijpc * r23z + aa * dt2dij_2;
              double F34_0 = fcjlpc * r34x + aa * dt2djl_0;
              double F34_1 = fcjlpc * r34y + aa * dt2djl_1;
              double F34_2 = fcjlpc * r34z + aa * dt2djl_2;
              const double F31_0 = fcjkpc * rjk_x;
              const double F31_1 = fcjkpc * rjk_y;
              const double F31_2 = fcjkpc * rjk_z;
              const double F24_0 = fcilpc * ril_x;
              const double F24_1 = fcilpc * ril_y;
              const double F24_2 = fcilpc * ril_z;

              // f1(k), f2(i), f3(j), f4(l)
              double fk0 = -F12_0 - F31_0, fk1 = -F12_1 - F31_1, fk2 = -F12_2 - F31_2;
              double fi0 = F23_0 + F12_0 + F24_0;
              double fi1 = F23_1 + F12_1 + F24_1;
              double fi2 = F23_2 + F12_2 + F24_2;
              double fj0 = -F23_0 + F34_0 + F31_0;
              double fj1 = -F23_1 + F34_1 + F31_1;
              double fj2 = -F23_2 + F34_2 + F31_2;
              double fl0 = -F34_0 - F24_0, fl1 = -F34_1 - F24_1, fl2 = -F34_2 - F24_2;

              // w-coordination forces
              const double pf_w21 = VA * Tij * sin2om * (1.0 - tspjik) * (1.0 - tspijl) * dw21 * w34 / r21mag;
              fi0 -= pf_w21 * r21x;
              fi1 -= pf_w21 * r21y;
              fi2 -= pf_w21 * r21z;
              fk0 += pf_w21 * r21x;
              fk1 += pf_w21 * r21y;
              fk2 += pf_w21 * r21z;
              const double pf_w34 = VA * Tij * sin2om * (1.0 - tspjik) * (1.0 - tspijl) * w21 * dw34 / r34mag;
              fj0 -= pf_w34 * r34x;
              fj1 -= pf_w34 * r34y;
              fj2 -= pf_w34 * r34z;
              size_t cell_l, p_l;
              buf.nbh.get(ll, cell_l, p_l);
              const double hl2 = 0.5;
              _fx += hl2 * fi0;
              _fy += hl2 * fi1;
              _fz += hl2 * fi2;
              atomic_add_contribution(cells[cell_j][field::fx][p_j], hl2 * fj0);
              atomic_add_contribution(cells[cell_j][field::fy][p_j], hl2 * fj1);
              atomic_add_contribution(cells[cell_j][field::fz][p_j], hl2 * fj2);
              const size_t cell_k = rebo_cell_k[kpos], p_k = rebo_p_k[kpos];
              atomic_add_contribution(cells[cell_k][field::fx][p_k], hl2 * fk0);
              atomic_add_contribution(cells[cell_k][field::fy][p_k], hl2 * fk1);
              atomic_add_contribution(cells[cell_k][field::fz][p_k], hl2 * fk2);
              fl0 += pf_w34 * r34x;
              fl1 += pf_w34 * r34y;
              fl2 += pf_w34 * r34z; // consolidate w34 on l
              atomic_add_contribution(cells[cell_l][field::fx][p_l], hl2 * fl0);
              atomic_add_contribution(cells[cell_l][field::fy][p_l], hl2 * fl1);
              atomic_add_contribution(cells[cell_l][field::fz][p_l], hl2 * fl2);
            }
          }
        }

        // ===================================================================
        // Tij dN3 forces (×0.5)
        // dN3[0]: d(Tij)/d(Nij)     → k-loop (dwik)
        // dN3[1]: d(Tij)/d(Nji)     → l-loop (dwjl)
        // dN3[2]: d(Tij)/d(Nijconj) → k-loop (NconjtmpI) + l-loop (NconjtmpJ)
        // ===================================================================
        if (std::fabs(Tij) > TOL && std::fabs(Etmp_Tij) > TOL)
        {
          for (int kpos = 0; kpos < rebo_count; ++kpos)
          {
            const int kk = rebo_idx[kpos];
            if (kk == jj)
              continue;
            const int ktype = buf.ext.type[kk];
            const double rik = rebo_r[kpos];
            double dwik;
            const double wik = rebo_w[kpos];
            dwik = rebo_dw[kpos];
            const double rikx = -buf.drx[kk], riky = -buf.dry[kk], rikz = -buf.drz[kk];
            const size_t cell_k = rebo_cell_k[kpos], p_k = rebo_p_k[kpos];

            const double pf0 = 0.5 * VA * dN3_tij[0] * dwik * Etmp_Tij / rik;
            _fx -= pf0 * rikx;
            _fy -= pf0 * riky;
            _fz -= pf0 * rikz;
            atomic_add_contribution(cells[cell_k][field::fx][p_k], pf0 * rikx);
            atomic_add_contribution(cells[cell_k][field::fy][p_k], pf0 * riky);
            atomic_add_contribution(cells[cell_k][field::fz][p_k], pf0 * rikz);

            if (ktype == 0)
            {
              const double SpNki = rebo_SpNk[kpos];
              const double dNki = rebo_dNk[kpos];
              const double pf2a = 0.5 * VA * dN3_tij[2] * Etmp_Tij * 2.0 * NconjtmpI * dwik * SpNki / rik;
              _fx -= pf2a * rikx;
              _fy -= pf2a * riky;
              _fz -= pf2a * rikz;
              atomic_add_contribution(cells[cell_k][field::fx][p_k], pf2a * rikx);
              atomic_add_contribution(cells[cell_k][field::fy][p_k], pf2a * riky);
              atomic_add_contribution(cells[cell_k][field::fz][p_k], pf2a * rikz);

              if (std::fabs(dNki) > TOL)
              {
                for (int mm = 0; mm < jnum; ++mm)
                {
                  if (mm == kk)
                    continue;
                  const int mtype = buf.ext.type[mm];
                  const double kmx = buf.drx[mm] - buf.drx[kk];
                  const double kmy = buf.dry[mm] - buf.dry[kk];
                  const double kmz = buf.drz[mm] - buf.drz[kk];
                  const double rkmsq = kmx * kmx + kmy * kmy + kmz * kmz;
                  if (rkmsq >= p.rcmaxsq[ktype][mtype])
                    continue;
                  const double rkm = std::sqrt(rkmsq);
                  if (rkm <= TOL)
                    continue;
                  double dwkm;
                  rebo_Sp(rkm, p.rcmin[ktype][mtype], p.rcmax[ktype][mtype], dwkm);
                  const double pf2b = 0.5 * VA * dN3_tij[2] * Etmp_Tij * 2.0 * NconjtmpI * wik * dNki * dwkm / rkm;
                  size_t cell_m, p_m;
                  buf.nbh.get(mm, cell_m, p_m);
                  atomic_add_contribution(cells[cell_k][field::fx][p_k], pf2b * kmx);
                  atomic_add_contribution(cells[cell_k][field::fy][p_k], pf2b * kmy);
                  atomic_add_contribution(cells[cell_k][field::fz][p_k], pf2b * kmz);
                  atomic_add_contribution(cells[cell_m][field::fx][p_m], -pf2b * kmx);
                  atomic_add_contribution(cells[cell_m][field::fy][p_m], -pf2b * kmy);
                  atomic_add_contribution(cells[cell_m][field::fz][p_m], -pf2b * kmz);
                }
              }
            }
          }
          for (int lpos = 0; lpos < l_count; ++lpos)
          {
            const int ll = l_idx[lpos];
            const int ltype = buf.ext.type[ll];
            const double jlx = l_jlx[lpos], jly = l_jly[lpos], jlz = l_jlz[lpos];
            const double rjl = l_rjl[lpos];
            const double wjl = l_wjl[lpos];
            const double dwjl = l_dwjl[lpos];
            const size_t cell_l = l_cell[lpos], p_l = l_p[lpos];

            // jlx = r_l-r_j = -rjl_LAMMPS → fj += +pf*jlx, fl += -pf*jlx
            const double pf1 = 0.5 * VA * dN3_tij[1] * dwjl * Etmp_Tij / rjl;
            atomic_add_contribution(cells[cell_j][field::fx][p_j], pf1 * jlx);
            atomic_add_contribution(cells[cell_j][field::fy][p_j], pf1 * jly);
            atomic_add_contribution(cells[cell_j][field::fz][p_j], pf1 * jlz);
            atomic_add_contribution(cells[cell_l][field::fx][p_l], -pf1 * jlx);
            atomic_add_contribution(cells[cell_l][field::fy][p_l], -pf1 * jly);
            atomic_add_contribution(cells[cell_l][field::fz][p_l], -pf1 * jlz);

            if (ltype == 0)
            {
              const double SpNlj = l_SpNlj[lpos];
              const double dNlj = l_dNlj[lpos];
              const double pf2a = 0.5 * VA * dN3_tij[2] * Etmp_Tij * 2.0 * NconjtmpJ * dwjl * SpNlj / rjl;
              atomic_add_contribution(cells[cell_j][field::fx][p_j], pf2a * jlx);
              atomic_add_contribution(cells[cell_j][field::fy][p_j], pf2a * jly);
              atomic_add_contribution(cells[cell_j][field::fz][p_j], pf2a * jlz);
              atomic_add_contribution(cells[cell_l][field::fx][p_l], -pf2a * jlx);
              atomic_add_contribution(cells[cell_l][field::fy][p_l], -pf2a * jly);
              atomic_add_contribution(cells[cell_l][field::fz][p_l], -pf2a * jlz);

              if (std::fabs(dNlj) > TOL)
              {
                // lmx = r_m-r_l = -rlm_LAMMPS → fl += +pf*lmx, fm += -pf*lmx
                for (int mm = 0; mm < jnum; ++mm)
                {
                  if (mm == ll || mm == jj)
                    continue;
                  const int mtype = buf.ext.type[mm];
                  const double lmx = buf.drx[mm] - buf.drx[ll];
                  const double lmy = buf.dry[mm] - buf.dry[ll];
                  const double lmz = buf.drz[mm] - buf.drz[ll];
                  const double rlmsq = lmx * lmx + lmy * lmy + lmz * lmz;
                  if (rlmsq >= p.rcmaxsq[ltype][mtype])
                    continue;
                  const double rlm = std::sqrt(rlmsq);
                  if (rlm <= TOL)
                    continue;
                  double dwlm;
                  rebo_Sp(rlm, p.rcmin[ltype][mtype], p.rcmax[ltype][mtype], dwlm);
                  const double pf2b = 0.5 * VA * dN3_tij[2] * Etmp_Tij * 2.0 * NconjtmpJ * wjl * dNlj * dwlm / rlm;
                  size_t cell_m, p_m;
                  buf.nbh.get(mm, cell_m, p_m);
                  atomic_add_contribution(cells[cell_l][field::fx][p_l], pf2b * lmx);
                  atomic_add_contribution(cells[cell_l][field::fy][p_l], pf2b * lmy);
                  atomic_add_contribution(cells[cell_l][field::fz][p_l], pf2b * lmz);
                  atomic_add_contribution(cells[cell_m][field::fx][p_m], -pf2b * lmx);
                  atomic_add_contribution(cells[cell_m][field::fy][p_m], -pf2b * lmy);
                  atomic_add_contribution(cells[cell_m][field::fz][p_m], -pf2b * lmz);
                }
                // m = i: LAMMPS l-loop excludes only j, so i is included; i is absent from
                // own buffer so must be handled separately (rln = r_l - r_i = buf.dr[ll])
                {
                  const double rlim = std::sqrt(buf.d2[ll]);
                  if (rlim > TOL)
                  {
                    double dwlim;
                    rebo_Sp(rlim, p.rcmin[ltype][itype], p.rcmax[ltype][itype], dwlim);
                    const double limx = -buf.drx[ll], limy = -buf.dry[ll], limz = -buf.drz[ll];
                    const double pf2b = 0.5 * VA * dN3_tij[2] * Etmp_Tij * 2.0 * NconjtmpJ * wjl * dNlj * dwlim / rlim;
                    atomic_add_contribution(cells[cell_l][field::fx][p_l], pf2b * limx);
                    atomic_add_contribution(cells[cell_l][field::fy][p_l], pf2b * limy);
                    atomic_add_contribution(cells[cell_l][field::fz][p_l], pf2b * limz);
                    _fx -= pf2b * limx;
                    _fy -= pf2b * limy;
                    _fz -= pf2b * limz;
                  }
                }
              }
            }
          }
        }

        // ===================================================================
        // Full bond order and energy
        // ===================================================================
        const double bij = 0.5 * (pij + pji) + piRC + Tij * Etmp_Tij;
        _en += 0.5 * (bij * VA);

        // --- Direct pair force (×0.5 for full-list); VR's contribution is
        // handled by REBOPairwiseForceOp. ---
        const double fpair = -(bij * dVA) / rij;
        _fx += 0.5 * delx * fpair;
        _fy += 0.5 * dely * fpair;
        _fz += 0.5 * delz * fpair;
        atomic_add_contribution(cells[cell_j][field::fx][p_j], -0.5 * delx * fpair);
        atomic_add_contribution(cells[cell_j][field::fy][p_j], -0.5 * dely * fpair);
        atomic_add_contribution(cells[cell_j][field::fz][p_j], -0.5 * delz * fpair);

      } // end j-loop

      atomic_add_contribution(en, _en);
      atomic_add_contribution(fx, _fx);
      atomic_add_contribution(fy, _fy);
      atomic_add_contribution(fz, _fz);
    }
  };

} // namespace exaStamp

// ComputePairTraits specialization must live in the exanb namespace
// where the primary template is defined.
namespace exanb
{
  template <> struct ComputePairTraits<exaStamp::REBOPairwiseForceOp>
  {
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = true;
    static inline constexpr bool HasParticleContextStart = false;
    static inline constexpr bool HasParticleContext = false;
    static inline constexpr bool HasParticleContextStop = false;
    static inline constexpr bool RequiresNbhOptionalData = true;
  };

  template <class NijcFieldT, class NijhFieldT, class NconjFieldT> struct ComputePairTraits<exaStamp::REBOManyBodyForceOp<NijcFieldT, NijhFieldT, NconjFieldT>>
  {
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = true;
    static inline constexpr bool HasParticleContextStart = false;
    static inline constexpr bool HasParticleContext = false;
    static inline constexpr bool HasParticleContextStop = false;
    static inline constexpr bool RequiresNbhOptionalData = true;
  };
} // namespace exanb
