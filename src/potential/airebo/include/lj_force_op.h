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
#include <onika/physics/units.h>
#include <exaStamp/unit_system.h>
#include <exanb/core/concurent_add_contributions.h>

#include "airebo_params.h"
#include "rebo_force_op_cached.h"   // rebo_Sp, rebo_Sp2, rebo_gSpline, spline helpers

namespace exaStamp
{
  using namespace exanb;

  // ===========================================================================
  // LJForceOp — Lennard-Jones non-bonded term of AIREBO (Stuart et al. 2000)
  //
  // For each pair (i,j) within rcLJmax[i][j]:
  //   1. Compute cij = 1 - best, where best = max REBO-path connectivity
  //      (2-path: direct wij; 3-path: i-k-j; 4-path: i-k-m-j).
  //      Skip the pair if cij == 0 (fully bonded).
  //   2. VLJ = 4*eps*((sig/r)^12-(sig/r)^6) * slw  (slw = outer cutoff switching)
  //   3. Str = Sp2(r, rcLJmin, rcLJmax)              (near-bond switching)
  //   4. When Str > 0: compute bij = 0.5*(pij+pji)+piRC+Tij via bondorderLJ
  //      (same as REBO bond-order but with wij evaluated at rcmin, angles at rij).
  //      Stb = Sp2(bij, bLJmin, bLJmax).
  //   5. E = cij*VLJ*(Str*Stb + (1-Str)).
  //      Force = -(d(Str)*(Stb-1)*cij*VLJ + dVLJ*(Str*Stb+1-Str)*cij)/rij.
  //   6. Path forces from dcij (npath = 2, 3, or 4).
  //   7. Bond-order forces from bij (VA_lj = Str*cij*VLJ, LAMMPS convention).
  //
  // Full-list: all contributions × 0.5.
  // ===========================================================================

  template<class NijcFieldT, class NijhFieldT, class NconjFieldT>  
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) LJForceOp
  {
    const AireboParamsRO* m_params = nullptr;
    NijcFieldT  m_nijc_field  = {};   // read neighbour's nC
    NijhFieldT  m_nijh_field  = {};   // read neighbour's nH
    NconjFieldT m_nconj_field = {};   // read neighbour's Nconj

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int type, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, nC_central, nH_central, Nconj_central, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int type, Mat3d& virial, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, nC_central, nH_central, Nconj_central, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t n, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int type, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells,
                           GridCellLocksT locks, ParticleLockT& lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, nC_central, nH_central, Nconj_central, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class Mat3dT,
             class GridCellLocksT, class ParticleLockT>
    inline void operator()(int jnum, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int type, Mat3dT& virial, double nC_central, double nH_central, double Nconj_central, CellParticlesT cells,
                           GridCellLocksT locks, ParticleLockT& lock_a) const
    {
      static constexpr double TOL   = 1.0e-9;
      static constexpr double thmin = -1.0, thmax = -0.995;
      static constexpr double conv  = 1.0; //EXASTAMP_CONST_QUANTITY( 1. * eV );

      double _en = 0.0, _fx = 0.0, _fy = 0.0, _fz = 0.0;
      const int itype = type;
      const AireboParamsRO& p = *m_params;

      const double cut3rebo = 3.0 * p.rcmax[0][0];

      // -----------------------------------------------------------------------
      // Main j-loop
      // -----------------------------------------------------------------------
      for (int jj = 0; jj < jnum; ++jj) {
        const int    jtype = buf.ext.type[jj];
        const double rsq   = buf.d2[jj];
        if (rsq >= p.rcLJmaxsq[itype][jtype]) continue;

        const double delx = -buf.drx[jj];
        const double dely = -buf.dry[jj];
        const double delz = -buf.drz[jj];
        const double rij  = std::sqrt(rsq);

        // ===================================================================
        // Step 1: Connectivity screening — find best REBO path i→j
        // ===================================================================
        double best   = 0.0;
        int best_npath = 0;

        // npath=2 storage
        double best_dwij2 = 0.0;

        // npath=3 storage (best 3-body path i-k-j)
        int    best_kk3   = -1;
        double best_delikx3 = 0.0, best_deliky3 = 0.0, best_delikz3 = 0.0;
        double best_rik3 = 0.0, best_wik3 = 0.0, best_dwik3 = 0.0;
        double best_deljkx3 = 0.0, best_deljky3 = 0.0, best_deljkz3 = 0.0;
        double best_rkj3 = 0.0, best_wkj3 = 0.0, best_dwkj3 = 0.0;

        // npath=4 storage (best 4-body path i-k-m-j)
        int    best_kk4 = -1, best_mm4 = -1;
        double best_delikx4 = 0.0, best_deliky4 = 0.0, best_delikz4 = 0.0;
        double best_rik4 = 0.0, best_wik4 = 0.0, best_dwik4 = 0.0;
        double best_delkmx4 = 0.0, best_delkmy4 = 0.0, best_delkmz4 = 0.0;
        double best_rkm4 = 0.0, best_wkm4 = 0.0, best_dwkm4 = 0.0;
        double best_deljmx4 = 0.0, best_deljmy4 = 0.0, best_deljmz4 = 0.0;
        double best_rmj4 = 0.0, best_wmj4 = 0.0, best_dwmj4 = 0.0;

        bool done = false;

        // 2-path (direct bond)
        if (rij < p.rcmax[itype][jtype]) {
          double dw2;
          best = rebo_Sp(rij, p.rcmin[itype][jtype], p.rcmax[itype][jtype], dw2);
          best_npath = 2; best_dwij2 = dw2;
          if (best >= 1.0) done = true;
        }

        // 3- and 4-path
        if (!done && rij < cut3rebo) {
          for (int kk = 0; kk < jnum && !done; ++kk) {
            if (kk == jj) continue;
            const int ktype = buf.ext.type[kk];
            const double rik = std::sqrt(buf.d2[kk]);
            double dwik;
            const double wik = rebo_Sp(rik, p.rcmin[itype][ktype], p.rcmax[itype][ktype], dwik);
            if (wik <= 0.0) continue;

            if (wik > best) {
              // 3-body path i-k-j
              const double dkjx = buf.drx[jj]-buf.drx[kk];
              const double dkjy = buf.dry[jj]-buf.dry[kk];
              const double dkjz = buf.drz[jj]-buf.drz[kk];
              const double rkj  = std::sqrt(dkjx*dkjx + dkjy*dkjy + dkjz*dkjz);
              double dwkj;
              const double wkj = rebo_Sp(rkj, p.rcmin[ktype][jtype], p.rcmax[ktype][jtype], dwkj);
              if (wik*wkj > best) {
                best = wik*wkj; best_npath = 3; best_kk3 = kk;
                best_delikx3 = -buf.drx[kk]; best_deliky3 = -buf.dry[kk]; best_delikz3 = -buf.drz[kk];
                best_rik3 = rik; best_wik3 = wik; best_dwik3 = dwik;
                best_deljkx3 = dkjx; best_deljky3 = dkjy; best_deljkz3 = dkjz;
                best_rkj3 = rkj; best_wkj3 = wkj; best_dwkj3 = dwkj;
                if (best >= 1.0) { done = true; break; }
              }

              // 4-body path i-k-m-j
              for (int mm = 0; mm < jnum && !done; ++mm) {
                if (mm == kk || mm == jj) continue;
                const int mtype = buf.ext.type[mm];
                const double kmx = buf.drx[mm]-buf.drx[kk];
                const double kmy = buf.dry[mm]-buf.dry[kk];
                const double kmz = buf.drz[mm]-buf.drz[kk];
                const double rkm = std::sqrt(kmx*kmx + kmy*kmy + kmz*kmz);
                double dwkm;
                const double wkm = rebo_Sp(rkm, p.rcmin[ktype][mtype], p.rcmax[ktype][mtype], dwkm);
                if (wik*wkm <= best) continue;
                const double djmx = buf.drx[jj]-buf.drx[mm];
                const double djmy = buf.dry[jj]-buf.dry[mm];
                const double djmz = buf.drz[jj]-buf.drz[mm];
                const double rmj  = std::sqrt(djmx*djmx + djmy*djmy + djmz*djmz);
                double dwmj;
                const double wmj = rebo_Sp(rmj, p.rcmin[mtype][jtype], p.rcmax[mtype][jtype], dwmj);
                if (wik*wkm*wmj > best) {
                  best = wik*wkm*wmj; best_npath = 4; best_kk4 = kk; best_mm4 = mm;
                  best_delikx4 = -buf.drx[kk]; best_deliky4 = -buf.dry[kk]; best_delikz4 = -buf.drz[kk];
                  best_rik4 = rik; best_wik4 = wik; best_dwik4 = dwik;
                  best_delkmx4 = kmx; best_delkmy4 = kmy; best_delkmz4 = kmz;
                  best_rkm4 = rkm; best_wkm4 = wkm; best_dwkm4 = dwkm;
                  best_deljmx4 = djmx; best_deljmy4 = djmy; best_deljmz4 = djmz;
                  best_rmj4 = rmj; best_wmj4 = wmj; best_dwmj4 = dwmj;
                  if (best >= 1.0) { done = true; break; }
                }
              }
            }
          }
        }

        const double cij = 1.0 - best;
        if (cij <= TOL) continue;

        // ===================================================================
        // Step 2: LJ potential with outer switching
        // ===================================================================
        const double sig  = p.sigma[itype][jtype];
        double dslw;
        const double slw  = rebo_Sp2(rij, p.sigmin*sig, p.cutlj*sig, dslw);

        const double r2inv = 1.0/rsq;
        const double r6inv = r2inv*r2inv*r2inv;
        const double sig2  = sig*sig;
        const double sig6  = sig2*sig2*sig2;
        const double sig12 = sig6*sig6;
        const double eps   = p.epsilon[itype][jtype];
        const double vdw   =  4.0*eps * r6inv * (sig12*r6inv - sig6);
        const double dvdw  = -4.0*eps * r6inv * (12.0*sig12*r6inv - 6.0*sig6) / rij;

        double VLJ  = vdw * slw;
        double dVLJ = dvdw*slw + vdw*dslw;
        // Convert eV → exaStamp internal units; all force/energy terms derive from VLJ
        VLJ *= conv; dVLJ *= conv;

        double dStr;
        const double Str = rebo_Sp2(rij, p.rcLJmin[itype][jtype], p.rcLJmax[itype][jtype], dStr);

        size_t cell_j, p_j;
        buf.nbh.get(jj, cell_j, p_j);

        // ===================================================================
        // Step 3: Bond-order bij (bondorderLJ — only when Str > 0)
        // wij at scaled distance rcmin → wij_bo = 1, dwij_bo = 0
        // angles use actual rij geometry; lambda uses rcmin[i][j]
        // ===================================================================
        double Stb = 0.0, dStb = 0.0;

        if (Str > TOL) {
          // wij evaluated at rcmin (scaled geometry), always = 1
          const double NijC = nC_central - (jtype == 0 ? 1.0 : 0.0);
          const double NijH = nH_central - (jtype == 1 ? 1.0 : 0.0);
          const double Nij  = NijC + NijH;
          const double rij_bo = p.rcmin[itype][jtype]; // for lambda

          // j's coordination from precomputed grid fields
          const double NjC  = cells[cell_j][m_nijc_field][p_j];
          const double NjH  = cells[cell_j][m_nijh_field][p_j];
          const double NjiC = NjC - (itype == 0 ? 1.0 : 0.0);
          const double NjiH = NjH - (itype == 1 ? 1.0 : 0.0);
          const double Nji  = NjiC + NjiH;

          // NconjtmpI/J from Nconj grid fields — O(1) via inclusion-exclusion on j's contribution
          double dw_ij_rebo;
          const double wij_rebo = rebo_Sp(rij, p.rcmin[itype][jtype], p.rcmax[itype][jtype], dw_ij_rebo);
          double dSpNj_conj, dSpNi_conj;
          const double NconjtmpI = Nconj_central
              - (jtype == 0 ? wij_rebo * rebo_Sp(NjC+NjH-wij_rebo, p.Nmin, p.Nmax, dSpNj_conj) : 0.0);
          const double Nconj_j   = cells[cell_j][m_nconj_field][p_j];
          const double NconjtmpJ = Nconj_j
              - (itype == 0 ? wij_rebo * rebo_Sp(nC_central+nH_central-wij_rebo, p.Nmin, p.Nmax, dSpNi_conj) : 0.0);

          // ---- pij: i's contribution to bond order ----
          double Etmp_pij = 0.0, tmp3_pij = 0.0;
          for (int kk = 0; kk < jnum; ++kk) {
            if (kk == jj) continue;
            const int    ktype = buf.ext.type[kk];
            const double rik   = std::sqrt(buf.d2[kk]);
            double dS;
            const double wik = rebo_Sp(rik, p.rcmin[itype][ktype], p.rcmax[itype][ktype], dS);
            if (wik <= TOL) continue;

            // lambda uses rcmin[i][j] (scaled) instead of actual rij
            const double lam = (itype==1)
              ? 4.0*((p.rho[ktype][1]-rik) - (p.rho[jtype][1]-rij_bo)) : 0.0;
            const double exp_lam = std::exp(lam);

            double cosjik = (buf.drx[jj]*buf.drx[kk] + buf.dry[jj]*buf.dry[kk]
                             + buf.drz[jj]*buf.drz[kk]) / (rij*rik);
            cosjik = std::min(1.0, std::max(-1.0, cosjik));

            double dgdc, dgdN;
            const double g = rebo_gSpline(cosjik, Nij, itype, p, dgdc, dgdN);
            Etmp_pij += wik * g * exp_lam;
            tmp3_pij += wik * dgdN * exp_lam;
          }

          double dN2_pij[2];
          const double PijS = rebo_PijSpline(NijC, NijH, itype, jtype, p, dN2_pij);
          const double pij   = 1.0 / std::sqrt(std::max(1.0 + Etmp_pij + PijS, TOL));
          const double tmp_pij = -0.5 * pij*pij*pij;

          // ---- pji: j's contribution to bond order ----
          double Etmp_pji = 0.0, tmp3_pji = 0.0;
          for (int ll = 0; ll < jnum; ++ll) {
            if (ll == jj) continue;
            const int    ltype = buf.ext.type[ll];
            const double jlx = buf.drx[ll]-buf.drx[jj];
            const double jly = buf.dry[ll]-buf.dry[jj];
            const double jlz = buf.drz[ll]-buf.drz[jj];
            const double rjl = std::sqrt(jlx*jlx+jly*jly+jlz*jlz);
            double dwjl;
            const double wjl = rebo_Sp(rjl, p.rcmin[jtype][ltype], p.rcmax[jtype][ltype], dwjl);
            if (wjl <= TOL) continue;

            // lambda uses rcmin[i][j] (scaled) instead of actual rij
            const double lam = (jtype==1)
              ? 4.0*((p.rho[ltype][1]-rjl) - (p.rho[itype][1]-rij_bo)) : 0.0;
            const double exp_lam = std::exp(lam);

            double cosijl = (delx*jlx + dely*jly + delz*jlz) / (rij*rjl);
            cosijl = std::min(1.0, std::max(-1.0, cosijl));

            double dgdc, dgdN;
            const double g = rebo_gSpline(cosijl, Nji, jtype, p, dgdc, dgdN);
            Etmp_pji += wjl * g * exp_lam;
            tmp3_pji += wjl * dgdN * exp_lam;
          }

          double dN2_pji[2];
          const double PjiS = rebo_PijSpline(NjiC, NjiH, jtype, itype, p, dN2_pji);
          const double pji   = 1.0 / std::sqrt(std::max(1.0 + Etmp_pji + PjiS, TOL));
          const double tmp_pji = -0.5 * pji*pji*pji;

          // ---- piRC and Tij ----
          const double Nijconj = 1.0 + NconjtmpI*NconjtmpI + NconjtmpJ*NconjtmpJ;
          double dN3_pi[3];
          const double piRC = rebo_piRCSpline(Nij, Nji, Nijconj, itype, jtype, p, dN3_pi);

          double dN3_tij[3] = {0.0, 0.0, 0.0};
          double Tij = 0.0, Etmp_Tij = 0.0;
          if (itype == 0 && jtype == 0)
            Tij = rebo_TijSpline(Nij, Nji, Nijconj, p, dN3_tij);

          if (std::fabs(Tij) > TOL) {
            const double r32x = buf.drx[jj], r32y = buf.dry[jj], r32z = buf.drz[jj];
            const double r23x = delx, r23y = dely, r23z = delz;
            const double rij2 = rsq;
            for (int kk = 0; kk < jnum; ++kk) {
              if (kk == jj) continue;
              const double r21x = -buf.drx[kk], r21y = -buf.dry[kk], r21z = -buf.drz[kk];
              const double r21mag = std::sqrt(buf.d2[kk]);
              if (r21mag <= TOL) continue;
              const double cos321  = -(r21x*r32x+r21y*r32y+r21z*r32z)/(r21mag*rij);
              const double cos321c = std::min(1.0,std::max(-1.0,cos321));
              const double sin321  = std::sqrt(std::max(0.0,1.0-cos321c*cos321c));
              if (sin321 <= TOL) continue;
              double dcut321;
              const double rr_ik = 0.5*(rij2+buf.d2[kk]
                - (r21x-r23x)*(r21x-r23x)-(r21y-r23y)*(r21y-r23y)-(r21z-r23z)*(r21z-r23z))
                / (rij*r21mag);
              const double tspjik = rebo_Sp2(rr_ik, thmin, thmax, dcut321);
              double dw21;
              const double w21 = rebo_Sp(r21mag, p.rcmin[itype][buf.ext.type[kk]],
                                          p.rcmaxp[itype][buf.ext.type[kk]], dw21);
              for (int ll = 0; ll < jnum; ++ll) {
                if (ll == jj || ll == kk) continue;
                const double r34x = buf.drx[jj]-buf.drx[ll];
                const double r34y = buf.dry[jj]-buf.dry[ll];
                const double r34z = buf.drz[jj]-buf.drz[ll];
                const double r34mag = std::sqrt(r34x*r34x+r34y*r34y+r34z*r34z);
                if (r34mag <= TOL) continue;
                const double cos234  = (r32x*r34x+r32y*r34y+r32z*r34z)/(rij*r34mag);
                const double cos234c = std::min(1.0,std::max(-1.0,cos234));
                const double sin234  = std::sqrt(std::max(0.0,1.0-cos234c*cos234c));
                if (sin234 <= TOL) continue;
                double dcut234;
                const double rr_jl = 0.5*(rij2+r34mag*r34mag
                  - (r23x+r34x)*(r23x+r34x)-(r23y+r34y)*(r23y+r34y)-(r23z+r34z)*(r23z+r34z))
                  / (rij*r34mag);
                const double tspijl = rebo_Sp2(rr_jl, thmin, thmax, dcut234);
                double dw34;
                const double w34 = rebo_Sp(r34mag, p.rcmin[jtype][buf.ext.type[ll]],
                                            p.rcmaxp[jtype][buf.ext.type[ll]], dw34);
                const double cx321_0 = r32y*r21z-r32z*r21y;
                const double cx321_1 = r32z*r21x-r32x*r21z;
                const double cx321_2 = r32x*r21y-r32y*r21x;
                const double cx234_0 = r23y*r34z-r23z*r34y;
                const double cx234_1 = r23z*r34x-r23x*r34z;
                const double cx234_2 = r23x*r34y-r23y*r34x;
                const double cwnum = cx321_0*cx234_0+cx321_1*cx234_1+cx321_2*cx234_2;
                const double cwnom = r21mag*r34mag*rij*rij*sin321*sin234;
                const double om1234 = cwnum/cwnom;
                Etmp_Tij += (1.0-om1234*om1234)*w21*w34*(1.0-tspjik)*(1.0-tspijl);
              }
            }
          }

          const double bij = 0.5*(pij+pji) + piRC + Tij*Etmp_Tij;
          Stb  = rebo_Sp2(bij, p.bLJmin[itype][jtype], p.bLJmax[itype][jtype], dStb);

          // ---- Bond-order forces (LAMMPS convention: VA_lj = Str*cij*VLJ, no dStb) ----
          const double VA_lj = Str * cij * VLJ;
          if (std::fabs(VA_lj) > TOL) {

            // --- pij k-loop forces ---
            for (int kk = 0; kk < jnum; ++kk) {
              if (kk == jj) continue;
              const int    ktype = buf.ext.type[kk];
              const double rik   = std::sqrt(buf.d2[kk]);
              double dwik;
              const double wik = rebo_Sp(rik, p.rcmin[itype][ktype], p.rcmax[itype][ktype], dwik);
              if (wik <= TOL && dwik == 0.0) continue;

              const double lam = (itype==1)
                ? 4.0*((p.rho[ktype][1]-rik) - (p.rho[jtype][1]-rij_bo)) : 0.0;
              const double exp_lam = std::exp(lam);

              const double rikx = -buf.drx[kk], riky = -buf.dry[kk], rikz = -buf.drz[kk];

              double cosjik = (buf.drx[jj]*buf.drx[kk]+buf.dry[jj]*buf.dry[kk]
                               +buf.drz[jj]*buf.drz[kk]) / (rij*rik);
              cosjik = std::min(1.0, std::max(-1.0, cosjik));

              double dgdc, dgdN;
              const double g = rebo_gSpline(cosjik, Nij, itype, p, dgdc, dgdN);

              const double inv_rij = 1.0/rij, inv_rik = 1.0/rik;
              const double inv_rij2 = inv_rij*inv_rij, inv_rik2 = inv_rik*inv_rik;
              const double irr = inv_rij*inv_rik;

              const double dcdrix = (delx+rikx)*irr - cosjik*(delx*inv_rij2+rikx*inv_rik2);
              const double dcdriy = (dely+riky)*irr - cosjik*(dely*inv_rij2+riky*inv_rik2);
              const double dcdriz = (delz+rikz)*irr - cosjik*(delz*inv_rij2+rikz*inv_rik2);
              const double dcdrjx = -rikx*irr + cosjik*delx*inv_rij2;
              const double dcdrjy = -riky*irr + cosjik*dely*inv_rij2;
              const double dcdrjz = -rikz*irr + cosjik*delz*inv_rij2;
              const double dcdrkx = -delx*irr + cosjik*rikx*inv_rik2;
              const double dcdrky = -dely*irr + cosjik*riky*inv_rik2;
              const double dcdrkz = -delz*irr + cosjik*rikz*inv_rik2;

              const double pf_ang = VA_lj * 0.25 * tmp_pij * wik * dgdc * exp_lam;
              double fix = -pf_ang*dcdrix, fiy = -pf_ang*dcdriy, fiz = -pf_ang*dcdriz;
              double fjx = -pf_ang*dcdrjx, fjy = -pf_ang*dcdrjy, fjz = -pf_ang*dcdrjz;
              double fkx = -pf_ang*dcdrkx, fky = -pf_ang*dcdrky, fkz = -pf_ang*dcdrkz;

              // lambda forces: d(lam)/d(rij)=0 (rcmin is constant), so only rik terms
              if (itype == 1) {
                const double pf_lam = VA_lj * 0.25 * tmp_pij * wik * g * exp_lam * 4.0;
                fix -= pf_lam*(-rikx*inv_rik);
                fiy -= pf_lam*(-riky*inv_rik);
                fiz -= pf_lam*(-rikz*inv_rik);
                fkx -= pf_lam*(rikx*inv_rik);
                fky -= pf_lam*(riky*inv_rik);
                fkz -= pf_lam*(rikz*inv_rik);
              }

              // wik forces
              const double pf_wik = VA_lj * 0.25 * tmp_pij * dwik * g * exp_lam / rik;
              fix -= pf_wik*rikx; fiy -= pf_wik*riky; fiz -= pf_wik*rikz;
              fkx += pf_wik*rikx; fky += pf_wik*riky; fkz += pf_wik*rikz;

              // dNij forces (from tmp3_pij, via dwik)
              const double pf_dN = VA_lj * 0.25 * tmp_pij * tmp3_pij * dwik / rik;
              fix -= pf_dN*rikx; fiy -= pf_dN*riky; fiz -= pf_dN*rikz;
              fkx += pf_dN*rikx; fky += pf_dN*riky; fkz += pf_dN*rikz;

              // PijS forces
              const double pf_pij = VA_lj * 0.25 * tmp_pij * dN2_pij[ktype] * dwik / rik;
              fix -= pf_pij*rikx; fiy -= pf_pij*riky; fiz -= pf_pij*rikz;
              fkx += pf_pij*rikx; fky += pf_pij*riky; fkz += pf_pij*rikz;

              _fx += fix; _fy += fiy; _fz += fiz;
              atomic_add_contribution(cells[cell_j][field::fx][p_j], fjx);
              atomic_add_contribution(cells[cell_j][field::fy][p_j], fjy);
              atomic_add_contribution(cells[cell_j][field::fz][p_j], fjz);
              size_t cell_k, p_k;
              buf.nbh.get(kk, cell_k, p_k);
              atomic_add_contribution(cells[cell_k][field::fx][p_k], fkx);
              atomic_add_contribution(cells[cell_k][field::fy][p_k], fky);
              atomic_add_contribution(cells[cell_k][field::fz][p_k], fkz);
            }

            // --- pji l-loop forces ---
            for (int ll = 0; ll < jnum; ++ll) {
              if (ll == jj) continue;
              const int    ltype = buf.ext.type[ll];
              const double jlx = buf.drx[ll]-buf.drx[jj];
              const double jly = buf.dry[ll]-buf.dry[jj];
              const double jlz = buf.drz[ll]-buf.drz[jj];
              const double rjl = std::sqrt(jlx*jlx+jly*jly+jlz*jlz);
              double dwjl;
              const double wjl = rebo_Sp(rjl, p.rcmin[jtype][ltype], p.rcmax[jtype][ltype], dwjl);
              if (wjl <= TOL && dwjl == 0.0) continue;

              const double lam = (jtype==1)
                ? 4.0*((p.rho[ltype][1]-rjl) - (p.rho[itype][1]-rij_bo)) : 0.0;
              const double exp_lam = std::exp(lam);

              double cosijl = (delx*jlx+dely*jly+delz*jlz)/(rij*rjl);
              cosijl = std::min(1.0, std::max(-1.0, cosijl));

              double dgdc, dgdN;
              const double g = rebo_gSpline(cosijl, Nji, jtype, p, dgdc, dgdN);

              const double inv_rij = 1.0/rij, inv_rjl = 1.0/rjl;
              const double inv_rij2 = inv_rij*inv_rij, inv_rjl2 = inv_rjl*inv_rjl;
              const double irr = inv_rij*inv_rjl;

              const double dcix = jlx*irr - cosijl*delx*inv_rij2;
              const double dciy = jly*irr - cosijl*dely*inv_rij2;
              const double dciz = jlz*irr - cosijl*delz*inv_rij2;
              const double dcjx = (-delx-jlx)*irr + cosijl*(delx*inv_rij2+jlx*inv_rjl2);
              const double dcjy = (-dely-jly)*irr + cosijl*(dely*inv_rij2+jly*inv_rjl2);
              const double dcjz = (-delz-jlz)*irr + cosijl*(delz*inv_rij2+jlz*inv_rjl2);
              const double dclx = delx*irr - cosijl*jlx*inv_rjl2;
              const double dcly = dely*irr - cosijl*jly*inv_rjl2;
              const double dclz = delz*irr - cosijl*jlz*inv_rjl2;

              const double pf_ang = VA_lj * 0.25 * tmp_pji * wjl * dgdc * exp_lam;
              double fix = -pf_ang*dcix, fiy = -pf_ang*dciy, fiz = -pf_ang*dciz;
              double fjx = -pf_ang*dcjx, fjy = -pf_ang*dcjy, fjz = -pf_ang*dcjz;
              double flx = -pf_ang*dclx, fly = -pf_ang*dcly, flz = -pf_ang*dclz;

              // lambda forces: d(lam)/d(rij)=0, only jl terms
              if (jtype == 1) {
                const double pf_lam = VA_lj * 0.25 * tmp_pji * wjl * g * exp_lam * 4.0;
                fjx -= pf_lam*(jlx*inv_rjl);
                fjy -= pf_lam*(jly*inv_rjl);
                fjz -= pf_lam*(jlz*inv_rjl);
                flx -= pf_lam*(-jlx*inv_rjl);
                fly -= pf_lam*(-jly*inv_rjl);
                flz -= pf_lam*(-jlz*inv_rjl);
              }

              const double pf_wjl = VA_lj * 0.25 * tmp_pji * dwjl * g * exp_lam / rjl;
              fjx -= pf_wjl*jlx; fjy -= pf_wjl*jly; fjz -= pf_wjl*jlz;
              flx += pf_wjl*jlx; fly += pf_wjl*jly; flz += pf_wjl*jlz;

              const double pf_dN = VA_lj * 0.25 * tmp_pji * tmp3_pji * dwjl / rjl;
              fjx -= pf_dN*jlx; fjy -= pf_dN*jly; fjz -= pf_dN*jlz;
              flx += pf_dN*jlx; fly += pf_dN*jly; flz += pf_dN*jlz;

              const double pf_pji = VA_lj * 0.25 * tmp_pji * dN2_pji[ltype] * dwjl / rjl;
              fjx -= pf_pji*jlx; fjy -= pf_pji*jly; fjz -= pf_pji*jlz;
              flx += pf_pji*jlx; fly += pf_pji*jly; flz += pf_pji*jlz;

              _fx += fix; _fy += fiy; _fz += fiz;
              atomic_add_contribution(cells[cell_j][field::fx][p_j], fjx);
              atomic_add_contribution(cells[cell_j][field::fy][p_j], fjy);
              atomic_add_contribution(cells[cell_j][field::fz][p_j], fjz);
              size_t cell_l, p_l;
              buf.nbh.get(ll, cell_l, p_l);
              atomic_add_contribution(cells[cell_l][field::fx][p_l], flx);
              atomic_add_contribution(cells[cell_l][field::fy][p_l], fly);
              atomic_add_contribution(cells[cell_l][field::fz][p_l], flz);
            }

            // --- piRC k-loop forces (dN3[0] and dN3[2]) ---
            for (int kk = 0; kk < jnum; ++kk) {
              if (kk == jj) continue;
              const int    ktype = buf.ext.type[kk];
              const double rik   = std::sqrt(buf.d2[kk]);
              double dwik;
              const double wik = rebo_Sp(rik, p.rcmin[itype][ktype], p.rcmax[itype][ktype], dwik);
              const double rikx = -buf.drx[kk], riky = -buf.dry[kk], rikz = -buf.drz[kk];
              size_t cell_k, p_k;
              buf.nbh.get(kk, cell_k, p_k);

              const double pf0 = 0.5 * VA_lj * dN3_pi[0] * dwik / rik;
              _fx -= pf0*rikx; _fy -= pf0*riky; _fz -= pf0*rikz;
              atomic_add_contribution(cells[cell_k][field::fx][p_k], pf0*rikx);
              atomic_add_contribution(cells[cell_k][field::fy][p_k], pf0*riky);
              atomic_add_contribution(cells[cell_k][field::fz][p_k], pf0*rikz);

              if (ktype == 0) {
                const double Nki = cells[cell_k][m_nijc_field][p_k] + cells[cell_k][m_nijh_field][p_k] - wik;
                double dNki;
                const double SpNki = rebo_Sp(Nki, p.Nmin, p.Nmax, dNki);
                const double pf2a = 0.5*VA_lj*dN3_pi[2]*2.0*NconjtmpI*dwik*SpNki/rik;
                _fx -= pf2a*rikx; _fy -= pf2a*riky; _fz -= pf2a*rikz;
                atomic_add_contribution(cells[cell_k][field::fx][p_k], pf2a*rikx);
                atomic_add_contribution(cells[cell_k][field::fy][p_k], pf2a*riky);
                atomic_add_contribution(cells[cell_k][field::fz][p_k], pf2a*rikz);

                if (std::fabs(dNki) > TOL) {
                  for (int mm = 0; mm < jnum; ++mm) {
                    if (mm == kk) continue;
                    const int mtype = buf.ext.type[mm];
                    const double kmx = buf.drx[mm]-buf.drx[kk];
                    const double kmy = buf.dry[mm]-buf.dry[kk];
                    const double kmz = buf.drz[mm]-buf.drz[kk];
                    const double rkm = std::sqrt(kmx*kmx+kmy*kmy+kmz*kmz);
                    if (rkm <= TOL) continue;
                    double dwkm;
                    rebo_Sp(rkm, p.rcmin[ktype][mtype], p.rcmax[ktype][mtype], dwkm);
                    const double pf2b = 0.5*VA_lj*dN3_pi[2]*2.0*NconjtmpI*wik*dNki*dwkm/rkm;
                    size_t cell_m, p_m;
                    buf.nbh.get(mm, cell_m, p_m);
                    atomic_add_contribution(cells[cell_k][field::fx][p_k], -pf2b*kmx);
                    atomic_add_contribution(cells[cell_k][field::fy][p_k], -pf2b*kmy);
                    atomic_add_contribution(cells[cell_k][field::fz][p_k], -pf2b*kmz);
                    atomic_add_contribution(cells[cell_m][field::fx][p_m],  pf2b*kmx);
                    atomic_add_contribution(cells[cell_m][field::fy][p_m],  pf2b*kmy);
                    atomic_add_contribution(cells[cell_m][field::fz][p_m],  pf2b*kmz);
                  }
                }
              }
            }

            // --- piRC l-loop forces (dN3[1] and dN3[2]) ---
            for (int ll = 0; ll < jnum; ++ll) {
              if (ll == jj) continue;
              const int    ltype = buf.ext.type[ll];
              const double jlx = buf.drx[ll]-buf.drx[jj];
              const double jly = buf.dry[ll]-buf.dry[jj];
              const double jlz = buf.drz[ll]-buf.drz[jj];
              const double rjl = std::sqrt(jlx*jlx+jly*jly+jlz*jlz);
              double dwjl;
              const double wjl = rebo_Sp(rjl, p.rcmin[jtype][ltype], p.rcmax[jtype][ltype], dwjl);
              size_t cell_l, p_l;
              buf.nbh.get(ll, cell_l, p_l);

              const double pf1 = 0.5*VA_lj*dN3_pi[1]*dwjl/rjl;
              atomic_add_contribution(cells[cell_j][field::fx][p_j], -pf1*jlx);
              atomic_add_contribution(cells[cell_j][field::fy][p_j], -pf1*jly);
              atomic_add_contribution(cells[cell_j][field::fz][p_j], -pf1*jlz);
              atomic_add_contribution(cells[cell_l][field::fx][p_l],  pf1*jlx);
              atomic_add_contribution(cells[cell_l][field::fy][p_l],  pf1*jly);
              atomic_add_contribution(cells[cell_l][field::fz][p_l],  pf1*jlz);

              if (ltype == 0) {
                const double Nlj = cells[cell_l][m_nijc_field][p_l] + cells[cell_l][m_nijh_field][p_l] - wjl;
                double dNlj;
                const double SpNlj = rebo_Sp(Nlj, p.Nmin, p.Nmax, dNlj);
                const double pf2a = 0.5*VA_lj*dN3_pi[2]*2.0*NconjtmpJ*dwjl*SpNlj/rjl;
                atomic_add_contribution(cells[cell_j][field::fx][p_j], -pf2a*jlx);
                atomic_add_contribution(cells[cell_j][field::fy][p_j], -pf2a*jly);
                atomic_add_contribution(cells[cell_j][field::fz][p_j], -pf2a*jlz);
                atomic_add_contribution(cells[cell_l][field::fx][p_l],  pf2a*jlx);
                atomic_add_contribution(cells[cell_l][field::fy][p_l],  pf2a*jly);
                atomic_add_contribution(cells[cell_l][field::fz][p_l],  pf2a*jlz);

                if (std::fabs(dNlj) > TOL) {
                  for (int mm = 0; mm < jnum; ++mm) {
                    if (mm == ll || mm == jj) continue;
                    const int mtype = buf.ext.type[mm];
                    const double lmx = buf.drx[mm]-buf.drx[ll];
                    const double lmy = buf.dry[mm]-buf.dry[ll];
                    const double lmz = buf.drz[mm]-buf.drz[ll];
                    const double rlm = std::sqrt(lmx*lmx+lmy*lmy+lmz*lmz);
                    if (rlm <= TOL) continue;
                    double dwlm;
                    rebo_Sp(rlm, p.rcmin[ltype][mtype], p.rcmax[ltype][mtype], dwlm);
                    const double pf2b = 0.5*VA_lj*dN3_pi[2]*2.0*NconjtmpJ*wjl*dNlj*dwlm/rlm;
                    size_t cell_m, p_m;
                    buf.nbh.get(mm, cell_m, p_m);
                    atomic_add_contribution(cells[cell_l][field::fx][p_l], -pf2b*lmx);
                    atomic_add_contribution(cells[cell_l][field::fy][p_l], -pf2b*lmy);
                    atomic_add_contribution(cells[cell_l][field::fz][p_l], -pf2b*lmz);
                    atomic_add_contribution(cells[cell_m][field::fx][p_m],  pf2b*lmx);
                    atomic_add_contribution(cells[cell_m][field::fy][p_m],  pf2b*lmy);
                    atomic_add_contribution(cells[cell_m][field::fz][p_m],  pf2b*lmz);
                  }
                }
              }
            }

          } // end if (|VA_lj| > TOL)
        } // end Str > TOL block

        // ===================================================================
        // Step 4: Direct pair force and energy
        // ===================================================================
        const double fpair = -(dStr*(Stb-1.0)*cij*VLJ
                               + dVLJ*(Str*Stb + 1.0 - Str)*cij) / rij;

        _en += 0.5*(Str*Stb + 1.0-Str)*cij*VLJ;

        _fx += 0.5*delx*fpair; _fy += 0.5*dely*fpair; _fz += 0.5*delz*fpair;
        atomic_add_contribution(cells[cell_j][field::fx][p_j], -0.5*delx*fpair);
        atomic_add_contribution(cells[cell_j][field::fy][p_j], -0.5*dely*fpair);
        atomic_add_contribution(cells[cell_j][field::fz][p_j], -0.5*delz*fpair);

        // ===================================================================
        // Step 5: Path forces from dcij (×0.5 full-list)
        // ===================================================================
        if (cij < 1.0) {
          const double dC = VLJ * (Str*Stb + 1.0 - Str);

          if (best_npath == 2) {
            const double fp = 0.5 * dC * best_dwij2 / rij;
            _fx += delx*fp; _fy += dely*fp; _fz += delz*fp;
            atomic_add_contribution(cells[cell_j][field::fx][p_j], -delx*fp);
            atomic_add_contribution(cells[cell_j][field::fy][p_j], -dely*fp);
            atomic_add_contribution(cells[cell_j][field::fz][p_j], -delz*fp);

          } else if (best_npath == 3) {
            const double fp1 = 0.5 * dC * best_dwik3 * best_wkj3 / best_rik3;
            const double fp2 = 0.5 * dC * best_wik3 * best_dwkj3 / best_rkj3;
            const double fix = best_delikx3*fp1, fiy = best_deliky3*fp1, fiz = best_delikz3*fp1;
            const double fjx = best_deljkx3*fp2, fjy = best_deljky3*fp2, fjz = best_deljkz3*fp2;
            _fx += fix; _fy += fiy; _fz += fiz;
            atomic_add_contribution(cells[cell_j][field::fx][p_j],  fjx);
            atomic_add_contribution(cells[cell_j][field::fy][p_j],  fjy);
            atomic_add_contribution(cells[cell_j][field::fz][p_j],  fjz);
            size_t cell_k, p_k;
            buf.nbh.get(best_kk3, cell_k, p_k);
            atomic_add_contribution(cells[cell_k][field::fx][p_k], -(fix+fjx));
            atomic_add_contribution(cells[cell_k][field::fy][p_k], -(fiy+fjy));
            atomic_add_contribution(cells[cell_k][field::fz][p_k], -(fiz+fjz));

          } else if (best_npath == 4) {
            const double fp1 = 0.5 * dC * best_dwik4 * best_wkm4 * best_wmj4 / best_rik4;
            const double fp2 = 0.5 * dC * best_wik4 * best_dwkm4 * best_wmj4 / best_rkm4;
            const double fp3 = 0.5 * dC * best_wik4 * best_wkm4 * best_dwmj4 / best_rmj4;
            const double fix = best_delikx4*fp1, fiy = best_deliky4*fp1, fiz = best_delikz4*fp1;
            const double fjx = best_deljmx4*fp3, fjy = best_deljmy4*fp3, fjz = best_deljmz4*fp3;
            const double fkx = best_delkmx4*fp2 - fix;
            const double fky = best_delkmy4*fp2 - fiy;
            const double fkz = best_delkmz4*fp2 - fiz;
            const double fmx = -best_delkmx4*fp2 - fjx;
            const double fmy = -best_delkmy4*fp2 - fjy;
            const double fmz = -best_delkmz4*fp2 - fjz;
            _fx += fix; _fy += fiy; _fz += fiz;
            atomic_add_contribution(cells[cell_j][field::fx][p_j],  fjx);
            atomic_add_contribution(cells[cell_j][field::fy][p_j],  fjy);
            atomic_add_contribution(cells[cell_j][field::fz][p_j],  fjz);
            size_t cell_k, p_k;
            buf.nbh.get(best_kk4, cell_k, p_k);
            atomic_add_contribution(cells[cell_k][field::fx][p_k], fkx);
            atomic_add_contribution(cells[cell_k][field::fy][p_k], fky);
            atomic_add_contribution(cells[cell_k][field::fz][p_k], fkz);
            size_t cell_m, p_m;
            buf.nbh.get(best_mm4, cell_m, p_m);
            atomic_add_contribution(cells[cell_m][field::fx][p_m], fmx);
            atomic_add_contribution(cells[cell_m][field::fy][p_m], fmy);
            atomic_add_contribution(cells[cell_m][field::fz][p_m], fmz);
          }
        }

      } // end j-loop

      atomic_add_contribution(en, _en);
      atomic_add_contribution(fx, _fx);
      atomic_add_contribution(fy, _fy);
      atomic_add_contribution(fz, _fz);
    }
  };

} // namespace exaStamp
