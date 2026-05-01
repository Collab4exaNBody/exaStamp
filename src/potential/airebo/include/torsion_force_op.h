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
#include "airebo_force_op_cached.h"   // shared buffer, helpers (rebo_Sp, rebo_Sp2, ...)

namespace exaStamp
{
  using namespace exanb;

  // =========================================================================
  // TorsionForceOp — explicit torsion barrier of AIREBO  (Stuart et al. 2000)
  //
  // Computes for every dihedral k-i-j-l where i-j is a C-C REBO bond:
  //   E = epsilonT[ktype][ltype] * (256/405 * cw2^5 - 1/10) * w_ij * w_ik * w_jl
  //         * (1 - tsp_jik) * (1 - tsp_ijl)
  // where cw2 = 0.5*(1 - cos(omega_{kijl})).
  //
  // Reference: Stuart, Tutein, Harrison, J. Chem. Phys. 112 (2000) eq. 12.
  // Port of LAMMPS PairAIREBO::TORSION() to exaStamp pairwise framework.
  //
  // Conventions (matching rebo_force_op.h):
  //   buf.drx[jj] = r_j - r_i  (exaStamp: from i to j)
  //   del32 = r_j - r_i  (LAMMPS: del32)
  //   del23 = r_i - r_j  (LAMMPS: del23)
  //   del21 = r_i - r_k  (LAMMPS: del21)
  //   del34 = r_j - r_l  (LAMMPS: del34)
  //
  // Full-list: all energy and force contributions × 0.5.
  // Only C-C central bonds (itype==0, jtype==0); k and l can be C or H.
  // =========================================================================

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) TorsionForceOp
  {
    const AireboParamsRO* m_params = nullptr;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int type, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int type, Mat3d& virial, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t n, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int type, CellParticlesT cells,
                           GridCellLocksT locks, ParticleLockT& lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class Mat3dT,
             class GridCellLocksT, class ParticleLockT>
    inline void operator()(int jnum, ComputeBufferT& buf, double& en,
                           double& fx, double& fy, double& fz,
                           int itype, Mat3dT& virial, CellParticlesT cells,
                           GridCellLocksT locks, ParticleLockT& lock_a) const
    {
      // TORSION only applies when i is a C atom (type 0)
      if (itype != 0) return;

      const AireboParamsRO& p = *m_params;
      static constexpr double conv  = 1.0; //EXASTAMP_CONST_QUANTITY(1.*eV);
      static constexpr double TOL   = 1.0e-9;
      static constexpr double thmin = -1.0;
      static constexpr double thmax = -0.995;

      for (int jj = 0; jj < jnum; ++jj) {
        const int jtype = buf.ext.type[jj];
        if (jtype != 0) continue;   // only C-C central bonds

        // ---- i-j bond vectors (LAMMPS: del32 = rj-ri, del23 = ri-rj) ----
        const double del32x = buf.drx[jj], del32y = buf.dry[jj], del32z = buf.drz[jj];
        const double r32sq  = buf.d2[jj];
        const double r23    = std::sqrt(r32sq);          // r_ij = r23 = r32
        const double del23x = -del32x, del23y = -del32y, del23z = -del32z;

        double dw23;
        const double w23 = rebo_Sp(r23, p.rcmin[0][0], p.rcmax[0][0], dw23);
        if (w23 <= TOL && dw23 == 0.0) continue;

        size_t cell_j, p_j;
        buf.nbh.get(jj, cell_j, p_j);

        // ---- loop over k: REBO neighbour of i, k != j ----
        for (int kk = 0; kk < jnum; ++kk) {
          if (kk == jj) continue;
          const int ktype = buf.ext.type[kk];

          // del21 = ri - rk  (buf.drx[kk] = rk - ri)
          const double del21x = -buf.drx[kk], del21y = -buf.dry[kk], del21z = -buf.drz[kk];
          const double r21    = std::sqrt(buf.d2[kk]);
          if (r21 <= TOL) continue;

          // cos of angle k-i-j at vertex i (= cos321 in LAMMPS)
          const double cos321 = std::min(1.0, std::max(-1.0,
              -(del21x*del32x + del21y*del32y + del21z*del32z) / (r21*r23)));
          const double sin321 = std::sqrt(std::max(0.0, 1.0 - cos321*cos321));
          if (sin321 < TOL) continue;

          // deljk = rj - rk = del21 - del23
          const double deljkx = del21x - del23x;
          const double deljky = del21y - del23y;
          const double deljkz = del21z - del23z;
          const double rjk    = std::sqrt(deljkx*deljkx + deljky*deljky + deljkz*deljkz);

          // tspjik: planarity switch for k-i-j angle
          const double rij2   = r32sq;
          const double rik2   = buf.d2[kk];
          const double costmp_jik = 0.5*(rij2 + rik2 - rjk*rjk) / (r23*r21);
          double dtsjik;
          const double tspjik = rebo_Sp2(costmp_jik, thmin, thmax, dtsjik);
          dtsjik = -dtsjik;   // LAMMPS sign convention

          double dw21;
          const double w21 = rebo_Sp(r21, p.rcmin[0][ktype], p.rcmax[0][ktype], dw21);

          // ---- loop over l: REBO neighbour of j, l != i, l != k ----
          for (int ll = 0; ll < jnum; ++ll) {
            if (ll == jj || ll == kk) continue;
            const int ltype = buf.ext.type[ll];

            // del34 = rj - rl = del32 - (rl - ri) = del32 - buf.drx[ll]
            const double del34x = del32x - buf.drx[ll];
            const double del34y = del32y - buf.dry[ll];
            const double del34z = del32z - buf.drz[ll];
            const double r34sq  = del34x*del34x + del34y*del34y + del34z*del34z;
            const double r34    = std::sqrt(r34sq);
            if (r34 <= TOL) continue;

            // cos of angle i-j-l at vertex j (= cos234 in LAMMPS)
            const double cos234 = std::min(1.0, std::max(-1.0,
                (del32x*del34x + del32y*del34y + del32z*del34z) / (r23*r34)));
            const double sin234 = std::sqrt(std::max(0.0, 1.0 - cos234*cos234));
            if (sin234 < TOL) continue;

            double dw34;
            const double w34 = rebo_Sp(r34, p.rcmin[0][ltype], p.rcmax[0][ltype], dw34);

            // delil = ri - rl = del23 + del34
            const double delilx = del23x + del34x;
            const double delily = del23y + del34y;
            const double delilz = del23z + del34z;
            const double ril2   = delilx*delilx + delily*delily + delilz*delilz;
            const double ril    = std::sqrt(ril2);

            // tspijl: planarity switch for i-j-l angle
            const double rjl2   = r34sq;
            const double costmp_ijl = 0.5*(rij2 + rjl2 - ril2) / (r23*r34);
            double dtsijl;
            const double tspijl = rebo_Sp2(costmp_ijl, thmin, thmax, dtsijl);
            dtsijl = -dtsijl;

            // ---- dihedral angle ----
            // cross321 = del32 × del21
            const double cx321_0 = del32y*del21z - del32z*del21y;
            const double cx321_1 = del32z*del21x - del32x*del21z;
            const double cx321_2 = del32x*del21y - del32y*del21x;
            const double cross321mag = std::sqrt(cx321_0*cx321_0 + cx321_1*cx321_1 + cx321_2*cx321_2);

            // cross234 = del23 × del34
            const double cx234_0 = del23y*del34z - del23z*del34y;
            const double cx234_1 = del23z*del34x - del23x*del34z;
            const double cx234_2 = del23x*del34y - del23y*del34x;
            const double cross234mag = std::sqrt(cx234_0*cx234_0 + cx234_1*cx234_1 + cx234_2*cx234_2);

            const double cwnum = cx321_0*cx234_0 + cx321_1*cx234_1 + cx321_2*cx234_2;
            const double cwnom = r21*r34*r23*r23*sin321*sin234;
            const double cw    = cwnum / cwnom;

            // ---- torsion energy ----
            const double cw2    = 0.5*(1.0 - cw);
            const double cw2p2  = cw2*cw2;
            const double cw2p4  = cw2p2*cw2p2;
            const double cw2p5  = cw2p4*cw2;
            const double ekijl  = conv * p.epsilonT[ktype][ltype];
            const double Ec     = 256.0*ekijl/405.0;
            const double Vtors  = Ec*cw2p5 - ekijl/10.0;
            const double w_prod = w21*w23*w34;
            const double tsp_prod = (1.0-tspjik)*(1.0-tspijl);

            en += 0.5 * Vtors * w_prod * tsp_prod;

            // ---- force prefactor: dvpdcw = -(dVtors/dcw) * w21*w23*w34*tsp_prod ----
            //   dVtors/dcw = Ec*5*cw2^4*(-0.5) = -2.5*Ec*cw2^4
            //   dvpdcw = 2.5*Ec*cw2^4 * w_prod * tsp_prod
            const double dvpdcw = 2.5*Ec*cw2p4 * w_prod * tsp_prod;
            const double dcwdn  =  1.0/cwnom;
            const double dcwddn = -cwnum/(cwnom*cwnom);

            // ---- d(cwnum)/d(del23), d(cwnum)/d(del21), d(cwnum)/d(del34) ----
            // (identical to LAMMPS dndij, dndik, dndjl)
            const double dndij_0 = cx234_1*del21z - cx234_2*del21y + del34y*cx321_2 - del34z*cx321_1;
            const double dndij_1 = cx234_2*del21x - cx234_0*del21z + del34z*cx321_0 - del34x*cx321_2;
            const double dndij_2 = cx234_0*del21y - cx234_1*del21x + del34x*cx321_1 - del34y*cx321_0;

            const double dndik_0 = del23y*cx234_2 - del23z*cx234_1;
            const double dndik_1 = del23z*cx234_0 - del23x*cx234_2;
            const double dndik_2 = del23x*cx234_1 - del23y*cx234_0;

            const double dndjl_0 = cx321_1*del23z - cx321_2*del23y;
            const double dndjl_1 = cx321_2*del23x - cx321_0*del23z;
            const double dndjl_2 = cx321_0*del23y - cx321_1*del23x;

            // ---- d(cwnom)/d(|bondlength|): ddndij, ddndik, ddndjk, ddndjl, ddndil ----
            // cwnom = (r21*r23*sin321) * (r34*r23*sin234)
            // Group xi = r21*r23*sin321 (cross321mag), xj = r34*r23*sin234 (cross234mag)
            const double dcidij = (r23*r23 - r21*r21 + rjk*rjk) / (2.0*r23*r23*r21);
            const double dcidik = (r21*r21 - r23*r23 + rjk*rjk) / (2.0*r23*r21*r21);
            const double dcidjk = -rjk / (r23*r21);

            const double dsidij = (-cos321/sin321)*dcidij;
            const double dsidik = (-cos321/sin321)*dcidik;
            const double dsidjk = (-cos321/sin321)*dcidjk;

            // d(r21*r23*sin321)/d(r23), d(r21)/d(r21), d(rjk)
            const double dxidij = r21*sin321 + r23*r21*dsidij;
            const double dxidik = r23*sin321 + r23*r21*dsidik;
            const double dxidjk = r23*r21*dsidjk;

            const double dcjdji = (r23*r23 - r34*r34 + ril2) / (2.0*r23*r23*r34);
            const double dcjdjl = (r34*r34 - r23*r23 + ril2) / (2.0*r23*r34*r34);
            const double dcjdil = -ril / (r23*r34);

            const double dsjdji = (-cos234/sin234)*dcjdji;
            const double dsjdjl = (-cos234/sin234)*dcjdjl;
            const double dsjdil = (-cos234/sin234)*dcjdil;

            // d(r34*r23*sin234)/d(r23), d(r34), d(ril)
            const double dxjdji = r34*sin234 + r23*r34*dsjdji;
            const double dxjdjl = r23*sin234 + r23*r34*dsjdjl;
            const double dxjdil = r23*r34*dsjdil;

            const double ddndij = dxidij*cross234mag + cross321mag*dxjdji;
            const double ddndik = dxidik*cross234mag;
            const double ddndjk = dxidjk*cross234mag;
            const double ddndjl = cross321mag*dxjdjl;
            const double ddndil = cross321mag*dxjdil;

            // ---- accumulate forces on i, j, k, l ----
            double fi0 = 0, fi1 = 0, fi2 = 0;
            double fj0 = 0, fj1 = 0, fj2 = 0;
            double fk0 = 0, fk1 = 0, fk2 = 0;
            double fl0 = 0, fl1 = 0, fl2 = 0;

            // -- cw-variation forces (dihedral geometry) --
            // del23 direction: affects i (+) and j (-)
            double Ftmp;
            Ftmp = dvpdcw*(dcwdn*dndij_0 + dcwddn*ddndij*del23x/r23);
            fi0 += Ftmp; fj0 -= Ftmp;
            Ftmp = dvpdcw*(dcwdn*dndij_1 + dcwddn*ddndij*del23y/r23);
            fi1 += Ftmp; fj1 -= Ftmp;
            Ftmp = dvpdcw*(dcwdn*dndij_2 + dcwddn*ddndij*del23z/r23);
            fi2 += Ftmp; fj2 -= Ftmp;

            // del21 direction: affects i (+) and k (-)
            Ftmp = dvpdcw*(dcwdn*dndik_0 + dcwddn*ddndik*del21x/r21);
            fi0 += Ftmp; fk0 -= Ftmp;
            Ftmp = dvpdcw*(dcwdn*dndik_1 + dcwddn*ddndik*del21y/r21);
            fi1 += Ftmp; fk1 -= Ftmp;
            Ftmp = dvpdcw*(dcwdn*dndik_2 + dcwddn*ddndik*del21z/r21);
            fi2 += Ftmp; fk2 -= Ftmp;

            // deljk direction: affects j (+) and k (-)
            Ftmp = dvpdcw*dcwddn*ddndjk*deljkx/rjk;
            fj0 += Ftmp; fk0 -= Ftmp;
            Ftmp = dvpdcw*dcwddn*ddndjk*deljky/rjk;
            fj1 += Ftmp; fk1 -= Ftmp;
            Ftmp = dvpdcw*dcwddn*ddndjk*deljkz/rjk;
            fj2 += Ftmp; fk2 -= Ftmp;

            // del34 direction: affects j (+) and l (-)
            Ftmp = dvpdcw*(dcwdn*dndjl_0 + dcwddn*ddndjl*del34x/r34);
            fj0 += Ftmp; fl0 -= Ftmp;
            Ftmp = dvpdcw*(dcwdn*dndjl_1 + dcwddn*ddndjl*del34y/r34);
            fj1 += Ftmp; fl1 -= Ftmp;
            Ftmp = dvpdcw*(dcwdn*dndjl_2 + dcwddn*ddndjl*del34z/r34);
            fj2 += Ftmp; fl2 -= Ftmp;

            // delil direction: affects i (+) and l (-)
            Ftmp = dvpdcw*dcwddn*ddndil*delilx/ril;
            fi0 += Ftmp; fl0 -= Ftmp;
            Ftmp = dvpdcw*dcwddn*ddndil*delily/ril;
            fi1 += Ftmp; fl1 -= Ftmp;
            Ftmp = dvpdcw*dcwddn*ddndil*delilz/ril;
            fi2 += Ftmp; fl2 -= Ftmp;

            // -- coordination forces (w derivatives) --
            double fpair;
            fpair = Vtors * dw21 * w23 * w34 * tsp_prod / r21;
            fi0 -= del21x*fpair; fi1 -= del21y*fpair; fi2 -= del21z*fpair;
            fk0 += del21x*fpair; fk1 += del21y*fpair; fk2 += del21z*fpair;

            fpair = Vtors * w21 * dw23 * w34 * tsp_prod / r23;
            fi0 -= del23x*fpair; fi1 -= del23y*fpair; fi2 -= del23z*fpair;
            fj0 += del23x*fpair; fj1 += del23y*fpair; fj2 += del23z*fpair;

            fpair = Vtors * w21 * w23 * dw34 * tsp_prod / r34;
            fj0 -= del34x*fpair; fj1 -= del34y*fpair; fj2 -= del34z*fpair;
            fl0 += del34x*fpair; fl1 += del34y*fpair; fl2 += del34z*fpair;

            // -- tsp_jik forces (angle-cutoff derivatives) --
            // fcpc = -Vtors*w21*w23*w34*dtsjik*(1-tspijl)  [dtsjik already negated]
            const double fcpc_jik = -Vtors*w21*w23*w34*dtsjik*(1.0-tspijl);

            fpair = fcpc_jik*dcidij/r23;
            fi0 += fpair*del23x; fi1 += fpair*del23y; fi2 += fpair*del23z;
            fj0 -= fpair*del23x; fj1 -= fpair*del23y; fj2 -= fpair*del23z;

            fpair = fcpc_jik*dcidik/r21;
            fi0 += fpair*del21x; fi1 += fpair*del21y; fi2 += fpair*del21z;
            fk0 -= fpair*del21x; fk1 -= fpair*del21y; fk2 -= fpair*del21z;

            fpair = fcpc_jik*dcidjk/rjk;
            fj0 += fpair*deljkx; fj1 += fpair*deljky; fj2 += fpair*deljkz;
            fk0 -= fpair*deljkx; fk1 -= fpair*deljky; fk2 -= fpair*deljkz;

            // -- tsp_ijl forces (angle-cutoff derivatives) --
            const double fcpc_ijl = -Vtors*w21*w23*w34*(1.0-tspjik)*dtsijl;

            fpair = fcpc_ijl*dcjdji/r23;
            fi0 += fpair*del23x; fi1 += fpair*del23y; fi2 += fpair*del23z;
            fj0 -= fpair*del23x; fj1 -= fpair*del23y; fj2 -= fpair*del23z;

            fpair = fcpc_ijl*dcjdjl/r34;
            fj0 += fpair*del34x; fj1 += fpair*del34y; fj2 += fpair*del34z;
            fl0 -= fpair*del34x; fl1 -= fpair*del34y; fl2 -= fpair*del34z;

            fpair = fcpc_ijl*dcjdil/ril;
            fi0 += fpair*delilx; fi1 += fpair*delily; fi2 += fpair*delilz;
            fl0 -= fpair*delilx; fl1 -= fpair*delily; fl2 -= fpair*delilz;

            // ---- scatter forces (full-list × 0.5) ----
            fx += 0.5*fi0; fy += 0.5*fi1; fz += 0.5*fi2;

            size_t cell_k, p_k;
            buf.nbh.get(kk, cell_k, p_k);
            size_t cell_l, p_l;
            buf.nbh.get(ll, cell_l, p_l);

            atomic_add_contribution(cells[cell_j][field::fx][p_j], 0.5*fj0);
            atomic_add_contribution(cells[cell_j][field::fy][p_j], 0.5*fj1);
            atomic_add_contribution(cells[cell_j][field::fz][p_j], 0.5*fj2);
            atomic_add_contribution(cells[cell_k][field::fx][p_k], 0.5*fk0);
            atomic_add_contribution(cells[cell_k][field::fy][p_k], 0.5*fk1);
            atomic_add_contribution(cells[cell_k][field::fz][p_k], 0.5*fk2);
            atomic_add_contribution(cells[cell_l][field::fx][p_l], 0.5*fl0);
            atomic_add_contribution(cells[cell_l][field::fy][p_l], 0.5*fl1);
            atomic_add_contribution(cells[cell_l][field::fz][p_l], 0.5*fl2);
          } // ll (l)
        } // kk (k)
      } // jj (j)
    }
  };

} // namespace exaStamp
