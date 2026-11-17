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

// Flat tables of a linear (YE3T / PACE-compatible .yace) atomic cluster expansion, built by ye3t_init, read by the
// kernels of ye3t_kernels.h through Ye3tView (raw pointers into device-addressable vectors).
//
// Model, for a central atom i of element mu (ye3t-lammps pair_style ye3t, .yace family):
//   radial base g_0 = fc, g_k = (1 - T_k(x)) / 2 fc (k >= 1) with x = 1 - 2 (exp(-lambda (r/rc - 1)) - 1) / (exp(lambda) - 1),
//     fc = (1 + cos(pi r / rc)) / 2, times (1 + cos(pi (r - rc + dcut) / dcut)) / 2 within dcut of rc  (ChebExpCos)
//   contracted radials R_nl = sum_k c_nlk g_k
//   g_k tabulated as uniform cubic Hermite splines, h = rc / n, n = floor(rc / deltaSplineBins), nodes (values and
//     derivatives) at r = h, 2h .. n h = rc, interval k = floor(r / h) >= 1; below h they are evaluated directly; R_nl
//     contracted from the interpolated g_k (the spline being linear in its nodes, this is the spline of R_nl of PACE)
//   sources (m >= 0): A_{mu_j k} = sum_j g_k(r_ij) (rank 1), A_{mu_j n l m} = sum_j R_nl(r_ij) Y_lm(r_ij / r) with Y_lm the
//     complex Condon-Shortley harmonics normalised to Y_00 = 1
//   channels: a source, or (-1)^m conj of a source for -m
//   rho = sum over monomials c prod_f channel_f ^ p_f   (the ctilde rows of the functions, merged; real coefficients
//     symmetric under m -> -m, so that Re rho uses one monomial of each conjugate pair, products shared along a prefix tree)
//   E_i = E0[mu] + FS0[mu] Re rho

#include <onika/cuda/cuda.h>
#include <onika/memory/allocator.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace exaStamp
{
  namespace mlip_ye3t
  {
    static constexpr int kYe3tMaxLmax = 12;       // largest angular degree of a model (templated kernels up to this)
    static constexpr int kYe3tMaxRank = 32;       // largest rank (factors, powers expanded) of a monomial

    // one directed bond (mu_i, mu_j), offsets into the double pool
    struct Ye3tBond
    {
      double rc = 0.0, dcut = 0.0, lambda = 1.0;
      double scale = 0.0;         // n / rc
      int nint = 0;               // spline intervals n
      int nbase = 0;              // radial base functions g_k
      int nrad = 0;               // radial channels n of R_nl
      int lb = 0;                 // lmax of the bond
      int ncontr = 0;             // nrad (lb + 1), R_nl at n (lb + 1) + l
      int64_t bspl = 0;           // base spline coefficients [(nint+1)][nbase][4]
      int64_t crad = 0;           // contraction coefficients [ncontr][nbase] (direct evaluation below the table)
    };

    struct Ye3tTables
    {
      std::vector<std::string> elements;
      int nel = 0;
      int lmax = 0;               // largest lmax of the bonds
      int nlh = 1;                // (lmax + 1)(lmax + 2) / 2
      int nbase_max = 1, ncontr_max = 1;
      int nsrc_max = 0, nch_max = 0;
      double rcut = 0.0;
      bool analytic = false;      // radial functions always evaluated directly (debug)

      onika::memory::CudaMMVector<Ye3tBond> bonds;  // [nel][nel]
      onika::memory::CudaMMVector<double> pool;     // splines, contraction coefficients
      onika::memory::CudaMMVector<double> ylm_a, ylm_b;  // [nlh] Y_lm recurrence coefficients

      // per element
      onika::memory::CudaMMVector<double> e0, fs0, rho_lim;
      onika::memory::CudaMMVector<int> src0;        // [nel+1] sources of element mu: src0[mu] .. src0[mu+1]
      onika::memory::CudaMMVector<int> srcb;        // [nel][nel+1] local source range of neighbour element nu: srcb[mu*(nel+1)+nu] ..
      onika::memory::CudaMMVector<int> ch0;         // [nel+1]
      onika::memory::CudaMMVector<int> no0;         // [nel+1] product tree nodes of element mu

      // sources (local index q = global - src0[mu]), sorted by neighbour element
      onika::memory::CudaMMVector<int> src_kind;    // 0: radial base g_k, 1: R_nl Y_lm
      onika::memory::CudaMMVector<int> src_r;       // k, or n (lb + 1) + l
      onika::memory::CudaMMVector<int> src_lm;      // l (l + 1) / 2 + m (0 for the base)
      onika::memory::CudaMMVector<int> src_nu;      // neighbour element

      // channels: local source, transform (0: the source, +1 / -1: +- conj of it)
      onika::memory::CudaMMVector<int> ch_src, ch_tr;

      // monomials (one of each complex-conjugate pair, real part) as a prefix tree of their sorted factor lists (powers
      // expanded), pre-order: factor channel (local), depth (1 = first factor), coefficient of the monomial ending there
      onika::memory::CudaMMVector<int> node_ch, node_depth;
      onika::memory::CudaMMVector<double> node_coef;

      size_t bytes() const
      {
        return bonds.size() * sizeof(Ye3tBond) + pool.size() * 8 + node_coef.size() * 8
             + ( 2 * node_ch.size() + ch_src.size() + ch_tr.size() + 4 * src_kind.size() ) * 4;
      }
    };

    // raw view given to the kernels
    struct Ye3tView
    {
      int nel = 0, lmax = 0, nlh = 1, nbase_max = 1, ncontr_max = 1, nsrc_max = 0;
      bool analytic = false;
      const Ye3tBond* bonds = nullptr;
      const double* pool = nullptr;
      const double *ylm_a = nullptr, *ylm_b = nullptr;
      const double *e0 = nullptr, *fs0 = nullptr, *rho_lim = nullptr;
      const int *src0 = nullptr, *srcb = nullptr, *ch0 = nullptr, *no0 = nullptr;
      const int *src_kind = nullptr, *src_r = nullptr, *src_lm = nullptr, *src_nu = nullptr;
      const int *ch_src = nullptr, *ch_tr = nullptr;
      const int *node_ch = nullptr, *node_depth = nullptr;
      const double* node_coef = nullptr;
    };

    inline Ye3tView make_view(const Ye3tTables& T)
    {
      Ye3tView v;
      v.nel = T.nel; v.lmax = T.lmax; v.nlh = T.nlh; v.nbase_max = T.nbase_max; v.ncontr_max = T.ncontr_max; v.nsrc_max = T.nsrc_max;
      v.analytic = T.analytic;
      v.bonds = T.bonds.data(); v.pool = T.pool.data(); v.ylm_a = T.ylm_a.data(); v.ylm_b = T.ylm_b.data();
      v.e0 = T.e0.data(); v.fs0 = T.fs0.data(); v.rho_lim = T.rho_lim.data();
      v.src0 = T.src0.data(); v.srcb = T.srcb.data(); v.ch0 = T.ch0.data(); v.no0 = T.no0.data();
      v.src_kind = T.src_kind.data(); v.src_r = T.src_r.data(); v.src_lm = T.src_lm.data(); v.src_nu = T.src_nu.data();
      v.ch_src = T.ch_src.data(); v.ch_tr = T.ch_tr.data();
      v.node_ch = T.node_ch.data(); v.node_depth = T.node_depth.data(); v.node_coef = T.node_coef.data();
      return v;
    }

    // per-chunk buffers, kept from step to step
    struct Ye3tBuffers
    {
      onika::memory::CudaMMVector<int> order;       // central atoms ordered by element
      onika::memory::CudaMMVector<double> A, W;     // [2][nsrc_max][chunk]: sources, dE/d(source) (re, im)
      onika::memory::CudaMMVector<double> C, G;     // [2][nch_max][chunk]: channels, d(E)/d(channel)
      onika::memory::CudaMMVector<double> scratch;  // block scratch
      onika::memory::CudaMMVector<int> flags;       // [1] density guard violations: rho >= limit
    };

  }
}
