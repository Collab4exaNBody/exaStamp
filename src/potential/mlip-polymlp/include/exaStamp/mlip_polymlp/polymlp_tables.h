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

// Flat tables of a polynomial machine-learning potential (A. Seko, pypolymlp / lammps-polymlp-package), built by
// polymlp_init from the file written by scripts/polymlp_export, read by the kernels of polymlp_kernels.h.
//
// Model, for a central atom i of type t (pair type tp = tp(t, t_j) of its neighbours j):
//   radial functions f_n(r) = exp(-beta_n (r - mu_n)^2) fc(r), fc = (1 + cos(pi r / rc)) / 2
//   sources (m >= 0): A_{nu n} = sum_{j of type nu} f_n(r_ij) (pair features), A_{nu n l m} = sum_j f_n(r_ij) Y_lm(r_ij / r)
//     (gtinv), Y_lm complex Condon-Shortley normalised to Y_00 = 1, r_ij = r_j - r_i
//   order parameters (channels): a = s A or s conj(A), s a real factor (orthonormal harmonics 1 / sqrt(4 pi), polymlp's
//     bond vector r_i - r_j: (-1)^l, negative m: (-1)^m conj)
//   features d_k = Re sum_terms c prod_factors a  (pair: d_k = a_k); the products are shared along a prefix tree of the
//     sorted factor lists, walked depth first (forward values, then reverse-mode adjoints) with a stack of depth max factors
//   E_i = sum_u w_u prod_{factors} d_k   (up to max_p factors)

#include <onika/memory/allocator.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace exaStamp
{
  namespace mlip_polymlp
  {
    static constexpr int kPolymlpMaxLmax = 20;      // largest angular degree
    static constexpr int kPolymlpMaxFactors = 8;    // largest number of factors of a feature term or polynomial term
    static constexpr int kPolymlpMaxRadial = 64;    // largest number of gaussians

    struct PolymlpTables
    {
      std::vector<std::string> elements;
      int nel = 0;
      double rc = 0.0;
      int ng = 0;                 // gaussians
      int lmax = 0, nlh = 1;
      int nsrc_max = 0, nfeat_max = 0, nch_max = 0;
      bool gtinv = true;

      onika::memory::CudaMMVector<double> beta, mu;       // [ng]
      onika::memory::CudaMMVector<double> ylm_a, ylm_b;   // [nlh]

      // per centre type t
      onika::memory::CudaMMVector<int> src0;      // [nel+1]
      onika::memory::CudaMMVector<int> srcb;      // [nel][nel+1] local source range of neighbour type nu
      onika::memory::CudaMMVector<int> ch0;       // [nel+1]
      onika::memory::CudaMMVector<int> ft0;       // [nel+1] features
      onika::memory::CudaMMVector<int> po0;       // [nel+1] polynomial terms

      // sources (local index q = global - src0[t]), sorted by neighbour type
      onika::memory::CudaMMVector<int> src_n, src_lm, src_nu;
      // channels: local source, conjugate (0/1), factor
      onika::memory::CudaMMVector<int> ch_src, ch_conj;
      onika::memory::CudaMMVector<double> ch_scale;
      // feature terms as a prefix tree of their (sorted) channel factors, in pre-order: node channel (local), depth
      // (1 = a factor), terms ending at the node nto[n] .. nto[n+1]: local feature, coefficient
      onika::memory::CudaMMVector<int> no0;       // [nel+1] nodes of type t
      onika::memory::CudaMMVector<int> node_ch, node_depth;
      onika::memory::CudaMMVector<int64_t> nto;   // [nnodes + 1]
      onika::memory::CudaMMVector<int> nt_feat;
      onika::memory::CudaMMVector<double> nt_coef;
      // polynomial terms: coefficient, factors pfo[u] .. pfo[u+1] (local features)
      onika::memory::CudaMMVector<double> pcoef;
      onika::memory::CudaMMVector<int64_t> pfo;
      onika::memory::CudaMMVector<int> pfac;

      size_t bytes() const
      {
        return 8 * ( nt_coef.size() + pcoef.size() + nto.size() + pfo.size() + ch_scale.size() )
             + 4 * ( node_ch.size() * 2 + nt_feat.size() + pfac.size() + src_n.size() * 3 + ch_src.size() * 2 );
      }
    };

    struct PolymlpView
    {
      int nel = 0, ng = 0, lmax = 0, nlh = 1, nsrc_max = 0, nfeat_max = 0, nch_max = 0;
      double rc = 0.0;
      bool gtinv = true;
      const double *beta = nullptr, *mu = nullptr, *ylm_a = nullptr, *ylm_b = nullptr;
      const int *src0 = nullptr, *srcb = nullptr, *ch0 = nullptr, *ft0 = nullptr, *po0 = nullptr, *no0 = nullptr;
      const int *src_n = nullptr, *src_lm = nullptr, *src_nu = nullptr, *ch_src = nullptr, *ch_conj = nullptr;
      const double* ch_scale = nullptr;
      const int *node_ch = nullptr, *node_depth = nullptr, *nt_feat = nullptr;
      const int64_t *nto = nullptr, *pfo = nullptr;
      const double *nt_coef = nullptr, *pcoef = nullptr;
      const int* pfac = nullptr;
    };

    inline PolymlpView make_view(const PolymlpTables& T)
    {
      PolymlpView v;
      v.nel = T.nel; v.ng = T.ng; v.lmax = T.lmax; v.nlh = T.nlh; v.nsrc_max = T.nsrc_max; v.nfeat_max = T.nfeat_max; v.nch_max = T.nch_max; v.rc = T.rc; v.gtinv = T.gtinv;
      v.beta = T.beta.data(); v.mu = T.mu.data(); v.ylm_a = T.ylm_a.data(); v.ylm_b = T.ylm_b.data();
      v.src0 = T.src0.data(); v.srcb = T.srcb.data(); v.ch0 = T.ch0.data(); v.ft0 = T.ft0.data(); v.po0 = T.po0.data(); v.no0 = T.no0.data();
      v.src_n = T.src_n.data(); v.src_lm = T.src_lm.data(); v.src_nu = T.src_nu.data(); v.ch_src = T.ch_src.data(); v.ch_conj = T.ch_conj.data();
      v.ch_scale = T.ch_scale.data();
      v.node_ch = T.node_ch.data(); v.node_depth = T.node_depth.data(); v.nt_feat = T.nt_feat.data(); v.nto = T.nto.data();
      v.nt_coef = T.nt_coef.data(); v.pfo = T.pfo.data(); v.pcoef = T.pcoef.data(); v.pfac = T.pfac.data();
      return v;
    }

    struct PolymlpBuffers
    {
      onika::memory::CudaMMVector<int> order;
      onika::memory::CudaMMVector<double> A, W;     // [2][nsrc_max][chunk]
      onika::memory::CudaMMVector<double> D, GD;    // [nfeat_max][chunk] features, dE/d(feature)
      onika::memory::CudaMMVector<double> C, G;     // [2][nch_max][chunk] channels, dE/d(channel)
      onika::memory::CudaMMVector<double> scratch;
    };
  }
}
