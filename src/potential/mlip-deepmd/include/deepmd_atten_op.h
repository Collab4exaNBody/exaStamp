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

// DeePMD se_atten / DPA1 force functor (host), one central atom i per call. The math (Zhang et al., DPA-1, 2022):
//
//   neighbors : the nnei nearest atoms within rcut, any type (empty slots counted, p = nnei - n)
//   Rh_j      = ( sw(r) (1/r, d/r^2) - davg[ti] ) / dstd[ti],  s_j = Rh_j[0],  u_j = Rh_j[1:4] / |Rh_j[1:4]|
//   G_j       = concat: emb_{ti,tj}(s_j) (type embeddings folded into the first layer bias by the exporter)
//               strip : g(s_j) ( 1 + T_{ti,tj} sw_j ) (smooth) or g(s_j) ( 1 + T_{ti,tj} )
//   attention layers, X = G :
//     q, k, v  = X W_in + b_in (each normalized per row if 'normalize'),  L_ij = scaling q_i . k_j
//     smooth   : L'_ij = ( L_ij + c ) sw_i sw_j - c,  A'_ij = softmax_j over the nnei slots (empty slots: logit -c),
//                A_ij = A'_ij sw_i sw_j
//     otherwise: A_ij = softmax over the real neighbors
//     dotr     : A_ij *= u_i . u_j
//     X       <- LayerNorm( X + ( A v ) W_out + b_out )
//   gr        = ( sum_j X_j x Rh_j + p Xpad x Rh_pad ) / nnei,  D = gr gr[:axis]^T,  descriptor = ( D, tebd_i )
//   E_i       = fitting(descriptor) + bias[ti]
// Empty slots are type independent (the exporter refuses models where they are not): their row Rh_pad = -davg/dstd
// gets no attention (masked) and its final feature Xpad is precomputed by deepmd_init (ctx.pad).
//
// Backward: hand-written vector-Jacobian products through the fitting network, the descriptor, the attention layers
// (layer norm, output projection, attention weights, softmax, smooth switch factors, angular factors, row
// normalizations, input projection) and the embedding networks, to dE/dRh_j, dE/dsw_j and dE/du_j, then dE/dd_j.
// Per-atom scratch O(nnei^2 + nnei (ng + 3 hidden)) per layer: thread-local buffers, host only (CudaCompatible false).

#include <cmath>
#include <vector>
#include <algorithm>
#include <exanb/compute/compute_pair_buffer.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <exanb/compute/compute_pair_traits.h>
#include <exanb/core/concurent_add_contributions.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>

#include "deepmd_model.h"
#include "deepmd_force_op.h"

namespace exaStamp
{
  using namespace exanb;

  // per-thread scratch of one central atom, grown on demand
  struct DeepmdAttenWork
  {
    // neighbors
    std::vector<int> nb, tj;                 // kept slots of the compute buffer, their model types
    std::vector<double> rh, u, vn, sw, dsw;  // [n][4], [n][3], |Rh[1:4]|, switch, its radial derivative
    std::vector<DeepmdForceOp::EnvRow> env;
    // per layer : input, normalized q k v and norms, raw logits, softmax, layer norm
    std::vector<double> X;                   // [nlayer + 1][n][ng] (X[0] = G, X[L] = output of layer L - 1)
    std::vector<double> q, k, v, nq, nk, nv; // [nlayer][n][hd], [nlayer][n]
    std::vector<double> lg, ab;              // [nlayer][n][n]
    std::vector<double> xh, sig;             // [nlayer][n][ng], [nlayer][n]
    // backward
    std::vector<double> dX, dXin, dO, dP, dA, dS, dW, du, dsw_acc, drh;
    std::vector<double> gr, dgr, desc, ddesc;
    std::vector<double> gs;                  // strip: g(s_j) [n][ng]
    DPMlpWork<DP_MAX_FIT_LAYERS, DP_MAX_FIT_WIDTH> fit;
    DPEmbeddingWork emb;
    std::vector<double> fg0, fg1;
    int nd = 0;                              // descriptor size
  };

  // attention model view
  struct DeepmdAttenView
  {
    const DPAttnLayer* layers = nullptr;
    int nlayers = 0, hidden = 0, tebd_dim = 0;
    double scaling = 1.0, shift = 20.0;
    bool normalize = true, smooth = true, dotr = true, strip = false, concat_output_tebd = true;
    const double* tebd = nullptr;            // [ntypes][tebd_dim]
    const double* strip_table = nullptr;     // [ntypes or ntypes^2][ng]
  };

  struct DeepmdAttenOp
  {
    DeepmdForceOp base;   // model view, environment rows, fitting networks
    DeepmdAttenView a;

    template <class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }
    template <class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3d &virial, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }
    template <class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    // forward part: neighbors, attention layers, descriptor in W.desc (W.nd values); the backward part's state stays in W
    template <class ComputeBufferT>
    inline void forward(size_t nn, ComputeBufferT &buf, int type, DeepmdAttenWork& W) const
    {
      const DeepmdModelView& m = base.m;
      const int ti = m.species_type[type];
      const int nt = m.ntypes, ng = m.ng, m2 = m.axis_neuron, nnei = m.nnei, hd = a.hidden, NL = a.nlayers;
      const double inv_nnei = 1.0 / nnei;

      // ---- the nnei nearest neighbors within rcut, any type
      const double rc2 = m.rcut * m.rcut;
      W.nb.clear();
      for(size_t jj = 0; jj < nn; ++jj) if( buf.d2[jj] < rc2 && buf.d2[jj] > 0.0 ) W.nb.push_back( int(jj) );
      if( int(W.nb.size()) > nnei )
      {
        std::nth_element( W.nb.begin(), W.nb.begin() + nnei, W.nb.end(), [&](int x, int y){ return buf.d2[x] < buf.d2[y]; } );
        W.nb.resize( nnei );
      }
      const int n = W.nb.size();
      const int npad = nnei - n;

      auto grow = [](std::vector<double>& v, size_t s) { if( v.size() < s ) v.resize(s); };
      W.tj.resize(n); W.env.resize(n);
      grow(W.rh, size_t(n) * 4); grow(W.u, size_t(n) * 3); grow(W.vn, n); grow(W.sw, n); grow(W.dsw, n);
      grow(W.X, size_t(NL + 1) * n * ng);
      grow(W.q, size_t(NL) * n * hd); grow(W.k, size_t(NL) * n * hd); grow(W.v, size_t(NL) * n * hd);
      grow(W.nq, size_t(NL) * n); grow(W.nk, size_t(NL) * n); grow(W.nv, size_t(NL) * n);
      grow(W.lg, size_t(NL) * n * n); grow(W.ab, size_t(NL) * n * n);
      grow(W.xh, size_t(NL) * n * ng); grow(W.sig, size_t(NL) * n);
      grow(W.gs, size_t(n) * ng);

      // ---- environment rows, unit vectors, switches, embeddings
      for(int q = 0; q < n; q++)
      {
        const int jj = W.nb[q];
        const int tj = m.species_type[ buf.nbh_pt[jj][field::type] ];
        W.tj[q] = tj;
        const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
        DeepmdForceOp::EnvRow& e = W.env[q];
        base.env_row(d, buf.d2[jj], ti, tj, e);
        for(int c = 0; c < 4; c++) W.rh[q*4 + c] = e.rh[c];
        const double vn = std::sqrt( e.rh[1]*e.rh[1] + e.rh[2]*e.rh[2] + e.rh[3]*e.rh[3] );
        W.vn[q] = vn;
        const double inv = 1.0 / std::max( vn , 1e-12 );
        for(int c = 0; c < 3; c++) W.u[q*3 + c] = e.rh[1 + c] * inv;
        W.sw[q] = e.w;
        W.dsw[q] = e.dw;
        double* G = W.X.data() + size_t(q) * ng;
        if( !a.strip )
        {
          const double* g = dp_mlp_forward( base.embedding_net(ti, tj), m.layers, m.weights, e.rh, W.emb );
          for(int c = 0; c < ng; c++) G[c] = g[c];
        }
        else
        {
          const double* g = dp_mlp_forward( m.embedding[0], m.layers, m.weights, e.rh, W.emb );
          const double* T = a.strip_table + size_t( m.type_one_side ? tj : ti * nt + tj ) * ng;
          const double f = a.smooth ? e.w : 1.0;
          for(int c = 0; c < ng; c++) { W.gs[size_t(q)*ng + c] = g[c]; G[c] = g[c] * ( 1.0 + T[c] * f ); }
        }
      }

      // ---- attention layers (forward, intermediates kept)
      for(int l = 0; l < NL; l++)
      {
        const DPAttnLayer& L = a.layers[l];
        const double *Win = m.weights + L.in_w, *bin = m.weights + L.in_b, *Wout = m.weights + L.out_w, *bout = m.weights + L.out_b;
        const double *lnw = m.weights + L.ln_w, *lnb = m.weights + L.ln_b;
        const double* X = W.X.data() + size_t(l) * n * ng;
        double* Xo = W.X.data() + size_t(l + 1) * n * ng;
        double *Q = W.q.data() + size_t(l) * n * hd, *K = W.k.data() + size_t(l) * n * hd, *V = W.v.data() + size_t(l) * n * hd;
        double *NQ = W.nq.data() + size_t(l) * n, *NK = W.nk.data() + size_t(l) * n, *NV = W.nv.data() + size_t(l) * n;
        double *LG = W.lg.data() + size_t(l) * n * n, *AB = W.ab.data() + size_t(l) * n * n;
        double *XH = W.xh.data() + size_t(l) * n * ng, *SG = W.sig.data() + size_t(l) * n;
        for(int i = 0; i < n; i++)
        {
          double* qi = Q + size_t(i) * hd; double* ki = K + size_t(i) * hd; double* vi = V + size_t(i) * hd;
          for(int c = 0; c < hd; c++) { qi[c] = bin[c]; ki[c] = bin[hd + c]; vi[c] = bin[2*hd + c]; }
          const double* xi = X + size_t(i) * ng;
          for(int a2 = 0; a2 < ng; a2++)
          {
            const double x = xi[a2];
            const double* w = Win + size_t(a2) * 3 * hd;
            for(int c = 0; c < hd; c++) { qi[c] += x * w[c]; ki[c] += x * w[hd + c]; vi[c] += x * w[2*hd + c]; }
          }
          if( a.normalize )
          {
            double sq = 0, sk = 0, sv = 0;
            for(int c = 0; c < hd; c++) { sq += qi[c]*qi[c]; sk += ki[c]*ki[c]; sv += vi[c]*vi[c]; }
            NQ[i] = std::sqrt(sq); NK[i] = std::sqrt(sk); NV[i] = std::sqrt(sv);
            for(int c = 0; c < hd; c++) { qi[c] /= NQ[i]; ki[c] /= NK[i]; vi[c] /= NV[i]; }
          }
        }
        for(int i = 0; i < n; i++)
        {
          const double* qi = Q + size_t(i) * hd;
          double* li = LG + size_t(i) * n;
          double* ai = AB + size_t(i) * n;
          double mx = -1e300;
          for(int j = 0; j < n; j++)
          {
            const double* kj = K + size_t(j) * hd;
            double s = 0.0;
            for(int c = 0; c < hd; c++) s += qi[c] * kj[c];
            s *= a.scaling;
            li[j] = s;
            const double lp = a.smooth ? ( s + a.shift ) * W.sw[i] * W.sw[j] - a.shift : s;
            ai[j] = lp;
            mx = std::max( mx , lp );
          }
          if( a.smooth && npad > 0 ) mx = std::max( mx , -a.shift );
          double z = 0.0;
          for(int j = 0; j < n; j++) { ai[j] = std::exp( ai[j] - mx ); z += ai[j]; }
          if( a.smooth ) z += npad * std::exp( -a.shift - mx );
          for(int j = 0; j < n; j++) ai[j] /= z;    // softmax A'
        }
        // O = A v with A = A' [sw_i sw_j] [u_i . u_j], output projection, residual, layer norm
        std::vector<double>& O = W.dO; grow(O, size_t(n) * hd);
        for(int i = 0; i < n; i++)
        {
          double* oi = O.data() + size_t(i) * hd;
          for(int c = 0; c < hd; c++) oi[c] = 0.0;
          const double* ai = AB + size_t(i) * n;
          for(int j = 0; j < n; j++)
          {
            double w = ai[j];
            if( a.smooth ) w *= W.sw[i] * W.sw[j];
            if( a.dotr ) w *= W.u[i*3]*W.u[j*3] + W.u[i*3+1]*W.u[j*3+1] + W.u[i*3+2]*W.u[j*3+2];
            const double* vj = V + size_t(j) * hd;
            for(int c = 0; c < hd; c++) oi[c] += w * vj[c];
          }
          double* yo = Xo + size_t(i) * ng;
          const double* xi = X + size_t(i) * ng;
          for(int c = 0; c < ng; c++) yo[c] = xi[c] + bout[c];
          for(int h = 0; h < hd; h++)
          {
            const double o = oi[h];
            const double* w = Wout + size_t(h) * ng;
            for(int c = 0; c < ng; c++) yo[c] += o * w[c];
          }
          double mu = 0.0, var = 0.0;
          for(int c = 0; c < ng; c++) mu += yo[c];
          mu /= ng;
          for(int c = 0; c < ng; c++) var += ( yo[c] - mu ) * ( yo[c] - mu );
          var /= ng;
          const double sg = std::sqrt( var + L.ln_eps );
          SG[i] = sg;
          double* xh = XH + size_t(i) * ng;
          for(int c = 0; c < ng; c++) { xh[c] = ( yo[c] - mu ) / sg; yo[c] = xh[c] * lnw[c] + lnb[c]; }
        }
      }
      const double* Xf = W.X.data() + size_t(NL) * n * ng;

      // ---- descriptor, fitting
      grow(W.gr, size_t(ng) * 4);
      for(int k = 0; k < ng * 4; k++) W.gr[k] = npad * m.pad[ size_t( ti * nt ) * ng * 4 + k ];
      for(int q = 0; q < n; q++)
      {
        const double* x = Xf + size_t(q) * ng;
        const double* r = W.rh.data() + q * 4;
        for(int c = 0; c < ng; c++) for(int k = 0; k < 4; k++) W.gr[c*4 + k] += x[c] * r[k];
      }
      for(int k = 0; k < ng * 4; k++) W.gr[k] *= inv_nnei;
      const int nd = ng * m2 + ( a.concat_output_tebd ? a.tebd_dim : 0 );
      grow(W.desc, nd); grow(W.ddesc, nd);
      for(int k = 0; k < ng * m2; k++)
      {
        const double* ga = W.gr.data() + ( k / m2 ) * 4;
        const double* gb = W.gr.data() + ( k % m2 ) * 4;
        W.desc[k] = ga[0] * gb[0] + ga[1] * gb[1] + ga[2] * gb[2] + ga[3] * gb[3];
      }
      if( a.concat_output_tebd ) for(int c = 0; c < a.tebd_dim; c++) W.desc[ng * m2 + c] = a.tebd[ ti * a.tebd_dim + c ];
      W.nd = nd;
    }

    // backward part, from dE/ddescriptor in W.ddesc: forces on the neighbors (atomic), the central atom force fi and
    // virial vir accumulated
    template <bool compute_virial, class ComputeBufferT, class CellParticlesT>
    inline void backward(ComputeBufferT &buf, int type, CellParticlesT cells, DeepmdAttenWork& W, Vec3d& fi, Mat3d& vir) const
    {
      const DeepmdModelView& m = base.m;
      const int ti = m.species_type[type];
      const int nt = m.ntypes, ng = m.ng, m2 = m.axis_neuron, nnei = m.nnei, hd = a.hidden, NL = a.nlayers;
      const double inv_nnei = 1.0 / nnei;
      const int n = W.nb.size();
      const double* Xf = W.X.data() + size_t(NL) * n * ng;
      auto grow = [](std::vector<double>& v, size_t s) { if( v.size() < s ) v.resize(s); };
      (void) nnei;
      // dE/dgr (unnormalized sum), dE/dX final and dE/dRh from gr
      grow(W.dgr, size_t(ng) * 4);
      for(int k = 0; k < ng * 4; k++)
      {
        const int c = k / 4, kk = k % 4;
        double s = 0.0;
        for(int b = 0; b < m2; b++) s += W.ddesc[c*m2 + b] * W.gr[b*4 + kk];
        if( c < m2 ) for(int c2 = 0; c2 < ng; c2++) s += W.ddesc[c2*m2 + c] * W.gr[c2*4 + kk];
        W.dgr[k] = s * inv_nnei;
      }
      grow(W.dX, size_t(n) * ng); grow(W.drh, size_t(n) * 4); grow(W.du, size_t(n) * 3); grow(W.dsw_acc, n);
      for(int q = 0; q < n; q++)
      {
        const double* r = W.rh.data() + q * 4;
        const double* x = Xf + size_t(q) * ng;
        double* dx = W.dX.data() + size_t(q) * ng;
        double* dr = W.drh.data() + q * 4;
        for(int k = 0; k < 4; k++) dr[k] = 0.0;
        for(int c = 0; c < ng; c++)
        {
          const double* g = W.dgr.data() + c * 4;
          dx[c] = g[0]*r[0] + g[1]*r[1] + g[2]*r[2] + g[3]*r[3];
          for(int k = 0; k < 4; k++) dr[k] += g[k] * x[c];
        }
        W.du[q*3] = W.du[q*3+1] = W.du[q*3+2] = 0.0;
        W.dsw_acc[q] = 0.0;
      }

      // ---- attention layers (backward)
      grow(W.dXin, size_t(n) * ng); grow(W.dO, size_t(n) * hd); grow(W.dP, size_t(n) * 3 * hd);
      grow(W.dA, size_t(n) * n); grow(W.dS, size_t(n) * n); grow(W.dW, size_t(n) * n);
      for(int l = NL - 1; l >= 0; l--)
      {
        const DPAttnLayer& L = a.layers[l];
        const double *Win = m.weights + L.in_w, *Wout = m.weights + L.out_w, *lnw = m.weights + L.ln_w;
        const double *Q = W.q.data() + size_t(l) * n * hd, *K = W.k.data() + size_t(l) * n * hd, *V = W.v.data() + size_t(l) * n * hd;
        const double *NQ = W.nq.data() + size_t(l) * n, *NK = W.nk.data() + size_t(l) * n, *NV = W.nv.data() + size_t(l) * n;
        const double *LG = W.lg.data() + size_t(l) * n * n, *AB = W.ab.data() + size_t(l) * n * n;
        const double *XH = W.xh.data() + size_t(l) * n * ng, *SG = W.sig.data() + size_t(l) * n;
        double* dX = W.dX.data();      // dE/d(output of layer l) on entry, dE/d(input of layer l) on exit
        double* dY = W.dXin.data();
        // layer norm backward: dY = ( dxh - mean(dxh) - xh mean(dxh xh) ) / sig, dxh = dX lnw
        for(int i = 0; i < n; i++)
        {
          const double* xh = XH + size_t(i) * ng;
          const double* dx = dX + size_t(i) * ng;
          double* dy = dY + size_t(i) * ng;
          double m1 = 0.0, m2s = 0.0;
          for(int c = 0; c < ng; c++) { const double g = dx[c] * lnw[c]; dy[c] = g; m1 += g; m2s += g * xh[c]; }
          m1 /= ng; m2s /= ng;
          for(int c = 0; c < ng; c++) dy[c] = ( dy[c] - m1 - xh[c] * m2s ) / SG[i];
        }
        // residual: dX_in = dY (+ the attention path below); dO = dY W_out^T
        for(int i = 0; i < n; i++)
        {
          const double* dy = dY + size_t(i) * ng;
          double* dO = W.dO.data() + size_t(i) * hd;
          for(int h = 0; h < hd; h++)
          {
            const double* w = Wout + size_t(h) * ng;
            double s = 0.0;
            for(int c = 0; c < ng; c++) s += dy[c] * w[c];
            dO[h] = s;
          }
        }
        // A_ij = A'_ij S_ij Om_ij : dA_ij = dO_i . v_j, dv_j = sum_i A_ij dO_i
        double* dPq = W.dP.data();      // [n][3 hd] : dq, dk, dv
        for(int k = 0; k < n * 3 * hd; k++) dPq[k] = 0.0;
        for(int i = 0; i < n; i++)
        {
          const double* dO = W.dO.data() + size_t(i) * hd;
          const double* ai = AB + size_t(i) * n;
          for(int j = 0; j < n; j++)
          {
            const double* vj = V + size_t(j) * hd;
            double s = 0.0;
            for(int h = 0; h < hd; h++) s += dO[h] * vj[h];
            const double S = a.smooth ? W.sw[i] * W.sw[j] : 1.0;
            const double Om = a.dotr ? W.u[i*3]*W.u[j*3] + W.u[i*3+1]*W.u[j*3+1] + W.u[i*3+2]*W.u[j*3+2] : 1.0;
            const double A = ai[j] * S * Om;
            double* dvj = dPq + size_t(j) * 3 * hd + 2 * hd;
            for(int h = 0; h < hd; h++) dvj[h] += A * dO[h];
            W.dA[size_t(i)*n + j] = s * S * Om;                        // dE/dA'
            W.dS[size_t(i)*n + j] = a.smooth ? s * ai[j] * Om : 0.0;   // dE/dS (post-softmax factor)
            W.dW[size_t(i)*n + j] = a.dotr ? s * ai[j] * S : 0.0;      // dE/dOm
          }
        }
        // softmax backward: dL'_ij = A'_ij ( dA'_ij - sum_k dA'_ik A'_ik ) (empty slots: dA' = 0); smooth logits
        for(int i = 0; i < n; i++)
        {
          const double* ai = AB + size_t(i) * n;
          double* da = W.dA.data() + size_t(i) * n;
          double t = 0.0;
          for(int j = 0; j < n; j++) t += da[j] * ai[j];
          for(int j = 0; j < n; j++)
          {
            double dl = ai[j] * ( da[j] - t );   // dE/dL'
            if( a.smooth )
            {
              const double S = W.sw[i] * W.sw[j];
              W.dS[size_t(i)*n + j] += dl * ( LG[size_t(i)*n + j] + a.shift );
              dl *= S;
            }
            da[j] = dl * a.scaling;              // dE/d(q_i . k_j)
          }
        }
        // q_i . k_j : dq_i = sum_j dL_ij k_j, dk_j = sum_i dL_ij q_i
        for(int i = 0; i < n; i++)
        {
          const double* dl = W.dA.data() + size_t(i) * n;
          const double* qi = Q + size_t(i) * hd;
          double* dqi = dPq + size_t(i) * 3 * hd;
          for(int j = 0; j < n; j++)
          {
            const double g = dl[j];
            const double* kj = K + size_t(j) * hd;
            double* dkj = dPq + size_t(j) * 3 * hd + hd;
            for(int h = 0; h < hd; h++) { dqi[h] += g * kj[h]; dkj[h] += g * qi[h]; }
          }
        }
        // S_ij = sw_i sw_j, Om_ij = u_i . u_j
        for(int i = 0; i < n; i++)
          for(int j = 0; j < n; j++)
          {
            const double ds = W.dS[size_t(i)*n + j];
            if( ds != 0.0 ) { W.dsw_acc[i] += ds * W.sw[j]; W.dsw_acc[j] += ds * W.sw[i]; }
            const double dw = W.dW[size_t(i)*n + j];
            if( dw != 0.0 ) for(int c = 0; c < 3; c++) { W.du[i*3 + c] += dw * W.u[j*3 + c]; W.du[j*3 + c] += dw * W.u[i*3 + c]; }
          }
        // row normalizations: dx = ( dxh - xh (xh . dxh) ) / |x|
        if( a.normalize )
          for(int i = 0; i < n; i++)
          {
            const double* xs[3] = { Q + size_t(i) * hd, K + size_t(i) * hd, V + size_t(i) * hd };
            const double nr[3] = { NQ[i], NK[i], NV[i] };
            for(int t = 0; t < 3; t++)
            {
              double* g = dPq + size_t(i) * 3 * hd + t * hd;
              double p = 0.0;
              for(int h = 0; h < hd; h++) p += xs[t][h] * g[h];
              for(int h = 0; h < hd; h++) g[h] = ( g[h] - xs[t][h] * p ) / nr[t];
            }
          }
        // input projection: dX_in = dY + dP W_in^T
        for(int i = 0; i < n; i++)
        {
          const double* g = dPq + size_t(i) * 3 * hd;
          const double* dy = dY + size_t(i) * ng;
          double* dx = dX + size_t(i) * ng;
          for(int c = 0; c < ng; c++)
          {
            const double* w = Win + size_t(c) * 3 * hd;
            double s = dy[c];
            for(int h = 0; h < 3 * hd; h++) s += g[h] * w[h];
            dx[c] = s;
          }
        }
      }

      // ---- embeddings backward, forces
      const double conv = m.conv_energy;
      double g0[DP_MAX_EMB_WIDTH], g1[DP_MAX_EMB_WIDTH];
      for(int q = 0; q < n; q++)
      {
        const int jj = W.nb[q], tj = W.tj[q];
        const DeepmdForceOp::EnvRow& e = W.env[q];
        const double* dG = W.dX.data() + size_t(q) * ng;
        double* drh = W.drh.data() + q * 4;
        double ds = 0.0;
        if( !a.strip )
        {
          const DPNet& emb = base.embedding_net(ti, tj);
          dp_mlp_forward( emb, m.layers, m.weights, e.rh, W.emb );
          for(int c = 0; c < ng; c++) g0[c] = dG[c];
          dp_mlp_backward( emb, m.layers, m.weights, e.rh, W.emb, g0, &ds, g0, g1 );
        }
        else
        {
          const DPNet& emb = m.embedding[0];
          dp_mlp_forward( emb, m.layers, m.weights, e.rh, W.emb );
          const double* T = a.strip_table + size_t( m.type_one_side ? tj : ti * nt + tj ) * ng;
          const double* gs = W.gs.data() + size_t(q) * ng;
          const double f = a.smooth ? e.w : 1.0;
          for(int c = 0; c < ng; c++)
          {
            g0[c] = dG[c] * ( 1.0 + T[c] * f );
            if( a.smooth ) W.dsw_acc[q] += dG[c] * gs[c] * T[c];
          }
          dp_mlp_backward( emb, m.layers, m.weights, e.rh, W.emb, g0, &ds, g0, g1 );
        }
        drh[0] += ds;
        // u = Rh[1:4] / max( |Rh[1:4]| , 1e-12 ) (below the floor, near rcut: u = Rh[1:4] / 1e-12)
        if( a.dotr )
        {
          const double* u = W.u.data() + q * 3;
          const double* du = W.du.data() + q * 3;
          if( W.vn[q] > 1e-12 )
          {
            const double p = u[0]*du[0] + u[1]*du[1] + u[2]*du[2];
            for(int c = 0; c < 3; c++) drh[1 + c] += ( du[c] - u[c] * p ) / W.vn[q];
          }
          else for(int c = 0; c < 3; c++) drh[1 + c] += du[c] / 1e-12;
        }
        const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
        Vec3d g = base.env_row_grad(d, e, drh, ti, tj);
        g = g + d * ( W.dsw_acc[q] * e.dw / e.r );
        g = g * conv;
        const Vec3d fj = g * -1.0;
        fi = fi + g;
        if constexpr ( compute_virial ) { vir = vir + tensor(fj, d); }
        size_t cell_j = 0, p_j = 0;
        buf.nbh.get(jj, cell_j, p_j);
        atomic_add_contribution(cells[cell_j][field::fx][p_j], fj.x);
        atomic_add_contribution(cells[cell_j][field::fy][p_j], fj.y);
        atomic_add_contribution(cells[cell_j][field::fz][p_j], fj.z);
      }
    }

    template <class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t nn, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3dT &virial, CellParticlesT cells, GridCellLocksT, ParticleLockT &) const
    {
      static constexpr bool compute_virial = std::is_same_v<Mat3dT, Mat3d>;
      static thread_local DeepmdAttenWork W;
      const DeepmdModelView& m = base.m;
      const int ti = m.species_type[type];
      forward(nn, buf, type, W);
      const int nd = W.nd;
      const DPNet& fit = m.fitting[ m.fitting_mixed_types ? 0 : ti ];
      const double e_atom = dp_mlp_forward( fit, m.layers, m.weights, W.desc.data(), W.fit )[0] + m.energy_bias[ti];
      if( W.fg0.size() < size_t( std::max(nd, DP_MAX_FIT_WIDTH) ) ) { W.fg0.resize( std::max(nd, DP_MAX_FIT_WIDTH) ); W.fg1.resize( std::max(nd, DP_MAX_FIT_WIDTH) ); }
      const double one = 1.0;
      dp_mlp_backward( fit, m.layers, m.weights, W.desc.data(), W.fit, &one, W.ddesc.data(), W.fg0.data(), W.fg1.data() );
      Vec3d fi = { 0.0, 0.0, 0.0 };
      Mat3d vir = { 0.,0.,0., 0.,0.,0., 0.,0.,0. };
      backward<compute_virial>(buf, type, cells, W, fi, vir);
      atomic_add_contribution(fx, fi.x);
      atomic_add_contribution(fy, fi.y);
      atomic_add_contribution(fz, fi.z);
      en += e_atom * m.conv_energy;
      if constexpr ( compute_virial ) { virial += vir; }
    }
  };

}

namespace exanb
{
  template <> struct ComputePairTraits<exaStamp::DeepmdAttenOp>
  {
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = false;   // host only: O(nnei^2) per-atom scratch
    static inline constexpr bool HasParticleContextStart = false;
    static inline constexpr bool HasParticleContextStop = false;
    static inline constexpr bool HasParticleContext = false;
    static inline constexpr bool RequiresNbhOptionalData = true;
    static inline constexpr bool ComputeEmptyNeighborhood = true;
  };
}
