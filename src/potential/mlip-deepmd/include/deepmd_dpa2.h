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

// DeePMD DPA2 (Zhang et al. 2024) evaluated natively on the host:
//
//   D0_i = [ repinit_i (se_atten block without attention), three-body_i (se_t_tebd block) ]   (pass A, neighbor buffers)
//   g1   = act( T D0 (+ Tt tebd_i) ),  g2_ij = act( Wg2 s_ij ),  h2_ij (normalized environment rows), repformer neighbors
//   layer: g1' = g1 + sum_m r_m u_m over u = [ act(Wself g1), conv, act(W1 [g1?, conv?, grrg(g2), drrd(g1_j)]) ] (deepmd's
//          order; conv = sum_j (Wp g2_ij) g1_j sw / nnei, or into the linear1 input as sum_j g2_ij (Wp g1_j) sw / nnei)
//          g2' = g2 + r0 act(W2 g2) + r1 LN( Wh attn(g2, h2) )  (not in the last layer), attn: multi-head, smooth, gated by
//          h2_a . h2_b, weights times h2_a . h2_b / sqrt(3)
//   descriptor = g1 (+ tebd), fitting network
//
// Phases (dpa2_compute): repformer node features g1 of every layer for all atoms (ghosts updated in between), the
// per-atom g2 of every layer stored; reverse phases fold dE/dg1 of the neighbors back to their owners, down to dE/dg1
// of layer 0 and dE/dD0 (pass C: the repinit / three-body backward from dE/dD0, neighbor buffers).

#include "deepmd_dpa3.h"
#include "deepmd_dpa3_driver.h"
#include "deepmd_force_op.h"
#include "deepmd_atten_op.h"

namespace exaStamp
{
  using namespace exanb;

  struct DPA2Layer
  {
    DPA3Lin linear1, linear2, self_mlp, proj, mapqk, mapv, head_map;
    size_t lnw = 0, lnb = 0;
    double lneps = 1e-5;
    bool chnnl2 = false, attn = false;
    std::vector<size_t> g1_res, g2_res;
  };

  struct DPA2Model
  {
    DPA3Model base;                      // weights W and the activation of the repformers
    std::vector<DPA2Layer> layers;
    DPA3Lin g2_embd, shape, tebd_tr;     // shape: none when identity
    std::vector<double> tebd, davg, dstd;
    int nt = 0, tdim = 0, g1d = 0, g2d = 0, axis = 0, nnei = 0, nh = 0, hid = 0, dim0 = 0, dim_rep = 0, dim_3b = 0;
    double rc = 0, rs = 0, prot = 0, s_sym = 1.0, shift = 20.0;
    bool out_conv = true, out_mlp = true, conv = true, grrg = true, drrd = true, gate = true, concat_tebd = false;
    std::shared_ptr<DeepmdContext> rep, three;
  };

  // stage profile (DPA2_PROF=1): wall time per evaluation; 'inner' stages are measured on OpenMP thread 0
  struct DPA2Prof
  {
    bool on = std::getenv("DPA2_PROF") != nullptr;
    int evals = 0;
    double t[24] = {};
    ~DPA2Prof()
    {
      if( !on || evals == 0 ) return;
      static const char* names[24] = { "pass A (repinit, 3-body)", "CSR, g1_0, g2 store", "forward layers", "fitting", "reverse layers", "pass C (backward)", "forces",
                                       "", "inner geometry+inputs", "inner fwd g1", "inner fwd linear2", "inner fwd attn maps", "inner fwd attn weights", "inner fwd head+LN",
                                       "inner vjp g1", "inner vjp linear2+LN", "inner vjp attn weights", "inner vjp attn maps", "inner vjp bonds+scatter", "inner repinit fwd", "inner 3-body fwd", "inner repinit bwd", "inner 3-body bwd", "" };
      std::fprintf( stderr, "[DPA2] profile, ms per evaluation (%d evaluations):\n", evals );
      for(int k = 0; k < 24; k++) if( t[k] > 0 ) std::fprintf( stderr, "  %-26s %10.3f\n", names[k], 1e3 * t[k] / evals );
    }
  };
  inline DPA2Prof& dpa2_prof() { static DPA2Prof p; return p; }
  struct DPA2Timer
  {
    bool on; double t0;
    inline DPA2Timer() : on( dpa2_prof().on && omp_get_thread_num() == 0 ), t0( on ? omp_get_wtime() : 0.0 ) {}
    inline void mark(int k) { if( on ) { const double t = omp_get_wtime(); dpa2_prof().t[k] += t - t0; t0 = t; } }
  };

  // ------------------------------------------------------------------ one atom of the repformers

  struct DPA2Geom
  {
    int n = 0, ti = 0;
    std::vector<int> j;
    std::vector<int64_t> bond;
    std::vector<double> d, r, q, w, dw, s, h2;
  };

  struct DPA2Work
  {
    std::vector<double> g1, g2;                      // inputs
    std::vector<double> a_self, d_self, P, c, hg2, hg1, mlp, a1, d1, a2, d2;
    std::vector<double> qk, v, L0, Ab, out, zh, xh, sig;
  };

  struct DPA2Scratch
  {
    DPA2Geom g;
    std::vector<DPA2Work> L;
    std::vector<int> order;
    std::vector<double> dg1, dg2, dg1_in, dg2_in, dgg1, gh2, gw, gs, tmp, tmp2, tmp3;
    std::vector<double> hh, ww, hw, aq, ak, av, ao, am, ad, dqh, dkh, dvh, dqk, dS, dHH, mx;   // attention, symmetrization
    std::vector<double> du, dz1, dmlp, dc, dM, gx;
  };

  inline void dpa2_geometry(const DPA2Model& M, int i, int ti, const int64_t* offset, const int* nbh, const double* rij, DPA2Scratch& S)
  {
    const DPA3Model& m = M.base;
    DPA2Geom& g = S.g;
    g.ti = ti;
    const int64_t b0 = offset[i], b1 = offset[i+1];
    const double rc2 = M.rc * M.rc;
    S.order.clear();
    for(int64_t b = b0; b < b1; b++)
    {
      const double* x = rij + 3 * b;
      const double d2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
      if( d2 < rc2 && d2 > 0.0 ) S.order.push_back( int( b - b0 ) );
    }
    std::stable_sort( S.order.begin(), S.order.end(), [&](int x, int y)
    {
      const double* a = rij + 3 * ( b0 + x ); const double* c = rij + 3 * ( b0 + y );
      return a[0]*a[0] + a[1]*a[1] + a[2]*a[2] < c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
    } );
    if( int(S.order.size()) > M.nnei ) S.order.resize( M.nnei );
    const int n = S.order.size();
    g.n = n;
    g.j.resize(n); g.bond.resize(n);
    dpa3_grow(g.d, size_t(n)*3); dpa3_grow(g.r, n); dpa3_grow(g.q, n); dpa3_grow(g.w, n); dpa3_grow(g.dw, n); dpa3_grow(g.s, n); dpa3_grow(g.h2, size_t(n)*3);
    const double* avg = M.davg.data() + 4 * ti;
    const double* sd = M.dstd.data() + 4 * ti;
    for(int k = 0; k < n; k++)
    {
      const int64_t b = b0 + S.order[k];
      g.bond[k] = b; g.j[k] = nbh[b];
      const double* x = rij + 3 * b;
      for(int c = 0; c < 3; c++) g.d[3*k + c] = x[c];
      const double r = std::sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
      g.r[k] = r; g.q[k] = r + M.prot;
      dpa3_switch( m, r, M.rs, M.rc, g.w[k], g.dw[k] );
      const double qq = g.q[k];
      for(int c = 0; c < 3; c++) g.h2[3*k + c] = ( g.w[k] * x[c] / ( qq * qq ) - avg[1 + c] ) / sd[1 + c];
      g.s[k] = ( g.w[k] / qq - avg[0] ) / sd[0];
    }
  }

  // y[n] = exp( y[n] )
  inline void dpa2_exp(double* y, size_t n)
  {
#   if defined(XSTAMP_DEEPMD_BLAS) && defined(__APPLE__)
    { const int nn = int(n); vvexp( y, y, &nn ); }
#   else
    for(size_t k = 0; k < n; k++) y[k] = std::exp( y[k] );
#   endif
  }
  // y[nout] (+)= x[nin] W, b added when given
  inline void dpa2_lin(const DPA3Model& m, const DPA3Lin& L, const double* x, double* y)
  {
    for(int o = 0; o < L.nout; o++) y[o] = L.b != size_t(-1) ? m.W[L.b + o] : 0.0;
    dpa3_gemm( false, false, 1, L.nout, L.nin, 1.0, x, L.nin, m.W.data() + L.w, L.nout, 1.0, y, L.nout );
  }
  // dx[nin] += W dy[nout]
  inline void dpa2_lin_t(const DPA3Model& m, const DPA3Lin& L, const double* dy, double* dx)
  {
    dpa3_gemm( false, false, L.nin, 1, L.nout, 1.0, m.W.data() + L.w, L.nout, dy, 1, 1.0, dx, 1 );
  }
  // sums over the neighbors: hw[n][3] = h2 sw s_sym
  inline void dpa2_hw(const DPA2Model& M, DPA2Scratch& S)
  {
    const DPA2Geom& g = S.g;
    dpa3_grow(S.hw, size_t(g.n) * 3);
    for(int k = 0; k < g.n; k++) { const double f = g.w[k] * M.s_sym; for(int x = 0; x < 3; x++) S.hw[3*k + x] = g.h2[3*k + x] * f; }
  }

  // initial g2 of the atom (into S.L[0].g2)
  inline void dpa2_init_g2(const DPA2Model& M, DPA2Scratch& S)
  {
    const DPA3Model& m = M.base;
    const int n = S.g.n, g2d = M.g2d;
    DPA2Work& w = S.L[0];
    dpa3_grow(w.g2, size_t(n) * g2d);
    for(int k = 0; k < n; k++)
    {
      double* y = w.g2.data() + size_t(k) * g2d;
      dpa3_lin( m, M.g2_embd, &S.g.s[k], y );
      for(int c = 0; c < g2d; c++) y[c] = dpa3_act( m, y[c] );
    }
  }

  // forward of repformer layer l: S.L[l] (g1, g2 inputs) -> S.L[l+1] (g1, g2), intermediates in S.L[l]
  inline void dpa2_layer_forward(const DPA2Model& M, int l, const double* Nl, DPA2Scratch& S)
  {
    const DPA3Model& m = M.base;
    const DPA2Layer& P = M.layers[l];
    const DPA2Geom& g = S.g;
    DPA2Work& w = S.L[l];
    DPA2Work& o = S.L[l+1];
    const int n = g.n, g1d = M.g1d, g2d = M.g2d, ax = M.axis, nh = M.nh, hid = M.hid;
    const double* W = m.W.data();
    const double* g1 = w.g1.data();
    const double* g2 = w.g2.data();
    DPA2Timer T;
    std::vector<const double*> upd;   // updates of g1 in deepmd's order
    dpa2_hw( M, S );
    // g1 of the neighbors
    std::vector<double>& Xg = S.tmp2; dpa3_grow(Xg, size_t(n) * g1d);
    for(int k = 0; k < n; k++) { const double* nj = Nl + size_t( g.j[k] ) * g1d; std::copy( nj, nj + g1d, Xg.data() + size_t(k) * g1d ); }
    const double* X1 = Xg.data();
    // self mlp
    if( M.out_mlp )
    {
      dpa3_grow(w.a_self, g1d); dpa3_grow(w.d_self, g1d);
      std::vector<double>& z = S.tmp; dpa3_grow(z, g1d);
      dpa2_lin( m, P.self_mlp, g1, z.data() );
      dpa3_actd( m, z.data(), g1d, w.a_self.data(), w.d_self.data() );
      upd.push_back( w.a_self.data() );
    }
    // linear1 input: [ g1 (if not out_mlp), conv (if not out_conv), grrg, drrd ]
    const int dmlp = ( M.out_mlp ? 0 : g1d ) + ( M.conv && !M.out_conv ? g2d : 0 ) + ( M.grrg ? ax * g2d : 0 ) + ( M.drrd ? ax * g1d : 0 );
    dpa3_grow(w.mlp, dmlp);
    int off = 0;
    if( !M.out_mlp ) { for(int c = 0; c < g1d; c++) w.mlp[c] = g1[c]; off += g1d; }
    if( M.conv )
    {
      const int pw = M.out_conv ? g1d : g2d;
      dpa3_grow(w.P, size_t(n) * pw);
      // out_conv: P = g2 Wp (g2d -> g1d); else P = g1_j Wp (g1d -> g2d)
      if( M.out_conv ) dpa3_gemm( false, false, n, g1d, g2d, 1.0, g2, g2d, W + P.proj.w, g1d, 0.0, w.P.data(), g1d );
      else
      {
        dpa3_gemm( false, false, n, g2d, g1d, 1.0, X1, g1d, W + P.proj.w, g2d, 0.0, w.P.data(), g2d );
      }
      dpa3_grow(w.c, pw);
      std::fill( w.c.begin(), w.c.begin() + pw, 0.0 );
      for(int k = 0; k < n; k++)
      {
        const double f = g.w[k] / M.nnei;
        const double* pk = w.P.data() + size_t(k) * pw;
        const double* xk = M.out_conv ? X1 + size_t(k) * g1d : g2 + size_t(k) * g2d;
        for(int c = 0; c < pw; c++) w.c[c] += pk[c] * xk[c] * f;
      }
      if( M.out_conv ) upd.push_back( w.c.data() );
      else { for(int c = 0; c < g2d; c++) w.mlp[off + c] = w.c[c]; off += g2d; }
    }
    if( M.grrg )
    {
      dpa3_grow(w.hg2, size_t(3) * g2d);
      dpa3_gemm( true, false, 3, g2d, n, 1.0, S.hw.data(), 3, g2, g2d, 0.0, w.hg2.data(), g2d );
      if( n == 0 ) std::fill( w.hg2.begin(), w.hg2.begin() + 3 * g2d, 0.0 );
      dpa3_grrg( w.hg2.data(), g2d, ax, w.mlp.data() + off ); off += ax * g2d;
    }
    if( M.drrd )
    {
      dpa3_grow(w.hg1, size_t(3) * g1d);
      dpa3_gemm( true, false, 3, g1d, n, 1.0, S.hw.data(), 3, X1, g1d, 0.0, w.hg1.data(), g1d );
      if( n == 0 ) std::fill( w.hg1.begin(), w.hg1.begin() + 3 * g1d, 0.0 );
      dpa3_grrg( w.hg1.data(), g1d, ax, w.mlp.data() + off ); off += ax * g1d;
    }
    dpa3_grow(w.a1, g1d); dpa3_grow(w.d1, g1d);
    {
      std::vector<double>& z = S.tmp; dpa3_grow(z, g1d);
      dpa2_lin( m, P.linear1, w.mlp.data(), z.data() );
      dpa3_actd( m, z.data(), g1d, w.a1.data(), w.d1.data() );
    }
    upd.push_back( w.a1.data() );
    o.g1.assign( g1, g1 + g1d );
    for(size_t u = 0; u < upd.size() && u < P.g1_res.size(); u++)
    {
      const double* r = W + P.g1_res[u];
      for(int c = 0; c < g1d; c++) o.g1[c] += r[c] * upd[u][c];
    }
    // g2 update
    o.g2.assign( g2, g2 + size_t(n) * g2d );
    T.mark(9);
    if( !P.chnnl2 ) return;
    dpa3_grow(w.a2, size_t(n) * g2d); dpa3_grow(w.d2, size_t(n) * g2d);
    {
      std::vector<double>& z = S.tmp; dpa3_grow(z, size_t(n) * g2d);
      for(int k = 0; k < n; k++) for(int c = 0; c < g2d; c++) z[size_t(k)*g2d + c] = W[P.linear2.b + c];
      dpa3_gemm( false, false, n, g2d, g2d, 1.0, g2, g2d, W + P.linear2.w, g2d, 1.0, z.data(), g2d );
      dpa3_actd( m, z.data(), size_t(n) * g2d, w.a2.data(), w.d2.data() );
    }
    const double* r0 = W + P.g2_res[0];
    for(size_t x = 0; x < size_t(n) * g2d; x++) o.g2[x] += r0[x % g2d] * w.a2[x];
    T.mark(10);
    if( !P.attn ) return;
    // attention: q, k (hid x 2 nh per edge), v (g2d x nh), weights per head
    const int nq = hid * 2 * nh, nv = g2d * nh;
    dpa3_grow(w.qk, size_t(n) * nq); dpa3_grow(w.v, size_t(n) * nv);
    dpa3_gemm( false, false, n, nq, g2d, 1.0, g2, g2d, W + P.mapqk.w, nq, 0.0, w.qk.data(), nq );
    dpa3_gemm( false, false, n, nv, g2d, 1.0, g2, g2d, W + P.mapv.w, nv, 0.0, w.v.data(), nv );
    T.mark(11);
    const size_t nn2 = size_t(n) * n;
    dpa3_grow(w.L0, size_t(nh) * nn2); dpa3_grow(w.Ab, size_t(nh) * nn2); dpa3_grow(w.out, size_t(n) * nv);
    std::fill( w.out.begin(), w.out.begin() + size_t(n) * nv, 0.0 );
    const int npad = M.nnei - n;
    const double ish = 1.0 / std::sqrt( double(hid) ), is3 = 1.0 / std::sqrt( 3.0 );
    // per pair: h2_a . h2_b, sw_a sw_b; per head: Q, K, V matrices of the head
    dpa3_grow(S.hh, nn2); dpa3_grow(S.ww, nn2); dpa3_grow(S.aq, size_t(n) * hid); dpa3_grow(S.ak, size_t(n) * hid);
    dpa3_grow(S.av, size_t(n) * g2d); dpa3_grow(S.ao, size_t(n) * g2d); dpa3_grow(S.am, nn2); dpa3_grow(S.mx, n);
    dpa3_gemm( false, true, n, n, 3, 1.0, g.h2.data(), 3, g.h2.data(), 3, 0.0, S.hh.data(), n );
    for(int a = 0; a < n; a++) for(int b = 0; b < n; b++) S.ww[size_t(a)*n + b] = g.w[a] * g.w[b];
    const double* hh = S.hh.data();
    const double* ww = S.ww.data();
    for(int h = 0; h < nh; h++)
    {
      double* L0 = w.L0.data() + size_t(h) * nn2;
      double* Ab = w.Ab.data() + size_t(h) * nn2;
      double* Am = S.am.data();
      for(int a = 0; a < n; a++)
      {
        const double* qa = w.qk.data() + size_t(a) * nq;
        const double* va = w.v.data() + size_t(a) * nv;
        double* Qa = S.aq.data() + size_t(a) * hid;
        double* Ka = S.ak.data() + size_t(a) * hid;
        double* Va = S.av.data() + size_t(a) * g2d;
        for(int c = 0; c < hid; c++) { Qa[c] = qa[c*2*nh + h]; Ka[c] = qa[c*2*nh + nh + h]; }
        for(int c = 0; c < g2d; c++) Va[c] = va[c*nh + h];
      }
      dpa3_gemm( false, true, n, n, hid, ish, S.aq.data(), hid, S.ak.data(), hid, 0.0, L0, n );
      if( M.gate ) for(size_t x = 0; x < nn2; x++) Ab[x] = ( L0[x] * hh[x] + M.shift ) * ww[x] - M.shift;
      else for(size_t x = 0; x < nn2; x++) Ab[x] = ( L0[x] + M.shift ) * ww[x] - M.shift;
      for(int a = 0; a < n; a++)
      {
        double* r = Ab + size_t(a) * n;
        double mx = npad > 0 ? -M.shift : -1e300;
        for(int b = 0; b < n; b++) mx = std::max( mx , r[b] );
        S.mx[a] = mx;
        for(int b = 0; b < n; b++) r[b] -= mx;
      }
      dpa2_exp( Ab, nn2 );
      for(int a = 0; a < n; a++)
      {
        double* r = Ab + size_t(a) * n;
        double z = npad > 0 ? npad * std::exp( -M.shift - S.mx[a] ) : 0.0;
#       pragma omp simd reduction(+:z)
        for(int b = 0; b < n; b++) z += r[b];
        const double iz = 1.0 / z;
        for(int b = 0; b < n; b++) r[b] *= iz;
      }
      // out[a][c nh + h] = sum_b A_ab v_b[c nh + h], A = Ab sw_a sw_b (h2_a . h2_b) / sqrt(3)
      for(size_t x = 0; x < nn2; x++) Am[x] = Ab[x] * ww[x] * hh[x] * is3;
      dpa3_gemm( false, false, n, g2d, n, 1.0, Am, n, S.av.data(), g2d, 0.0, S.ao.data(), g2d );
      for(int a = 0; a < n; a++)
      {
        double* oa = w.out.data() + size_t(a) * nv;
        const double* Oa = S.ao.data() + size_t(a) * g2d;
        for(int c = 0; c < g2d; c++) oa[c*nh + h] = Oa[c];
      }
    }
    T.mark(12);
    // head map, layer norm
    dpa3_grow(w.zh, size_t(n) * g2d); dpa3_grow(w.xh, size_t(n) * g2d); dpa3_grow(w.sig, n);
    for(int k = 0; k < n; k++) for(int c = 0; c < g2d; c++) w.zh[size_t(k)*g2d + c] = W[P.head_map.b + c];
    dpa3_gemm( false, false, n, g2d, nv, 1.0, w.out.data(), nv, W + P.head_map.w, g2d, 1.0, w.zh.data(), g2d );
    const double* r1 = W + P.g2_res[1];
    for(int k = 0; k < n; k++)
    {
      const double* z = w.zh.data() + size_t(k) * g2d;
      double mu = 0.0, var = 0.0;
#     pragma omp simd reduction(+:mu)
      for(int c = 0; c < g2d; c++) mu += z[c];
      mu /= g2d;
#     pragma omp simd reduction(+:var)
      for(int c = 0; c < g2d; c++) var += ( z[c] - mu ) * ( z[c] - mu );
      var /= g2d;
      const double sg = std::sqrt( var + P.lneps );
      w.sig[k] = sg;
      for(int c = 0; c < g2d; c++)
      {
        const double xh = ( z[c] - mu ) / sg;
        w.xh[size_t(k)*g2d + c] = xh;
        o.g2[size_t(k)*g2d + c] += r1[c] * ( xh * W[P.lnw + c] + W[P.lnb + c] );
      }
    }
    T.mark(13);
  }

  // VJP of layer l: S.dg1 (g1d), S.dg2 (n x g2d) of the outputs -> S.dg1_in, S.dg2_in (inputs), S.dgg1 (n x g1d, dE/dg1 of
  // the neighbors, assigned); geometry gradients added into S.gh2, S.gw if geom
  inline void dpa2_layer_vjp(const DPA2Model& M, int l, const double* Nl, DPA2Scratch& S, bool geom)
  {
    const DPA3Model& m = M.base;
    const DPA2Layer& P = M.layers[l];
    const DPA2Geom& g = S.g;
    const DPA2Work& w = S.L[l];
    const int n = g.n, g1d = M.g1d, g2d = M.g2d, ax = M.axis, nh = M.nh, hid = M.hid;
    const double* W = m.W.data();
    const double* g2 = w.g2.data();
    const double* dg1 = S.dg1.data();
    const double* dg2 = S.dg2.data();
    dpa3_grow(S.dg1_in, g1d); dpa3_grow(S.dg2_in, size_t(n) * g2d); dpa3_grow(S.dgg1, size_t(n) * g1d);
    double* dg1_in = S.dg1_in.data();
    double* dg2_in = S.dg2_in.data();
    double* dgg1 = S.dgg1.data();
    DPA2Timer T;
    for(int c = 0; c < g1d; c++) dg1_in[c] = dg1[c];
    for(size_t x = 0; x < size_t(n) * g2d; x++) dg2_in[x] = P.chnnl2 ? dg2[x] : dg2[x];
    std::fill( dgg1, dgg1 + size_t(n) * g1d, 0.0 );
    // residual weights of the g1 updates, in deepmd's order
    int ui = 0;
    auto res = [&](int u) { return W + P.g1_res[u]; };
    std::vector<double>& du = S.du; dpa3_grow(du, g1d);
    // g1 of the neighbors, symmetrization weights
    std::vector<double>& Xg = S.gx; dpa3_grow(Xg, size_t(n) * g1d);
    for(int k = 0; k < n; k++) { const double* nj = Nl + size_t( g.j[k] ) * g1d; std::copy( nj, nj + g1d, Xg.data() + size_t(k) * g1d ); }
    const double* X1 = Xg.data();
    // self mlp
    if( M.out_mlp )
    {
      const double* r = res(ui++);
      for(int c = 0; c < g1d; c++) du[c] = r[c] * dg1[c] * w.d_self[c];
      dpa2_lin_t( m, P.self_mlp, du.data(), dg1_in );
    }
    // conv as an update (out_conv)
    std::vector<double>& dc = S.dc; dpa3_grow(dc, std::max( g1d, g2d ));
    if( M.conv && M.out_conv )
    {
      const double* r = res(ui++);
      for(int c = 0; c < g1d; c++) dc[c] = r[c] * dg1[c];
    }
    // linear1
    const double* r1 = res(ui++);
    std::vector<double>& dz1 = S.dz1; dpa3_grow(dz1, g1d);
    for(int c = 0; c < g1d; c++) dz1[c] = r1[c] * dg1[c] * w.d1[c];
    const int dmlp = ( M.out_mlp ? 0 : g1d ) + ( M.conv && !M.out_conv ? g2d : 0 ) + ( M.grrg ? ax * g2d : 0 ) + ( M.drrd ? ax * g1d : 0 );
    std::vector<double>& dmlpv = S.dmlp; dmlpv.assign( dmlp, 0.0 );
    dpa2_lin_t( m, P.linear1, dz1.data(), dmlpv.data() );
    int off = 0;
    if( !M.out_mlp ) { for(int c = 0; c < g1d; c++) dg1_in[c] += dmlpv[c]; off += g1d; }
    if( M.conv && !M.out_conv ) { std::copy( dmlpv.begin() + off, dmlpv.begin() + off + g2d, dc.begin() ); off += g2d; }
    // conv backward: c = sum_k P_k x_k w_k / nnei
    if( M.conv )
    {
      const int pw = M.out_conv ? g1d : g2d;
      std::vector<double>& dP = S.tmp; dpa3_grow(dP, size_t(n) * pw);
      for(int k = 0; k < n; k++)
      {
        const double f = g.w[k] / M.nnei;
        const double* pk = w.P.data() + size_t(k) * pw;
        const double* xk = M.out_conv ? X1 + size_t(k) * g1d : g2 + size_t(k) * g2d;
        double* dx = M.out_conv ? dgg1 + size_t(k) * g1d : dg2_in + size_t(k) * g2d;
        double gwk = 0.0;
#       pragma omp simd reduction(+:gwk)
        for(int c = 0; c < pw; c++)
        {
          dP[size_t(k)*pw + c] = dc[c] * xk[c] * f;
          dx[c] += dc[c] * pk[c] * f;
          gwk += dc[c] * pk[c] * xk[c];
        }
        if( geom ) S.gw[k] += gwk / M.nnei;
      }
      // out_conv: P = g2 Wp -> dg2 += dP Wp^T ; else P = g1_j Wp -> dgg1 += dP Wp^T
      if( M.out_conv ) dpa3_gemm( false, true, n, g2d, g1d, 1.0, dP.data(), g1d, W + P.proj.w, g1d, 1.0, dg2_in, g2d );
      else dpa3_gemm( false, true, n, g1d, g2d, 1.0, dP.data(), g2d, W + P.proj.w, g2d, 1.0, dgg1, g1d );
    }
    // grrg / drrd backward
    auto grrg_bwd = [&](const double* Mh, int width, const double* dG, std::vector<double>& dM)
    {
      dM.assign( size_t(3) * width, 0.0 );
      for(int a = 0; a < ax; a++) for(int c = 0; c < width; c++)
      {
        const double gv = dG[size_t(a) * width + c] / 3.0;
        if( gv == 0.0 ) continue;
        for(int x = 0; x < 3; x++) { dM[x*width + c] += gv * Mh[x*width + a]; dM[x*width + a] += gv * Mh[x*width + c]; }
      }
    };
    // hg[3][width] = hw^T X: dX += hw dM, gh[n][3] = X dM^T
    double gh3[3 * 256];
    std::vector<double> ghv;
    auto sym_bwd = [&](const std::vector<double>& dM, int width, bool nei)
    {
      const double* X = nei ? X1 : g2;
      double* dX = nei ? dgg1 : dg2_in;
      dpa3_gemm( false, false, n, width, 3, 1.0, S.hw.data(), 3, dM.data(), width, 1.0, dX, width );
      if( !geom ) return;
      double* gh = size_t(n) <= 256 ? gh3 : ( ghv.resize( size_t(n) * 3 ), ghv.data() );
      dpa3_gemm( false, true, n, 3, width, 1.0, X, width, dM.data(), width, 0.0, gh, 3 );
      for(int k = 0; k < n; k++)
      {
        const double f = g.w[k] * M.s_sym;
        double gsw = 0.0;
        for(int x = 0; x < 3; x++) { S.gh2[3*k + x] += gh[3*k + x] * f; gsw += gh[3*k + x] * g.h2[3*k + x]; }
        S.gw[k] += gsw * M.s_sym;
      }
    };
    if( M.grrg ) { grrg_bwd( w.hg2.data(), g2d, dmlpv.data() + off, S.dM ); sym_bwd( S.dM, g2d, false ); off += ax * g2d; }
    if( M.drrd ) { grrg_bwd( w.hg1.data(), g1d, dmlpv.data() + off, S.dM ); sym_bwd( S.dM, g1d, true ); off += ax * g1d; }
    T.mark(14);
    if( !P.chnnl2 ) return;
    // g2 updates: linear2
    {
      const double* r0 = W + P.g2_res[0];
      std::vector<double>& dz = S.tmp; dpa3_grow(dz, size_t(n) * g2d);
      for(size_t x = 0; x < size_t(n) * g2d; x++) dz[x] = r0[x % g2d] * dg2[x] * w.d2[x];
      dpa3_gemm( false, true, n, g2d, g2d, 1.0, dz.data(), g2d, W + P.linear2.w, g2d, 1.0, dg2_in, g2d );
    }
    if( !P.attn ) return;
    // attention: layer norm, head map
    const int nq = hid * 2 * nh, nv = g2d * nh;
    const double* r1g = W + P.g2_res[1];
    std::vector<double>& dzh = S.tmp; dpa3_grow(dzh, size_t(n) * g2d);
    for(int k = 0; k < n; k++)
    {
      const double* xh = w.xh.data() + size_t(k) * g2d;
      double* dz = dzh.data() + size_t(k) * g2d;
      double m1 = 0.0, m2 = 0.0;
#     pragma omp simd reduction(+:m1,m2)
      for(int c = 0; c < g2d; c++) { const double gg = r1g[c] * dg2[size_t(k)*g2d + c] * W[P.lnw + c]; dz[c] = gg; m1 += gg; m2 += gg * xh[c]; }
      m1 /= g2d; m2 /= g2d;
      for(int c = 0; c < g2d; c++) dz[c] = ( dz[c] - m1 - xh[c] * m2 ) / w.sig[k];
    }
    std::vector<double>& dout = S.tmp2; dpa3_grow(dout, size_t(n) * nv);
    dpa3_gemm( false, true, n, nv, g2d, 1.0, dzh.data(), g2d, W + P.head_map.w, g2d, 0.0, dout.data(), nv );
    T.mark(15);
    std::vector<double>& dv = S.tmp3; dpa3_grow(dv, size_t(n) * nv);
    std::vector<double>& dqk = S.dqk; dpa3_grow(dqk, size_t(n) * nq);
    const size_t nn2 = size_t(n) * n;
    const double ish = 1.0 / std::sqrt( double(hid) ), is3 = 1.0 / std::sqrt( 3.0 );
    S.dS.assign( nn2, 0.0 ); S.dHH.assign( nn2, 0.0 );
    double* dS = S.dS.data();
    double* dHH = S.dHH.data();
    // Q, K, V, dE/dout of the head; ad: dE/dA, am: A, then dE/dL0 (with hh, ww of the forward)
    dpa3_grow(S.ad, nn2); dpa3_grow(S.dqh, size_t(n) * hid); dpa3_grow(S.dkh, size_t(n) * hid); dpa3_grow(S.dvh, size_t(n) * g2d);
    const double* hh = S.hh.data();
    const double* ww = S.ww.data();
    double* Am = S.am.data();
    double* dA = S.ad.data();
    double* Qh = S.aq.data();
    double* Kh = S.ak.data();
    double* Vh = S.av.data();
    double* Dh = S.ao.data();
    for(int h = 0; h < nh; h++)
    {
      const double* L0 = w.L0.data() + size_t(h) * nn2;
      const double* Ab = w.Ab.data() + size_t(h) * nn2;
      for(int a = 0; a < n; a++)
      {
        const double* qa = w.qk.data() + size_t(a) * nq;
        const double* va = w.v.data() + size_t(a) * nv;
        const double* da = dout.data() + size_t(a) * nv;
        for(int c = 0; c < hid; c++) { Qh[size_t(a)*hid + c] = qa[c*2*nh + h]; Kh[size_t(a)*hid + c] = qa[c*2*nh + nh + h]; }
        for(int c = 0; c < g2d; c++) { Vh[size_t(a)*g2d + c] = va[c*nh + h]; Dh[size_t(a)*g2d + c] = da[c*nh + h]; }
      }
      for(size_t x = 0; x < nn2; x++) Am[x] = Ab[x] * ww[x] * hh[x] * is3;
      dpa3_gemm( false, true, n, n, g2d, 1.0, Dh, g2d, Vh, g2d, 0.0, dA, n );                     // dE/dA
      dpa3_gemm( true, false, n, g2d, n, 1.0, Am, n, Dh, g2d, 0.0, S.dvh.data(), g2d );           // dE/dv of the head
      for(int a = 0; a < n; a++)
      {
        const size_t ra = size_t(a) * n;
        double t = 0.0;
#       pragma omp simd reduction(+:t)
        for(int b = 0; b < n; b++) t += dA[ra + b] * ww[ra + b] * hh[ra + b] * is3 * Ab[ra + b];
#       pragma omp simd
        for(int b = 0; b < n; b++)
        {
          const size_t x = ra + b;
          const double dAx = dA[x], ab = Ab[x], S2 = ww[x], hx = hh[x], l0 = L0[x];
          dS[x] += dAx * ab * hx * is3;
          dHH[x] += dAx * ab * S2 * is3;
          const double dlp = ab * ( dAx * S2 * hx * is3 - t );                   // dE/dL'
          const double lg = M.gate ? l0 * hx : l0;
          dS[x] += dlp * ( lg + M.shift );
          const double dlg = dlp * S2;
          dHH[x] += M.gate ? dlg * l0 : 0.0;
          dA[x] = ( M.gate ? dlg * hx : dlg ) * ish;                             // dE/dL0, from here on
        }
      }
      dpa3_gemm( false, false, n, hid, n, 1.0, dA, n, Kh, hid, 0.0, S.dqh.data(), hid );          // dE/dQ
      dpa3_gemm( true, false, n, hid, n, 1.0, dA, n, Qh, hid, 0.0, S.dkh.data(), hid );           // dE/dK
      for(int a = 0; a < n; a++)
      {
        double* dqa = dqk.data() + size_t(a) * nq;
        double* dva = dv.data() + size_t(a) * nv;
        for(int c = 0; c < hid; c++) { dqa[c*2*nh + h] = S.dqh[size_t(a)*hid + c]; dqa[c*2*nh + nh + h] = S.dkh[size_t(a)*hid + c]; }
        for(int c = 0; c < g2d; c++) dva[c*nh + h] = S.dvh[size_t(a)*g2d + c];
      }
    }
    T.mark(16);
    dpa3_gemm( false, true, n, g2d, nv, 1.0, dv.data(), nv, W + P.mapv.w, nv, 1.0, dg2_in, g2d );
    dpa3_gemm( false, true, n, g2d, nq, 1.0, dqk.data(), nq, W + P.mapqk.w, nq, 1.0, dg2_in, g2d );
    if( geom )
      for(int a = 0; a < n; a++) for(int b = 0; b < n; b++)
      {
        const double ds = dS[size_t(a)*n + b];
        S.gw[a] += ds * g.w[b]; S.gw[b] += ds * g.w[a];
        const double dh = dHH[size_t(a)*n + b];
        for(int x = 0; x < 3; x++) { S.gh2[3*a + x] += dh * g.h2[3*b + x]; S.gh2[3*b + x] += dh * g.h2[3*a + x]; }
      }
    T.mark(17);
  }

  // VJP of the initial g2 into the geometry gradients (S.gs)
  inline void dpa2_init_g2_vjp(const DPA2Model& M, const double* dg2, DPA2Scratch& S)
  {
    const DPA3Model& m = M.base;
    const double* W = m.W.data();
    for(int k = 0; k < S.g.n; k++)
    {
      double gx = 0.0;
      for(int c = 0; c < M.g2d; c++)
      {
        const double z = W[M.g2_embd.w + c] * S.g.s[k] + ( M.g2_embd.b != size_t(-1) ? W[M.g2_embd.b + c] : 0.0 );
        gx += dg2[size_t(k)*M.g2d + c] * dpa3_dact( m, z ) * W[M.g2_embd.w + c];
      }
      S.gs[k] += gx;
    }
  }

  inline void dpa2_geom_to_bonds(const DPA2Model& M, const DPA2Scratch& S, double* gb)
  {
    const DPA2Geom& g = S.g;
    const double* sd = M.dstd.data() + 4 * g.ti;
    for(int k = 0; k < g.n; k++)
    {
      const double* d = g.d.data() + 3 * k;
      const double r = g.r[k], q = g.q[k], w = g.w[k], dw = g.dw[k];
      double v[3] = { 0.0, 0.0, 0.0 }, radial = 0.0;
      for(int x = 0; x < 3; x++)
      {
        const double gh = S.gh2[3*k + x] / sd[1 + x];
        v[x] += gh * w / ( q * q );
        radial += gh * d[x] * ( dw / ( q * q ) - 2.0 * w / ( q * q * q ) );
      }
      radial += S.gw[k] * dw;
      radial += S.gs[k] / sd[0] * ( dw / q - w / ( q * q ) );
      for(int x = 0; x < 3; x++) gb[ 3 * g.bond[k] + x ] += v[x] + radial * d[x] / r;
    }
  }

  // ------------------------------------------------------------------ passes A / C (repinit, three-body; neighbor buffers)

  struct DPA2BlockOp
  {
    DeepmdAttenOp rep;
    DeepmdForceOp three;
    bool has_three = false;
    int dim_rep = 0, dim_3b = 0, dim0 = 0;
    const size_t* poff = nullptr;
    double* D0 = nullptr;              // pass A: [n_nodes][dim0] out
    const double* dD0 = nullptr;       // pass C: [n_nodes][dim0] in
    bool backward = false;

    template <class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells) const
    { FakeMat3d virial; ComputePairOptionalLocks<false> locks; FakeParticleLock lock_a; this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a); }
    template <class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3d &virial, CellParticlesT cells) const
    { ComputePairOptionalLocks<false> locks; FakeParticleLock lock_a; this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a); }
    template <class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    { FakeMat3d virial; this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a); }

    template <class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t nn, ComputeBufferT &buf, double &, double &fx, double &fy, double &fz, int type, Mat3dT &virial, CellParticlesT cells, GridCellLocksT, ParticleLockT &) const
    {
      static constexpr bool compute_virial = std::is_same_v<Mat3dT, Mat3d>;
      static thread_local DeepmdAttenWork W;
      static thread_local std::vector<ComputeBufferT> b3;
      if( b3.empty() ) b3.resize(1);
      const size_t i = poff[buf.cell] + buf.part;
      DPA2Timer T;
      rep.forward( nn, buf, type, W );
      T.mark(19);
      if( has_three ) { b3[0] = buf; three.forward_desc( nn, b3[0], type ); }
      T.mark(20);
      if( !backward )
      {
        double* out = D0 + i * dim0;
        for(int c = 0; c < dim_rep; c++) out[c] = W.desc[c];
        for(int c = 0; c < dim_3b; c++) out[dim_rep + c] = b3[0].ext.desc[c];
        return;
      }
      const double* g = dD0 + i * dim0;
      for(int c = 0; c < dim_rep; c++) W.ddesc[c] = g[c];
      Vec3d fi = { 0.0, 0.0, 0.0 };
      Mat3d vir = { 0.,0.,0., 0.,0.,0., 0.,0.,0. };
      rep.template backward<compute_virial>( buf, type, cells, W, fi, vir );
      T.mark(21);
      if( has_three )
      {
        auto& x = b3[0].ext;
        for(int c = 0; c < dim_3b; c++) x.desc[c] = g[dim_rep + c];
        three.template backward_desc<compute_virial>( b3[0], type, cells );
        T.mark(22);
        fi = fi + Vec3d{ x.fi[0], x.fi[1], x.fi[2] };
        if constexpr ( compute_virial ) vir = vir + Mat3d{ x.vir[0], x.vir[1], x.vir[2], x.vir[3], x.vir[4], x.vir[5], x.vir[6], x.vir[7], x.vir[8] };
      }
      atomic_add_contribution(fx, fi.x);
      atomic_add_contribution(fy, fi.y);
      atomic_add_contribution(fz, fi.z);
      if constexpr ( compute_virial ) { virial += vir; }
    }
  };

  // ------------------------------------------------------------------ repformer phases

  struct DPA2Buffers
  {
    DPA3Buffers c;                                   // CSR bonds, N / G per layer, forces, energies, locks
    std::vector<double> D0, dD0, z0;                 // [n_nodes][dim0], [n_nodes][dim0], [n_nodes][g1d] (pre-activation of g1_0)
    std::vector<int64_t> goff;                       // [inum + 1] g2 values per central atom
    std::vector< std::vector<double> > SG2;          // [layer][..] g2 inputs
    std::vector<double> dSG2;
  };

  // after pass A (B.D0 of the central atoms): repformers, fitting, reverse phases; B.dD0 of the central atoms, energies,
  // repformer forces and virials into the grid
  template<class GridT>
  inline void dpa2_compute( GridT& grid, const GridChunkNeighbors& chunk_neighbors, const Domain& domain, const DPA2Model& M,
                            const DeepmdModelView& fit, const GhostCommunicationScheme* scheme, MPI_Comm comm, DPA2Buffers& BB )
  {
    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;
    static constexpr double conv_energy = EXASTAMP_CONST_QUANTITY( 1. * eV );
    const DPA3Model& m = M.base;
    DPA3Buffers& B = BB.c;
    const size_t ncells = grid.number_of_cells();
    const size_t* poff = grid.cell_particle_offset_data();
    const size_t n_nodes = grid.number_of_particles();
    auto* cells = grid.cells();
    const IJK dims = grid.dimension();
    const Mat3d xform = domain.xform();
    const bool identity = domain.xform_is_identity();
    const double rc2 = M.rc * M.rc;
    const int g1d = M.g1d, g2d = M.g2d, NL = M.layers.size();
    DPA2Timer TW;

    // CSR bonds within the repformer cutoff
    B.type.resize( n_nodes );
    B.offset.assign( n_nodes + 1 , 0 );
    size_t ncentral = 0;
#   pragma omp parallel for schedule(dynamic) reduction(+:ncentral)
    for(size_t c = 0; c < ncells; c++)
    {
      const size_t np = cells[c].size();
      const auto* __restrict__ ty = cells[c][field::type];
      for(size_t p = 0; p < np; p++) B.type[poff[c] + p] = fit.species_type[ ty[p] ];
      if( grid.is_ghost_cell(c) ) continue;
      ncentral += np;
      const auto* __restrict__ rx = cells[c][field::rx];
      const auto* __restrict__ ry = cells[c][field::ry];
      const auto* __restrict__ rz = cells[c][field::rz];
      apply_cell_particle_neighbors( grid, chunk_neighbors, c, grid_index_to_ijk(dims, c), std::false_type{},
        [&](unsigned int p, size_t cb, unsigned int pb, size_t)
        {
          Vec3d d = { cells[cb][field::rx][pb] - rx[p], cells[cb][field::ry][pb] - ry[p], cells[cb][field::rz][pb] - rz[p] };
          if( !identity ) d = xform * d;
          if( norm2(d) < rc2 ) ++ B.offset[ poff[c] + p + 1 ];
        } );
    }
    B.ilist.resize( ncentral );
    { size_t k = 0; for(size_t c = 0; c < ncells; c++) { if( grid.is_ghost_cell(c) ) continue; for(size_t p = 0; p < cells[c].size(); p++) B.ilist[k++] = int( poff[c] + p ); } }
    for(size_t i = 0; i < n_nodes; i++) B.offset[i+1] += B.offset[i];
    const size_t nbonds = B.offset[n_nodes];
    B.nbh.resize( nbonds ); B.rij.resize( 3 * nbonds );
#   pragma omp parallel for schedule(dynamic)
    for(size_t c = 0; c < ncells; c++)
    {
      if( grid.is_ghost_cell(c) ) continue;
      const auto* __restrict__ rx = cells[c][field::rx];
      const auto* __restrict__ ry = cells[c][field::ry];
      const auto* __restrict__ rz = cells[c][field::rz];
      int64_t q = 0; unsigned int plast = ~0u;
      apply_cell_particle_neighbors( grid, chunk_neighbors, c, grid_index_to_ijk(dims, c), std::false_type{},
        [&](unsigned int p, size_t cb, unsigned int pb, size_t)
        {
          if( p != plast ) { q = B.offset[ poff[c] + p ]; plast = p; }
          Vec3d d = { cells[cb][field::rx][pb] - rx[p], cells[cb][field::ry][pb] - ry[p], cells[cb][field::rz][pb] - rz[p] };
          if( !identity ) d = xform * d;
          if( norm2(d) < rc2 ) { B.nbh[q] = int( poff[cb] + pb ); B.rij[3*q] = d.x; B.rij[3*q+1] = d.y; B.rij[3*q+2] = d.z; ++ q; }
        } );
    }
    const int inum = ncentral;
    const int64_t* offset = B.offset.data();
    const int* nbh = B.nbh.data();
    const double* rij = B.rij.data();
    auto exchange = [&](std::vector<double>& buf, bool to_ghost)
    {
      if( scheme == nullptr ) fatal_error() << "DPA2: no ghost communication scheme (a ghost update must run before deepmd_force)" << std::endl;
      GhostFlatBuffer fb;
      fb.m_data = buf.data(); fb.m_row_size = g1d; fb.m_row_stride = g1d; fb.m_nrows = n_nodes; fb.m_type = GhostFlatBuffer::ValueType::Double;
      if( to_ghost ) flat_buffer_ghost_exchange<true>( comm, *scheme, poff, fb, 323 );
      else flat_buffer_ghost_exchange<false>( comm, *scheme, poff, fb, 324 );
    };

    // g1 of layer 0 = act( T D0 (+ Tt tebd) ) of the central atoms, then the ghosts
    B.N.resize( NL + 1 ); B.G.resize( NL + 1 );
    for(auto& v : B.N) v.assign( n_nodes * g1d , 0.0 );
    BB.z0.assign( n_nodes * g1d, 0.0 );
#   pragma omp parallel for schedule(static)
    for(int q = 0; q < inum; q++)
    {
      const int i = B.ilist[q];
      double* z = BB.z0.data() + size_t(i) * g1d;
      const double* D0 = BB.D0.data() + size_t(i) * M.dim0;
      if( M.shape.valid() ) dpa3_lin( m, M.shape, D0, z, false, false ); else for(int c = 0; c < g1d; c++) z[c] = D0[c];
      if( M.tebd_tr.valid() ) dpa3_lin( m, M.tebd_tr, M.tebd.data() + size_t( B.type[i] ) * M.tdim, z, true, false );
      for(int c = 0; c < g1d; c++) B.N[0][size_t(i)*g1d + c] = dpa3_act( m, z[c] );
    }
    exchange( B.N[0], true );

    // per-atom g2 store
    std::vector<int> cn( inum );
#   pragma omp parallel
    {
      static thread_local DPA2Scratch S;
#     pragma omp for schedule(dynamic, 16)
      for(int q = 0; q < inum; q++) { dpa2_geometry( M, B.ilist[q], B.type[B.ilist[q]], offset, nbh, rij, S ); cn[q] = S.g.n; }
    }
    BB.goff.assign( inum + 1, 0 );
    for(int q = 0; q < inum; q++) BB.goff[q+1] = BB.goff[q] + int64_t(cn[q]) * g2d;
    BB.SG2.resize( NL );
    for(auto& v : BB.SG2) v.resize( BB.goff[inum] );
    BB.dSG2.assign( BB.goff[inum], 0.0 );
    auto inputs = [&](int q, int l, DPA2Scratch& S)
    {
      const int i = B.ilist[q];
      if( int(S.L.size()) < NL + 1 ) S.L.resize( NL + 1 );
      dpa2_geometry( M, i, B.type[i], offset, nbh, rij, S );
      DPA2Work& w = S.L[l];
      w.g1.assign( B.N[l].begin() + size_t(i) * g1d, B.N[l].begin() + size_t(i+1) * g1d );
      if( l == 0 )
      {
        dpa2_init_g2( M, S );
        std::copy( w.g2.begin(), w.g2.begin() + ( BB.goff[q+1] - BB.goff[q] ), BB.SG2[0].begin() + BB.goff[q] );
      }
      else w.g2.assign( BB.SG2[l].begin() + BB.goff[q], BB.SG2[l].begin() + BB.goff[q+1] );
    };

    TW.mark(1);
    // forward layers
    for(int l = 0; l < NL; l++)
    {
#     pragma omp parallel
      {
        static thread_local DPA2Scratch S;
#       pragma omp for schedule(dynamic, 1)
        for(int q = 0; q < inum; q++)
        {
          const int i = B.ilist[q];
          DPA2Timer T;
          inputs( q, l, S );
          T.mark(8);
          dpa2_layer_forward( M, l, B.N[l].data(), S );
          for(int c = 0; c < g1d; c++) B.N[l+1][size_t(i)*g1d + c] = S.L[l+1].g1[c];
          if( l + 1 < NL ) std::copy( S.L[l+1].g2.begin(), S.L[l+1].g2.begin() + ( BB.goff[q+1] - BB.goff[q] ), BB.SG2[l+1].begin() + BB.goff[q] );
        }
      }
      exchange( B.N[l+1], true );
    }

    TW.mark(2);
    // fitting
    B.eatom.assign( n_nodes, 0.0 );
    B.G[NL].assign( n_nodes * g1d, 0.0 );
    const int din = g1d + ( M.concat_tebd ? M.tdim : 0 );
#   pragma omp parallel
    {
      std::vector<double> x( din ), dx( din ), w0( std::max( din, DP_MAX_FIT_WIDTH ) ), w1( std::max( din, DP_MAX_FIT_WIDTH ) );
      auto wk = std::make_unique< DPMlpWork<DP_MAX_FIT_LAYERS, DP_MAX_FIT_WIDTH> >();
#     pragma omp for schedule(static)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q], ti = B.type[i];
        for(int c = 0; c < g1d; c++) x[c] = B.N[NL][size_t(i)*g1d + c];
        if( M.concat_tebd ) for(int c = 0; c < M.tdim; c++) x[g1d + c] = M.tebd[ size_t(ti) * M.tdim + c ];
        const DPNet& net = fit.fitting[ fit.fitting_mixed_types ? 0 : ti ];
        B.eatom[i] = dp_mlp_forward( net, fit.layers, fit.weights, x.data(), *wk )[0] + fit.energy_bias[ti];
        const double one = 1.0;
        dp_mlp_backward( net, fit.layers, fit.weights, x.data(), *wk, &one, dx.data(), w0.data(), w1.data() );
        for(int c = 0; c < g1d; c++) B.G[NL][size_t(i)*g1d + c] = dx[c];
      }
    }

    TW.mark(3);
    // reverse phases, down to dE/dg1 of layer 0
    B.gb.assign( 3 * nbonds, 0.0 );
    if( B.locks.size() != n_nodes ) B.locks.assign( n_nodes, 0 );
    for(int L = NL - 1; L >= 0; L--)
    {
      B.G[L].assign( n_nodes * g1d, 0.0 );
#     pragma omp parallel
      {
        static thread_local DPA2Scratch S;
#       pragma omp for schedule(dynamic, 1)
        for(int q = 0; q < inum; q++)
        {
          const int i = B.ilist[q];
          DPA2Timer T;
          inputs( q, L, S );
          T.mark(8);
          dpa2_layer_forward( M, L, B.N[L].data(), S );
          const int n = S.g.n;
          S.gh2.assign( size_t(n) * 3, 0.0 ); S.gw.assign( n, 0.0 ); S.gs.assign( n, 0.0 );
          dpa3_grow( S.dg1, g1d ); dpa3_grow( S.dg2, size_t(n) * g2d );
          for(int c = 0; c < g1d; c++) S.dg1[c] = B.G[L+1][size_t(i)*g1d + c];
          if( L < NL - 1 ) std::copy( BB.dSG2.begin() + BB.goff[q], BB.dSG2.begin() + BB.goff[q+1], S.dg2.begin() );
          else std::fill( S.dg2.begin(), S.dg2.begin() + size_t(n) * g2d, 0.0 );
          dpa2_layer_vjp( M, L, B.N[L].data(), S, true );
          T = DPA2Timer();
          {
            DPA3RowLock lk( &B.locks[i] );
            double* gi = B.G[L].data() + size_t(i) * g1d;
            for(int c = 0; c < g1d; c++) gi[c] += S.dg1_in[c];
          }
          for(int k = 0; k < n; k++)
          {
            const int j = S.g.j[k];
            DPA3RowLock lk( &B.locks[j] );
            double* gj = B.G[L].data() + size_t(j) * g1d;
            const double* dj = S.dgg1.data() + size_t(k) * g1d;
            for(int c = 0; c < g1d; c++) gj[c] += dj[c];
          }
          if( L > 0 ) std::copy( S.dg2_in.begin(), S.dg2_in.begin() + size_t(n) * g2d, BB.dSG2.begin() + BB.goff[q] );
          else dpa2_init_g2_vjp( M, S.dg2_in.data(), S );
          dpa2_geom_to_bonds( M, S, B.gb.data() );
          T.mark(18);
        }
      }
      exchange( B.G[L], false );
    }

    TW.mark(4);
    // dE/dD0 = T^T ( dE/dg1_0 act'(z0) ) of the central atoms
    BB.dD0.assign( n_nodes * M.dim0, 0.0 );
#   pragma omp parallel for schedule(static)
    for(int q = 0; q < inum; q++)
    {
      const int i = B.ilist[q];
      std::vector<double> dz( g1d );
      for(int c = 0; c < g1d; c++) dz[c] = B.G[0][size_t(i)*g1d + c] * dpa3_dact( m, BB.z0[size_t(i)*g1d + c] );
      double* dD = BB.dD0.data() + size_t(i) * M.dim0;
      if( M.shape.valid() ) dpa2_lin_t( m, M.shape, dz.data(), dD ); else for(int c = 0; c < g1d; c++) dD[c] = dz[c];
    }

    // repformer forces, per-atom virials, energies into the grid
    B.f.assign( 3 * n_nodes, 0.0 );
    B.vatom.assign( 9 * n_nodes, 0.0 );
#   pragma omp parallel for schedule(dynamic, 16)
    for(int q = 0; q < inum; q++)
    {
      const int i = B.ilist[q];
      double fi[3] = { 0, 0, 0 }, vi[9] = { 0,0,0, 0,0,0, 0,0,0 };
      for(int64_t b = offset[i]; b < offset[i+1]; b++)
      {
        const double* g = B.gb.data() + 3 * b;
        if( g[0] == 0.0 && g[1] == 0.0 && g[2] == 0.0 ) continue;
        const double* d = rij + 3 * b;
        for(int x = 0; x < 3; x++) fi[x] += g[x];
        for(int x = 0; x < 3; x++) for(int y = 0; y < 3; y++) vi[3*x + y] += -g[x] * d[y];
        DPA3RowLock l( &B.locks[ nbh[b] ] );
        for(int x = 0; x < 3; x++) B.f[3*size_t( nbh[b] ) + x] -= g[x];
      }
      DPA3RowLock l( &B.locks[i] );
      for(int x = 0; x < 3; x++) B.f[3*size_t(i) + x] += fi[x];
      for(int x = 0; x < 9; x++) B.vatom[9*size_t(i) + x] += vi[x];
    }
#   pragma omp parallel for schedule(static)
    for(size_t c = 0; c < ncells; c++)
    {
      const size_t np = cells[c].size();
      for(size_t p = 0; p < np; p++)
      {
        const size_t i = poff[c] + p;
        cells[c][field::fx][p] += B.f[3*i] * conv_energy;
        cells[c][field::fy][p] += B.f[3*i+1] * conv_energy;
        cells[c][field::fz][p] += B.f[3*i+2] * conv_energy;
        cells[c][field::ep][p] += B.eatom[i] * conv_energy;
        if constexpr ( has_virial_field )
        {
          const double* v = B.vatom.data() + 9*i;
          cells[c][field::virial][p] += Mat3d{ v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8] } * conv_energy;
        }
      }
    }
    TW.mark(6);
  }

}

namespace exanb
{
  template <> struct ComputePairTraits<exaStamp::DPA2BlockOp>
  {
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = false;
    static inline constexpr bool HasParticleContextStart = false;
    static inline constexpr bool HasParticleContextStop = false;
    static inline constexpr bool HasParticleContext = false;
    static inline constexpr bool RequiresNbhOptionalData = true;
    static inline constexpr bool ComputeEmptyNeighborhood = true;
  };
}
