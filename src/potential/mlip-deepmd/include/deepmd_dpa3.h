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

// DeePMD DPA3 (repflows, Zhang et al. 2025) evaluated natively on the host, a message-passing model: every layer reads
// the node features of the neighbors. Per layer (static neighbor selection; the dynamic one changes the normalizations):
//
//   node   n'  = n + r0 act(Wns n) + r1 act(Wsym [grrg(e sw, h2), grrg(n_j sw, h2)]) + r2 sum_j act(Wne [n, n_j, e]) sw * s_msg
//   edge   e'  = e + r0 act(Wes [n, n_j, e]) + r1 act(W2 sum_k a_sw_j a_sw_k act(W1 [a_jk, nc, ec_k, ec_j]) * s_ang)
//   angle  a'  = a + r0 act(Was [a_jk, nc, ec_k, ec_j])
//
// grrg(g, h) = M[:, :axis]^T M / 3, M = s_sym sum_j h2_j x g_j (3 x width); nc, ec: compressed node / edge features
// (first channels, a linear map, or all); edges: the nnei nearest atoms within rcut, angles: pairs of the a_sel nearest
// ones within a_rcut. Initial features: node = act(tebd[type]), edge = edge_embd(r or s), angle = angle_embd(cos / sqrt(pi)).
// Descriptor: the last node features (+ tebd), then the fitting network.
//
// Evaluation in phases (deepmd_force): node features of every layer N_l are kept for all atoms, ghosts updated after
// each forward layer (flat buffer ghost exchange); an atom's edge and angle features are recomputed from them when
// needed (no per-atom edge / angle storage). Reverse, for L = nlayers-1 .. 0: every atom recomputes its chain, runs the
// vector-Jacobian products from the top layer down to L with the node gradients of the layers above (already complete),
// adds dE/dn_i of layer L to its row of G_L and dE/dn_j to the neighbors' rows (row locks), and its geometry gradients
// of layer L (plus the initial embeddings for L = 0) to its bonds; G_L ghost rows are then folded into their owners.

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/particle_neighbors/chunk_neighbors_apply.h>
#include <exanb/mpi/flat_buffer_ghosts.h>
#include <exaStamp/unit_system.h>
#include <onika/log.h>

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include "deepmd_model.h"

#if defined(XSTAMP_DEEPMD_BLAS) && defined(__APPLE__)
extern "C" void vvexp(double* y, const double* x, const int* n);   // Accelerate vForce
#endif
#ifdef XSTAMP_DEEPMD_BLAS
extern "C" void cblas_dgemm(int layout, int transa, int transb, int M, int N, int K, double alpha, const double* A, int lda,
                            const double* B, int ldb, double beta, double* C, int ldc);
#endif

namespace exaStamp
{
  using namespace exanb;

  // C[M][N] = alpha op(A) op(B) + beta C, row major (op = transpose when ta / tb: A stored [K][M], B stored [N][K])
  inline void dpa3_gemm(bool ta, bool tb, int M, int N, int K, double alpha, const double* A, int lda, const double* B, int ldb,
                        double beta, double* C, int ldc)
  {
    if( M <= 0 || N <= 0 ) return;
#   ifdef XSTAMP_DEEPMD_BLAS
    cblas_dgemm( 101, ta ? 112 : 111, tb ? 112 : 111, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc );
#   else
    for(int i = 0; i < M; i++)
    {
      double* c = C + size_t(i) * ldc;
      if( beta == 0.0 ) for(int j = 0; j < N; j++) c[j] = 0.0; else if( beta != 1.0 ) for(int j = 0; j < N; j++) c[j] *= beta;
      for(int k = 0; k < K; k++)
      {
        const double a = alpha * ( ta ? A[ size_t(k) * lda + i ] : A[ size_t(i) * lda + k ] );
        if( a == 0.0 ) continue;
        if( !tb ) { const double* b = B + size_t(k) * ldb; for(int j = 0; j < N; j++) c[j] += a * b[j]; }
        else for(int j = 0; j < N; j++) c[j] += a * B[ size_t(j) * ldb + k ];
      }
    }
#   endif
  }

  // stage profile (DPA3_PROF=1): wall time per stage on OpenMP thread 0, printed at exit
  struct DPA3Prof
  {
    bool on = std::getenv("DPA3_PROF") != nullptr;
    double t[16] = {};
    ~DPA3Prof()
    {
      if( !on ) return;
      static const char* names[16] = { "geometry", "init features", "fwd node self+sym", "fwd edge products", "fwd node out", "fwd angle products",
                                       "fwd reduce+edge2", "fwd outputs", "vjp node+sym", "vjp edge", "vjp edge2", "vjp angle", "vjp compress",
                                       "init vjp + bonds", "fwd angle act", "driver other" };
      std::fprintf( stderr, "[DPA3] thread 0 profile, s:\n" );
      for(int k = 0; k < 16; k++) if( t[k] > 0 ) std::fprintf( stderr, "  %-20s %9.4f\n", names[k], t[k] );
    }
  };
  inline DPA3Prof& dpa3_prof() { static DPA3Prof p; return p; }
  struct DPA3Timer
  {
    bool on; double t0;
    inline DPA3Timer() : on( dpa3_prof().on && omp_get_thread_num() == 0 ), t0( on ? omp_get_wtime() : 0.0 ) {}
    inline void mark(int k) { if( on ) { const double t = omp_get_wtime(); dpa3_prof().t[k] += t - t0; t0 = t; } }
  };

  // ------------------------------------------------------------------ model

  struct DPA3Lin
  {
    size_t w = size_t(-1), b = size_t(-1);   // offsets in DPA3Model::W, w: nin x nout row-major, b: nout or none
    int nin = 0, nout = 0;
    inline bool valid() const { return w != size_t(-1); }
  };

  struct DPA3Layer
  {
    DPA3Lin node_self, node_sym, node_edge, edge_self, edge_angle1, edge_angle2, angle_self, compress_n, compress_e;
    size_t n_res[3] = { 0, 0, 0 }, e_res[2] = { 0, 0 }, a_res = 0;
  };

  enum DPA3Act : int { DPA3_SILUT = 0, DPA3_SILU, DPA3_TANH, DPA3_RELU, DPA3_GELU, DPA3_SIGMOID, DPA3_SOFTPLUS, DPA3_LINEAR };

  struct DPA3Model
  {
    std::vector<double> W;
    std::vector<DPA3Layer> layers;
    DPA3Lin edge_embd, angle_embd;
    std::vector<double> tebd, davg, dstd;          // [type][n_dim], [type][4], [type][4]
    int nt = 0, n_dim = 0, e_dim = 0, a_dim = 0, axis = 0, nnei = 0, a_sel = 0, nca = 0, eca = 0;
    int compress = 0;                              // 0: none (all channels), 1: split (first channels), 2: linear maps
    double rc = 0, rs = 0, arc = 0, ars = 0, prot = 0;
    bool exp_sw = false, use_dist = false, smooth_edge = true, concat_tebd = false;
    int act = DPA3_SILUT;
    double act_T = 3.0, act_slope = 0.0, act_const = 0.0;
    double s_sym = 1.0, s_msg = 1.0, s_ang = 1.0;  // normalizations of grrg, the node-edge message and the edge-angle sum
  };

  // activation and its derivative
  inline double dpa3_sig(double z) { return 1.0 / ( 1.0 + std::exp(-z) ); }
  inline double dpa3_act(const DPA3Model& m, double x)
  {
    switch( m.act )
    {
      case DPA3_SILUT : return x < m.act_T ? x * dpa3_sig(x) : std::tanh( m.act_slope * ( x - m.act_T ) ) + m.act_const;
      case DPA3_SILU : return x * dpa3_sig(x);
      case DPA3_TANH : return std::tanh(x);
      case DPA3_RELU : return x > 0.0 ? x : 0.0;
      case DPA3_GELU : { const double k = 0.7978845608028654; return 0.5 * x * ( 1.0 + std::tanh( k * ( x + 0.044715 * x*x*x ) ) ); }
      case DPA3_SIGMOID : return dpa3_sig(x);
      case DPA3_SOFTPLUS : return x > 30.0 ? x : std::log1p( std::exp(x) );
      default : return x;
    }
  }
  inline double dpa3_dact(const DPA3Model& m, double x)
  {
    switch( m.act )
    {
      case DPA3_SILUT :
        if( x < m.act_T ) { const double s = dpa3_sig(x); return s + x * s * ( 1.0 - s ); }
        else { const double t = std::tanh( m.act_slope * ( x - m.act_T ) ); return m.act_slope * ( 1.0 - t*t ); }
      case DPA3_SILU : { const double s = dpa3_sig(x); return s + x * s * ( 1.0 - s ); }
      case DPA3_TANH : { const double t = std::tanh(x); return 1.0 - t*t; }
      case DPA3_RELU : return x > 0.0 ? 1.0 : 0.0;
      case DPA3_GELU :
      {
        const double k = 0.7978845608028654, t = std::tanh( k * ( x + 0.044715 * x*x*x ) );
        return 0.5 * ( 1.0 + t ) + 0.5 * x * ( 1.0 - t*t ) * k * ( 1.0 + 3.0 * 0.044715 * x*x );
      }
      case DPA3_SIGMOID : { const double s = dpa3_sig(x); return s * ( 1.0 - s ); }
      case DPA3_SOFTPLUS : return dpa3_sig(x);
      default : return 1.0;
    }
  }

  // activations and derivatives of an array (one exponential per value, batched with Accelerate vvexp on macOS)
  inline void dpa3_actd(const DPA3Model& m, const double* z, size_t n, double* a, double* d)
  {
    if( m.act == DPA3_SILUT || m.act == DPA3_SILU || m.act == DPA3_SIGMOID )
    {
      static thread_local std::vector<double> e;
      if( e.size() < n ) e.resize(n);
      for(size_t k = 0; k < n; k++) e[k] = -z[k];
#     if defined(XSTAMP_DEEPMD_BLAS) && defined(__APPLE__)
      { const int nn = int(n); vvexp( e.data(), e.data(), &nn ); }
#     else
      for(size_t k = 0; k < n; k++) e[k] = std::exp( e[k] );
#     endif
      if( m.act == DPA3_SIGMOID )
      {
        for(size_t k = 0; k < n; k++) { const double sg = 1.0 / ( 1.0 + e[k] ); if( a ) a[k] = sg; if( d ) d[k] = sg * ( 1.0 - sg ); }
        return;
      }
      // silu (branch-free loops), then the values past the silut threshold
      if( a ) for(size_t k = 0; k < n; k++) a[k] = z[k] / ( 1.0 + e[k] );
      if( d ) for(size_t k = 0; k < n; k++) { const double sg = 1.0 / ( 1.0 + e[k] ); d[k] = sg + z[k] * sg * ( 1.0 - sg ); }
      if( m.act == DPA3_SILUT )
        for(size_t k = 0; k < n; k++)
          if( z[k] >= m.act_T )
          {
            const double t = std::tanh( m.act_slope * ( z[k] - m.act_T ) );
            if( a ) a[k] = t + m.act_const;
            if( d ) d[k] = m.act_slope * ( 1.0 - t*t );
          }
      return;
    }
    for(size_t k = 0; k < n; k++) { if( a ) a[k] = dpa3_act( m, z[k] ); if( d ) d[k] = dpa3_dact( m, z[k] ); }
  }

  // switch functions (polynomial or exponential) and their radial derivative
  inline void dpa3_switch(const DPA3Model& m, double r, double rs, double rc, double& w, double& dw)
  {
    if( m.exp_sw )
    {
      // w = exp( -exp( 20/rs ( min(r, rc) - rs ) ) )
      if( r >= rc ) { const double e = std::exp( 20.0 / rs * ( rc - rs ) ); w = std::exp(-e); dw = 0.0; return; }
      const double a = 20.0 / rs, e = std::exp( a * ( r - rs ) );
      w = std::exp(-e);
      dw = -w * e * a;
      return;
    }
    if( r <= rs ) { w = 1.0; dw = 0.0; return; }
    if( r >= rc ) { w = 0.0; dw = 0.0; return; }
    const double u = ( r - rs ) / ( rc - rs ), u2 = u * u;
    w = u2 * u * ( -6.0 * u2 + 15.0 * u - 10.0 ) + 1.0;
    dw = u2 * ( -30.0 * u2 + 60.0 * u - 30.0 ) / ( rc - rs );
  }

  // y[nout] (+)= x[nin] W (+ b)
  inline void dpa3_lin(const DPA3Model& m, const DPA3Lin& L, const double* x, double* y, bool accumulate = false, bool bias = true)
  {
    const double* w = m.W.data() + L.w;
    if( !accumulate ) for(int o = 0; o < L.nout; o++) y[o] = ( bias && L.b != size_t(-1) ) ? m.W[L.b + o] : 0.0;
    else if( bias && L.b != size_t(-1) ) for(int o = 0; o < L.nout; o++) y[o] += m.W[L.b + o];
    for(int i = 0; i < L.nin; i++)
    {
      const double xi = x[i];
      if( xi == 0.0 ) continue;
      const double* wi = w + size_t(i) * L.nout;
      for(int o = 0; o < L.nout; o++) y[o] += xi * wi[o];
    }
  }
  // rows r0 .. r0 + n of W times x: y[nout] += x[n] W[r0:r0+n]
  inline void dpa3_lin_rows(const DPA3Model& m, const DPA3Lin& L, int r0, int n, const double* x, double* y)
  {
    const double* w = m.W.data() + L.w + size_t(r0) * L.nout;
    for(int i = 0; i < n; i++)
    {
      const double xi = x[i];
      if( xi == 0.0 ) continue;
      const double* wi = w + size_t(i) * L.nout;
      for(int o = 0; o < L.nout; o++) y[o] += xi * wi[o];
    }
  }
  // transposed: dx[n] += W[r0:r0+n] dy[nout]
  inline void dpa3_lin_rows_t(const DPA3Model& m, const DPA3Lin& L, int r0, int n, const double* dy, double* dx)
  {
    const double* w = m.W.data() + L.w + size_t(r0) * L.nout;
    for(int i = 0; i < n; i++)
    {
      const double* wi = w + size_t(i) * L.nout;
      double s = 0.0;
      for(int o = 0; o < L.nout; o++) s += wi[o] * dy[o];
      dx[i] += s;
    }
  }

  // ------------------------------------------------------------------ one atom

  // geometry of one central atom: edges (nearest nnei within rcut), angle edges (first a_sel within a_rcut)
  struct DPA3Geom
  {
    int n = 0, na = 0, ti = 0;
    std::vector<int> j;          // [n] flat index of the neighbor
    std::vector<int64_t> bond;   // [n] its bond in the CSR list
    std::vector<double> d, r, q, w, dw, h2, s;   // [n][3], [n], [n], [n], [n], [n][3], [n] (normalized radial row)
    std::vector<double> u, aw, adw, cosv;        // angle edges (the first na edges): [na][3] d/(r+1e-6), switch, derivative, [na][na]
  };

  // features and intermediates of one layer of one atom
  struct DPA3Work
  {
    std::vector<double> n, e, a;                 // inputs: [n_dim], [n][e_dim], [na][na][a_dim]
    std::vector<double> z_self, z_sym, sym_in, hg_e, hg_n, z_msg, z_es, nc, ec, z_ea, red, full, z_e2, z_as;
    // activations and their derivatives of the pre-activations above
    std::vector<double> a_self, d_self, a_sym, d_sym, a_msg, d_msg, a_es, d_es, a_ea, d_ea, a_e2, d_e2, a_as, d_as;
  };

  // gradients of the geometry of one atom (one phase)
  struct DPA3GeomGrad
  {
    std::vector<double> gh2, gw, gr, gs, gu, gaw;   // [n][3], [n], [n], [n], [na][3], [na]
    void reset(int n, int na)
    {
      gh2.assign( size_t(n) * 3, 0.0 ); gw.assign( n, 0.0 ); gr.assign( n, 0.0 ); gs.assign( n, 0.0 );
      gu.assign( size_t(na) * 3, 0.0 ); gaw.assign( na, 0.0 );
    }
  };

  struct DPA3AtomScratch
  {
    DPA3Geom g;
    std::vector<DPA3Work> L;                       // layer inputs and intermediates, L[nlayers] holds the outputs
    std::vector<int> order;
    std::vector<double> X, dX, dZm, dZe, dZa1, dZa2, rs1, cs1, rs2, cs2;   // matrix product operands
    // backward
    std::vector<double> dn, de, da, dn_in, de_in, da_in, dz, dz2, dnc, dec, dnj, tmp, dsym, dhg_e, dhg_n;
    DPA3GeomGrad gg;
  };

  inline void dpa3_grow(std::vector<double>& v, size_t s) { if( v.size() < s ) v.resize(s); }

  // bonds of central atom i (CSR: offset, nbh, rij), its edge and angle geometry
  inline void dpa3_geometry(const DPA3Model& m, int i, int ti, const int64_t* offset, const int* nbh, const double* rij, DPA3AtomScratch& S)
  {
    DPA3Geom& g = S.g;
    g.ti = ti;
    const int64_t b0 = offset[i], b1 = offset[i+1];
    const double rc2 = m.rc * m.rc;
    S.order.clear();
    std::vector<double> d2s;
    for(int64_t b = b0; b < b1; b++)
    {
      const double* x = rij + 3 * b;
      const double d2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
      if( d2 < rc2 && d2 > 0.0 ) S.order.push_back( int( b - b0 ) );
    }
    // nearest first (deepmd's sorted neighbor list), the nnei nearest
    std::stable_sort( S.order.begin(), S.order.end(), [&](int x, int y)
    {
      const double* a = rij + 3 * ( b0 + x ); const double* c = rij + 3 * ( b0 + y );
      return a[0]*a[0] + a[1]*a[1] + a[2]*a[2] < c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
    } );
    if( int(S.order.size()) > m.nnei ) S.order.resize( m.nnei );
    const int n = S.order.size();
    g.n = n;
    g.j.resize(n); g.bond.resize(n);
    dpa3_grow(g.d, size_t(n)*3); dpa3_grow(g.r, n); dpa3_grow(g.q, n); dpa3_grow(g.w, n); dpa3_grow(g.dw, n); dpa3_grow(g.h2, size_t(n)*3); dpa3_grow(g.s, n);
    const double* avg = m.davg.data() + 4 * ti;
    const double* sd = m.dstd.data() + 4 * ti;
    for(int k = 0; k < n; k++)
    {
      const int64_t b = b0 + S.order[k];
      g.bond[k] = b;
      g.j[k] = nbh[b];
      const double* x = rij + 3 * b;
      for(int c = 0; c < 3; c++) g.d[3*k + c] = x[c];
      const double r = std::sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
      g.r[k] = r;
      g.q[k] = r + m.prot;
      dpa3_switch( m, r, m.rs, m.rc, g.w[k], g.dw[k] );
      const double qq = g.q[k];
      for(int c = 0; c < 3; c++) g.h2[3*k + c] = ( g.w[k] * x[c] / ( qq * qq ) - avg[1 + c] ) / sd[1 + c];
      g.s[k] = ( g.w[k] / qq - avg[0] ) / sd[0];
    }
    // angle edges: the first a_sel edges within a_rcut (edges are sorted: a prefix)
    int na = 0;
    while( na < n && na < m.a_sel && g.r[na] < m.arc ) na++;
    // (an edge past a_rcut ends the prefix: the following ones are farther)
    g.na = na;
    dpa3_grow(g.u, size_t(na)*3); dpa3_grow(g.aw, na); dpa3_grow(g.adw, na); dpa3_grow(g.cosv, size_t(na)*na);
    for(int k = 0; k < na; k++)
    {
      dpa3_switch( m, g.r[k], m.ars, m.arc, g.aw[k], g.adw[k] );
      for(int c = 0; c < 3; c++) g.u[3*k + c] = g.d[3*k + c] / ( g.r[k] + 1e-6 );
    }
    for(int a = 0; a < na; a++) for(int b = 0; b < na; b++)
      g.cosv[size_t(a)*na + b] = ( g.u[3*a]*g.u[3*b] + g.u[3*a+1]*g.u[3*b+1] + g.u[3*a+2]*g.u[3*b+2] ) * ( 1.0 - 1e-6 );
  }

  // initial edge and angle features of the atom into S.L[0] (node: from N0)
  inline void dpa3_init_features(const DPA3Model& m, const double* n0, DPA3AtomScratch& S)
  {
    const DPA3Geom& g = S.g;
    DPA3Work& w = S.L[0];
    const int ed = m.e_dim, ad = m.a_dim, n = g.n, na = g.na;
    w.n.assign( n0, n0 + m.n_dim );
    dpa3_grow(w.e, size_t(n) * ed);
    for(int k = 0; k < n; k++)
    {
      const double x = m.use_dist ? g.r[k] : g.s[k];
      double* e = w.e.data() + size_t(k) * ed;
      dpa3_lin( m, m.edge_embd, &x, e );
      if( !m.use_dist ) for(int c = 0; c < ed; c++) e[c] = dpa3_act( m, e[c] );
    }
    dpa3_grow(w.a, size_t(na) * na * ad);
    const double isp = 1.0 / std::sqrt( M_PI );
    for(int p = 0; p < na * na; p++)
    {
      const double x = g.cosv[p] * isp;
      dpa3_lin( m, m.angle_embd, &x, w.a.data() + size_t(p) * ad );
    }
  }

  // grrg(M) = M[:, :axis]^T M / 3, M: 3 x width
  inline void dpa3_grrg(const double* M, int width, int axis, double* out)
  {
    for(int a = 0; a < axis; a++) for(int c = 0; c < width; c++)
      out[size_t(a) * width + c] = ( M[a] * M[c] + M[width + a] * M[width + c] + M[2*width + a] * M[2*width + c] ) / 3.0;
  }

  // forward of layer l for the atom: S.L[l] (inputs) -> S.L[l+1] (n, e, a), intermediates kept in S.L[l];
  // Nl: node features of layer l of all atoms (ghosts included). for_vjp: only what the vector-Jacobian product reads
  // (activation derivatives, no outputs); else the outputs only (no derivatives). mode 0: outputs, 1: for the VJP,
  // 2: both
  inline void dpa3_layer_forward(const DPA3Model& m, int l, const double* Nl, DPA3AtomScratch& S, int mode = 0)
  {
    const bool fv = mode >= 1, outs = mode != 1;
    const bool last = l == int( m.layers.size() ) - 1;    // last layer: only the node output is used
    const DPA3Layer& P = m.layers[l];
    const DPA3Geom& g = S.g;
    DPA3Work& w = S.L[l];
    DPA3Work& o = S.L[l+1];
    const int nd = m.n_dim, ed = m.e_dim, ad = m.a_dim, n = g.n, na = g.na, ax = m.axis, nca = m.nca, eca = m.eca;
    const double* nvec = w.n.data();
    const double* W = m.W.data();
    DPA3Timer T;
    // node self
    dpa3_grow(w.z_self, nd); dpa3_grow(w.a_self, nd); dpa3_grow(w.d_self, nd);
    dpa3_lin( m, P.node_self, nvec, w.z_self.data() );
    dpa3_actd( m, w.z_self.data(), nd, outs ? w.a_self.data() : nullptr, fv ? w.d_self.data() : nullptr );
    // symmetrization: hg_e, hg_n (3 x width), grrg, linear
    dpa3_grow(w.hg_e, size_t(3) * ed); dpa3_grow(w.hg_n, size_t(3) * nd);
    std::fill( w.hg_e.begin(), w.hg_e.begin() + 3 * ed, 0.0 );
    std::fill( w.hg_n.begin(), w.hg_n.begin() + 3 * nd, 0.0 );
    for(int k = 0; k < n; k++)
    {
      const double* h = g.h2.data() + 3 * k;
      const double sw = g.w[k] * m.s_sym;
      const double* e = w.e.data() + size_t(k) * ed;
      const double* nj = Nl + size_t( g.j[k] ) * nd;
      for(int x = 0; x < 3; x++)
      {
        const double hs = h[x] * sw;
        for(int c = 0; c < ed; c++) w.hg_e[x*ed + c] += hs * e[c];
        for(int c = 0; c < nd; c++) w.hg_n[x*nd + c] += hs * nj[c];
      }
    }
    dpa3_grow(w.sym_in, size_t(ax) * ( ed + nd ));
    dpa3_grrg( w.hg_e.data(), ed, ax, w.sym_in.data() );
    dpa3_grrg( w.hg_n.data(), nd, ax, w.sym_in.data() + size_t(ax) * ed );
    dpa3_grow(w.z_sym, nd); dpa3_grow(w.a_sym, nd); dpa3_grow(w.d_sym, nd);
    dpa3_lin( m, P.node_sym, w.sym_in.data(), w.z_sym.data() );
    dpa3_actd( m, w.z_sym.data(), nd, outs ? w.a_sym.data() : nullptr, fv ? w.d_sym.data() : nullptr );
    T.mark(2);
    // node-edge message and edge self update: z = W[0:nd] n + W[nd:2nd] n_j + W[2nd:] e + b
    std::vector<double>& base_m = S.tmp; dpa3_grow(base_m, size_t(nd) + ed);
    for(int c = 0; c < nd; c++) base_m[c] = P.node_edge.b != size_t(-1) ? W[P.node_edge.b + c] : 0.0;
    dpa3_lin_rows( m, P.node_edge, 0, nd, nvec, base_m.data() );
    double* base_e = base_m.data() + nd;
    for(int c = 0; c < ed; c++) base_e[c] = P.edge_self.b != size_t(-1) ? W[P.edge_self.b + c] : 0.0;
    dpa3_lin_rows( m, P.edge_self, 0, nd, nvec, base_e );
    // X = [ n_j | e ] per edge, z = X W[nd:] + base (two matrix products)
    dpa3_grow(w.z_msg, size_t(n) * nd); dpa3_grow(w.z_es, size_t(n) * ed);
    std::vector<double>& X = S.X; dpa3_grow(X, size_t(n) * ( nd + ed ));
    for(int k = 0; k < n; k++)
    {
      const double* nj = Nl + size_t( g.j[k] ) * nd;
      double* xk = X.data() + size_t(k) * ( nd + ed );
      for(int c = 0; c < nd; c++) xk[c] = nj[c];
      for(int c = 0; c < ed; c++) xk[nd + c] = w.e[size_t(k)*ed + c];
      for(int c = 0; c < nd; c++) w.z_msg[size_t(k)*nd + c] = base_m[c];
      for(int c = 0; c < ed; c++) w.z_es[size_t(k)*ed + c] = base_e[c];
    }
    dpa3_gemm( false, false, n, nd, nd + ed, 1.0, X.data(), nd + ed, W + P.node_edge.w + size_t(nd) * nd, nd, 1.0, w.z_msg.data(), nd );
    if( !last ) dpa3_gemm( false, false, n, ed, nd + ed, 1.0, X.data(), nd + ed, W + P.edge_self.w + size_t(nd) * ed, ed, 1.0, w.z_es.data(), ed );
    dpa3_grow(w.a_msg, size_t(n) * nd); dpa3_grow(w.d_msg, size_t(n) * nd); dpa3_grow(w.a_es, size_t(n) * ed); dpa3_grow(w.d_es, size_t(n) * ed);
    dpa3_actd( m, w.z_msg.data(), size_t(n) * nd, w.a_msg.data(), fv ? w.d_msg.data() : nullptr );   // the VJP reads a_msg
    if( !last ) dpa3_actd( m, w.z_es.data(), size_t(n) * ed, outs ? w.a_es.data() : nullptr, fv ? w.d_es.data() : nullptr );
    T.mark(3);
    if( outs )
    {
    // node output
    const double *rn0 = W + P.n_res[0], *rn1 = W + P.n_res[1], *rn2 = W + P.n_res[2];
    o.n.resize(nd);
    std::vector<double> msg( nd, 0.0 );
    for(int k = 0; k < n; k++)
    {
      const double* am = w.a_msg.data() + size_t(k) * nd;
      const double wk = g.w[k];
      for(int c = 0; c < nd; c++) msg[c] += am[c] * wk;
    }
    for(int c = 0; c < nd; c++) o.n[c] = nvec[c] + rn0[c] * w.a_self[c] + rn1[c] * w.a_sym[c] + rn2[c] * msg[c] * m.s_msg;
    }
    T.mark(4);
    if( last ) return;
    // compressed features for the angles
    dpa3_grow(w.nc, nca); dpa3_grow(w.ec, size_t(na) * eca);
    if( m.compress == 2 ) dpa3_lin( m, P.compress_n, nvec, w.nc.data(), false, false );
    else for(int c = 0; c < nca; c++) w.nc[c] = nvec[c];
    for(int k = 0; k < na; k++)
    {
      const double* e = w.e.data() + size_t(k) * ed;
      if( m.compress == 2 ) dpa3_lin( m, P.compress_e, e, w.ec.data() + size_t(k) * eca, false, false );
      else for(int c = 0; c < eca; c++) w.ec[size_t(k)*eca + c] = e[c];
    }
    // angle infos z[j][k] = W[0:ad] a_jk + W[ad:ad+nca] nc + W[.. +eca] ec_k + W[.. +eca] ec_j + b, for the edge-angle and
    // angle-self linears
    auto angle_z = [&](const DPA3Lin& Lw, int nout, std::vector<double>& z)
    {
      dpa3_grow(z, size_t(na) * na * nout);
      std::vector<double> bn( nout ), pk( size_t(na) * nout, 0.0 ), pj( size_t(na) * nout, 0.0 );
      for(int c = 0; c < nout; c++) bn[c] = Lw.b != size_t(-1) ? W[Lw.b + c] : 0.0;
      dpa3_lin_rows( m, Lw, ad, nca, w.nc.data(), bn.data() );
      dpa3_gemm( false, false, na, nout, eca, 1.0, w.ec.data(), eca, W + Lw.w + size_t( ad + nca ) * nout, nout, 0.0, pk.data(), nout );
      dpa3_gemm( false, false, na, nout, eca, 1.0, w.ec.data(), eca, W + Lw.w + size_t( ad + nca + eca ) * nout, nout, 0.0, pj.data(), nout );
      for(int a = 0; a < na; a++) for(int b = 0; b < na; b++)
      {
        double* zz = z.data() + ( size_t(a) * na + b ) * nout;
        for(int c = 0; c < nout; c++) zz[c] = bn[c] + pk[size_t(b)*nout + c] + pj[size_t(a)*nout + c];
      }
      dpa3_gemm( false, false, na * na, nout, ad, 1.0, w.a.data(), ad, W + Lw.w, nout, 1.0, z.data(), nout );
    };
    angle_z( P.edge_angle1, ed, w.z_ea );
    angle_z( P.angle_self, ad, w.z_as );
    T.mark(5);
    {
      const size_t np2 = size_t(na) * na;
      dpa3_grow(w.a_ea, np2 * ed); dpa3_grow(w.d_ea, np2 * ed); dpa3_grow(w.a_as, np2 * ad); dpa3_grow(w.d_as, np2 * ad);
      dpa3_actd( m, w.z_ea.data(), np2 * ed, w.a_ea.data(), fv ? w.d_ea.data() : nullptr );   // the reduction and the VJP read a_ea
      if( fv ) dpa3_actd( m, w.z_as.data(), np2 * ad, nullptr, w.d_as.data() );
    }
    T.mark(14);
    // edge-angle reduction red_j = s_ang sum_k aw_j aw_k act(z_ea[j][k]), then the full edge input of linear2
    dpa3_grow(w.red, size_t(na) * ed);
    for(int a = 0; a < na; a++)
    {
      double* rr = w.red.data() + size_t(a) * ed;
      for(int c = 0; c < ed; c++) rr[c] = 0.0;
      for(int b = 0; b < na; b++)
      {
        const double f = g.aw[a] * g.aw[b] * m.s_ang;
        const double* aa = w.a_ea.data() + ( size_t(a) * na + b ) * ed;
        for(int c = 0; c < ed; c++) rr[c] += f * aa[c];
      }
    }
    dpa3_grow(w.full, size_t(n) * ed); dpa3_grow(w.z_e2, size_t(n) * ed);
    for(int k = 0; k < n; k++)
    {
      double* f = w.full.data() + size_t(k) * ed;
      if( k < na ) for(int c = 0; c < ed; c++) f[c] = w.red[size_t(k)*ed + c];
      else if( !m.smooth_edge ) for(int c = 0; c < ed; c++) f[c] = w.e[size_t(k)*ed + c];
      else for(int c = 0; c < ed; c++) f[c] = 0.0;
      double* z2 = w.z_e2.data() + size_t(k) * ed;
      for(int c = 0; c < ed; c++) z2[c] = P.edge_angle2.b != size_t(-1) ? W[P.edge_angle2.b + c] : 0.0;
    }
    dpa3_gemm( false, false, n, ed, ed, 1.0, w.full.data(), ed, W + P.edge_angle2.w, ed, 1.0, w.z_e2.data(), ed );
    dpa3_grow(w.a_e2, size_t(n) * ed); dpa3_grow(w.d_e2, size_t(n) * ed);
    dpa3_actd( m, w.z_e2.data(), size_t(n) * ed, outs ? w.a_e2.data() : nullptr, fv ? w.d_e2.data() : nullptr );
    if( !outs ) { T.mark(6); return; }
    T.mark(6);
    // edge and angle outputs
    const double *re0 = W + P.e_res[0], *re1 = W + P.e_res[1], *ra0 = W + P.a_res;
    dpa3_grow(o.e, size_t(n) * ed);
    for(int k = 0; k < n; k++)
      for(int c = 0; c < ed; c++)
      {
        const size_t x = size_t(k) * ed + c;
        o.e[x] = w.e[x] + re0[c] * w.a_es[x] + re1[c] * w.a_e2[x];
      }
    dpa3_grow(o.a, size_t(na) * na * ad);
    // angle outputs: a' = a + r0 act(z_as), activations computed here (not needed by the VJP)
    {
      const size_t np2 = size_t(na) * na * ad;
      dpa3_actd( m, w.z_as.data(), np2, w.a_as.data(), nullptr );
      for(size_t x = 0; x < np2; x++) o.a[x] = w.a[x] + ra0[x % ad] * w.a_as[x];
    }
    T.mark(7);
  }

  // vector-Jacobian product of layer l: (dn, de, da) of its outputs -> (dn_in, de_in, da_in) of its inputs (assigned),
  // dE/dn_j added into dnj [n][n_dim] (assigned), geometry gradients added into S.gg if geom
  inline void dpa3_layer_vjp(const DPA3Model& m, int l, const double* Nl, DPA3AtomScratch& S, bool geom)
  {
    const DPA3Layer& P = m.layers[l];
    const DPA3Geom& g = S.g;
    const DPA3Work& w = S.L[l];
    const int nd = m.n_dim, ed = m.e_dim, ad = m.a_dim, n = g.n, na = g.na, ax = m.axis, nca = m.nca, eca = m.eca;
    const double* W = m.W.data();
    const bool last = l == int( m.layers.size() ) - 1;    // edge / angle outputs unused: their gradients are zero
    const double *rn0 = W + P.n_res[0], *rn1 = W + P.n_res[1], *rn2 = W + P.n_res[2];
    const double *re0 = W + P.e_res[0], *re1 = W + P.e_res[1], *ra0 = W + P.a_res;
    const double* dn = S.dn.data();
    const double* de = S.de.data();
    const double* da = S.da.data();
    dpa3_grow(S.dn_in, nd); dpa3_grow(S.de_in, size_t(n) * ed); dpa3_grow(S.da_in, size_t(na) * na * ad);
    dpa3_grow(S.dnj, size_t(n) * nd);
    double* dn_in = S.dn_in.data();
    double* de_in = S.de_in.data();
    double* da_in = S.da_in.data();
    double* dnj = S.dnj.data();
    DPA3Timer T;
    // residual identities
    for(int c = 0; c < nd; c++) dn_in[c] = dn[c];
    for(size_t x = 0; x < size_t(n) * ed; x++) de_in[x] = de[x];
    for(size_t x = 0; x < size_t(na) * na * ad; x++) da_in[x] = da[x];
    std::fill( dnj, dnj + size_t(n) * nd, 0.0 );
    std::vector<double>& dz = S.dz; dpa3_grow(dz, std::max( { nd, ed, ad, ax * ( ed + nd ) } ) );
    // node self
    for(int c = 0; c < nd; c++) dz[c] = rn0[c] * dn[c] * w.d_self[c];
    dpa3_lin_rows_t( m, P.node_self, 0, nd, dz.data(), dn_in );
    // symmetrization
    for(int c = 0; c < nd; c++) dz[c] = rn1[c] * dn[c] * w.d_sym[c];
    std::vector<double>& dsym = S.dsym; dsym.assign( size_t(ax) * ( ed + nd ), 0.0 );
    dpa3_lin_rows_t( m, P.node_sym, 0, ax * ( ed + nd ), dz.data(), dsym.data() );
    // grrg backward: G[a][c] = sum_x M[x][a] M[x][c] / 3
    auto grrg_bwd = [&](const double* M, int width, const double* dG, std::vector<double>& dM)
    {
      dM.assign( size_t(3) * width, 0.0 );
      for(int a = 0; a < ax; a++) for(int c = 0; c < width; c++)
      {
        const double gv = dG[size_t(a) * width + c] / 3.0;
        if( gv == 0.0 ) continue;
        for(int x = 0; x < 3; x++) { dM[x*width + c] += gv * M[x*width + a]; dM[x*width + a] += gv * M[x*width + c]; }
      }
    };
    grrg_bwd( w.hg_e.data(), ed, dsym.data(), S.dhg_e );
    grrg_bwd( w.hg_n.data(), nd, dsym.data() + size_t(ax) * ed, S.dhg_n );
    for(int k = 0; k < n; k++)
    {
      const double* h = g.h2.data() + 3 * k;
      const double sw = g.w[k] * m.s_sym;
      const double* e = w.e.data() + size_t(k) * ed;
      const double* nj = Nl + size_t( g.j[k] ) * nd;
      double* dek = de_in + size_t(k) * ed;
      double* dnjk = dnj + size_t(k) * nd;
      double gsw = 0.0;
      for(int x = 0; x < 3; x++)
      {
        double gh = 0.0;
        const double* dhe = S.dhg_e.data() + x * ed;
        const double* dhn = S.dhg_n.data() + x * nd;
        for(int c = 0; c < ed; c++) { dek[c] += dhe[c] * h[x] * sw; gh += dhe[c] * e[c]; }
        for(int c = 0; c < nd; c++) { dnjk[c] += dhn[c] * h[x] * sw; gh += dhn[c] * nj[c]; }
        if( geom ) S.gg.gh2[3*k + x] += gh * sw;
        gsw += gh * h[x];
      }
      if( geom ) S.gg.gw[k] += gsw * m.s_sym;
    }
    // node-edge messages and edge self updates
    T.mark(8);
    // dZm, dZe per edge, then dX = dZm Wne[nd:]^T + dZe Wes[nd:]^T, X = [ n_j | e ]
    dpa3_grow(S.dZm, size_t(n) * nd); dpa3_grow(S.dZe, size_t(n) * ed); dpa3_grow(S.dX, size_t(n) * ( nd + ed ));
    std::vector<double> dbase_n( nd, 0.0 ), dbase_e( ed, 0.0 );
    for(int k = 0; k < n; k++)
    {
      const double* am = w.a_msg.data() + size_t(k) * nd;
      const double* dmk = w.d_msg.data() + size_t(k) * nd;
      const double* dek_ = w.d_es.data() + size_t(k) * ed;
      double* dzm = S.dZm.data() + size_t(k) * nd;
      double* dze = S.dZe.data() + size_t(k) * ed;
      double gsw = 0.0;
      for(int c = 0; c < nd; c++)
      {
        const double gm = rn2[c] * dn[c] * m.s_msg;
        gsw += gm * am[c];
        dzm[c] = gm * g.w[k] * dmk[c];
        dbase_n[c] += dzm[c];
      }
      if( geom ) S.gg.gw[k] += gsw;
      if( !last ) for(int c = 0; c < ed; c++) { dze[c] = re0[c] * de[size_t(k)*ed + c] * dek_[c]; dbase_e[c] += dze[c]; }
      else for(int c = 0; c < ed; c++) dze[c] = 0.0;
    }
    dpa3_gemm( false, true, n, nd + ed, nd, 1.0, S.dZm.data(), nd, W + P.node_edge.w + size_t(nd) * nd, nd, 0.0, S.dX.data(), nd + ed );
    dpa3_gemm( false, true, n, nd + ed, ed, 1.0, S.dZe.data(), ed, W + P.edge_self.w + size_t(nd) * ed, ed, 1.0, S.dX.data(), nd + ed );
    for(int k = 0; k < n; k++)
    {
      const double* dx = S.dX.data() + size_t(k) * ( nd + ed );
      double* dj = dnj + size_t(k) * nd;
      double* dek = de_in + size_t(k) * ed;
      for(int c = 0; c < nd; c++) dj[c] += dx[c];
      for(int c = 0; c < ed; c++) dek[c] += dx[nd + c];
    }
    dpa3_lin_rows_t( m, P.node_edge, 0, nd, dbase_n.data(), dn_in );
    dpa3_lin_rows_t( m, P.edge_self, 0, nd, dbase_e.data(), dn_in );
    T.mark(9);
    if( last ) { T.mark(12); return; }
    // edge-angle: e_ang = act(z2), z2 = full W2 + b2
    std::vector<double> dred( size_t(na) * ed, 0.0 );
    dpa3_grow(S.dZe, size_t(n) * ed);
    for(size_t x = 0; x < size_t(n) * ed; x++) S.dZe[x] = re1[x % ed] * de[x] * w.d_e2[x];
    dpa3_grow(S.dX, size_t(n) * ed);
    dpa3_gemm( false, true, n, ed, ed, 1.0, S.dZe.data(), ed, W + P.edge_angle2.w, ed, 0.0, S.dX.data(), ed );   // dfull
    for(int k = 0; k < n; k++)
    {
      const double* dfull = S.dX.data() + size_t(k) * ed;
      if( k < na ) for(int c = 0; c < ed; c++) dred[size_t(k)*ed + c] = dfull[c];
      else if( !m.smooth_edge ) for(int c = 0; c < ed; c++) de_in[size_t(k)*ed + c] += dfull[c];
    }
    // reduction red_j = s_ang sum_k aw_j aw_k act(z_ea[j][k]); angle self a_self = act(z_as)
    std::vector<double> dnc( nca, 0.0 ), dec( size_t(na) * eca, 0.0 );
    T.mark(10);
    // dZ of the two angle linears [na*na][out], then: da += dZ W[0:ad]^T, dnc += (sum dZ) W[ad:ad+nca]^T,
    // dec_b += (column sums) W[ik]^T, dec_a += (row sums) W[ij]^T
    const size_t np2 = size_t(na) * na;
    dpa3_grow(S.dZa1, np2 * ed); dpa3_grow(S.dZa2, np2 * ad);
    std::vector<double> rs1( size_t(na) * ed, 0.0 ), cs1( size_t(na) * ed, 0.0 ), rs2( size_t(na) * ad, 0.0 ), cs2( size_t(na) * ad, 0.0 );
    std::vector<double> t1( ed, 0.0 ), t2( ad, 0.0 );
    for(int a = 0; a < na; a++) for(int b = 0; b < na; b++)
    {
      const size_t p = size_t(a) * na + b;
      const double f = g.aw[a] * g.aw[b] * m.s_ang;
      const double* aa = w.a_ea.data() + p * ed;
      const double* da1 = w.d_ea.data() + p * ed;
      const double* dr = dred.data() + size_t(a) * ed;
      double* d1 = S.dZa1.data() + p * ed;
      double* r1 = rs1.data() + size_t(a) * ed;
      double* c1 = cs1.data() + size_t(b) * ed;
      double gf = 0.0;
      for(int c = 0; c < ed; c++) { gf += dr[c] * aa[c]; d1[c] = dr[c] * f * da1[c]; r1[c] += d1[c]; c1[c] += d1[c]; }
      if( geom ) { S.gg.gaw[a] += gf * g.aw[b] * m.s_ang; S.gg.gaw[b] += gf * g.aw[a] * m.s_ang; }
      const double* das = w.d_as.data() + p * ad;
      const double* dd = da + p * ad;
      double* d2 = S.dZa2.data() + p * ad;
      double* r2 = rs2.data() + size_t(a) * ad;
      double* c2 = cs2.data() + size_t(b) * ad;
      for(int c = 0; c < ad; c++) { d2[c] = ra0[c] * dd[c] * das[c]; r2[c] += d2[c]; c2[c] += d2[c]; }
    }
    for(int a = 0; a < na; a++) { for(int c = 0; c < ed; c++) t1[c] += rs1[size_t(a)*ed + c]; for(int c = 0; c < ad; c++) t2[c] += rs2[size_t(a)*ad + c]; }
    if( na > 0 )
    {
      dpa3_gemm( false, true, int(np2), ad, ed, 1.0, S.dZa1.data(), ed, W + P.edge_angle1.w, ed, 1.0, da_in, ad );
      dpa3_gemm( false, true, int(np2), ad, ad, 1.0, S.dZa2.data(), ad, W + P.angle_self.w, ad, 1.0, da_in, ad );
      dpa3_lin_rows_t( m, P.edge_angle1, ad, nca, t1.data(), dnc.data() );
      dpa3_lin_rows_t( m, P.angle_self, ad, nca, t2.data(), dnc.data() );
      dpa3_gemm( false, true, na, eca, ed, 1.0, cs1.data(), ed, W + P.edge_angle1.w + size_t( ad + nca ) * ed, ed, 1.0, dec.data(), eca );
      dpa3_gemm( false, true, na, eca, ed, 1.0, rs1.data(), ed, W + P.edge_angle1.w + size_t( ad + nca + eca ) * ed, ed, 1.0, dec.data(), eca );
      dpa3_gemm( false, true, na, eca, ad, 1.0, cs2.data(), ad, W + P.angle_self.w + size_t( ad + nca ) * ad, ad, 1.0, dec.data(), eca );
      dpa3_gemm( false, true, na, eca, ad, 1.0, rs2.data(), ad, W + P.angle_self.w + size_t( ad + nca + eca ) * ad, ad, 1.0, dec.data(), eca );
    }
    T.mark(11);
    // compressed features
    if( m.compress == 2 )
    {
      dpa3_lin_rows_t( m, P.compress_n, 0, nd, dnc.data(), dn_in );
      for(int k = 0; k < na; k++) dpa3_lin_rows_t( m, P.compress_e, 0, ed, dec.data() + size_t(k) * eca, de_in + size_t(k) * ed );
    }
    else
    {
      for(int c = 0; c < nca; c++) dn_in[c] += dnc[c];
      for(int k = 0; k < na; k++) for(int c = 0; c < eca; c++) de_in[size_t(k)*ed + c] += dec[size_t(k)*eca + c];
    }
    T.mark(12);
  }

  // VJP of the initial edge and angle embeddings (de, da of layer 0 inputs) into the geometry gradients
  inline void dpa3_init_vjp(const DPA3Model& m, const double* de, const double* da, DPA3AtomScratch& S)
  {
    const DPA3Geom& g = S.g;
    const int ed = m.e_dim, ad = m.a_dim, n = g.n, na = g.na;
    const double* W = m.W.data();
    for(int k = 0; k < n; k++)
    {
      const double x = m.use_dist ? g.r[k] : g.s[k];
      double gx = 0.0;
      for(int c = 0; c < ed; c++)
      {
        double gz = de[size_t(k)*ed + c];
        if( !m.use_dist )
        {
          double z = W[m.edge_embd.w + c] * x + ( m.edge_embd.b != size_t(-1) ? W[m.edge_embd.b + c] : 0.0 );
          gz *= dpa3_dact( m, z );
        }
        gx += gz * W[m.edge_embd.w + c];
      }
      if( m.use_dist ) S.gg.gr[k] += gx; else S.gg.gs[k] += gx;
    }
    const double isp = 1.0 / std::sqrt( M_PI );
    for(int a = 0; a < na; a++) for(int b = 0; b < na; b++)
    {
      const double* dd = da + ( size_t(a) * na + b ) * ad;
      double gx = 0.0;
      for(int c = 0; c < ad; c++) gx += dd[c] * W[m.angle_embd.w + c];
      const double gc = gx * isp * ( 1.0 - 1e-6 );   // dE/d(u_a . u_b)
      for(int x = 0; x < 3; x++) { S.gg.gu[3*a + x] += gc * g.u[3*b + x]; S.gg.gu[3*b + x] += gc * g.u[3*a + x]; }
    }
  }

  // geometry gradients of the atom -> dE/dd per bond (added into gb [nbonds][3])
  inline void dpa3_geom_to_bonds(const DPA3Model& m, const DPA3AtomScratch& S, double* gb)
  {
    const DPA3Geom& g = S.g;
    const double* sd = m.dstd.data() + 4 * g.ti;
    for(int k = 0; k < g.n; k++)
    {
      const double* d = g.d.data() + 3 * k;
      const double r = g.r[k], q = g.q[k], w = g.w[k], dw = g.dw[k];
      double v[3] = { 0.0, 0.0, 0.0 };
      // h2_x = w d_x / q^2 / sd_x (+ const)
      double radial = 0.0;
      for(int x = 0; x < 3; x++)
      {
        const double gh = S.gg.gh2[3*k + x] / sd[1 + x];
        v[x] += gh * w / ( q * q );
        radial += gh * d[x] * ( dw / ( q * q ) - 2.0 * w / ( q * q * q ) );
      }
      radial += S.gg.gw[k] * dw;                                     // switch
      radial += S.gg.gr[k];                                          // distance
      radial += S.gg.gs[k] / sd[0] * ( dw / q - w / ( q * q ) );     // normalized radial row
      if( k < g.na )
      {
        double aw, adw;
        dpa3_switch( m, r, m.ars, m.arc, aw, adw );
        radial += S.gg.gaw[k] * adw;
        // u = d / (r + eps): du/dd = I/(r+eps) - d d^T / ( r (r+eps)^2 )
        const double re = r + 1e-6;
        const double* gu = S.gg.gu.data() + 3 * k;
        const double gd = gu[0]*d[0] + gu[1]*d[1] + gu[2]*d[2];
        for(int x = 0; x < 3; x++) v[x] += gu[x] / re;
        radial -= gd / ( re * re );
      }
      for(int x = 0; x < 3; x++) gb[ 3 * g.bond[k] + x ] += v[x] + radial * d[x] / r;
    }
  }

}
