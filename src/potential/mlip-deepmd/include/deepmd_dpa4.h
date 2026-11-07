/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

// DPA4 (SeZM) descriptor, host only: SO(3) node features x[D][C] (D = (lmax+1)^2), l <= 3, m <= 1, one or several focus streams, one
// attention head, as the presets of deepmd (exaStamp/scripts/python/deepmd_dpa4.py is the numpy specification).
//
//   per edge j -> i : quaternion q(d_ij), Wigner blocks D(q), envelope, Bessel or Gaussian radial basis, radial features
//   x0 = FiLM(type embedding; env seed) on l = 0, geometric initial embedding on l > 0
//   block: xp = pre_mix(optional norm(x));  edge message: rotate xp_j to the edge frame, radial mixer, SO(2) layers, rotate back
//          out = gate_i sum_j alpha_ij msg_ij;  out += grid(out, xp_i);  x1 = x + optional norm(post_mix(out));
//          x2 = x1 + ffn(norm(x1))
//   descriptor = x[l = 0] + read-out(x), GLU fitting network
//
// Phases (dpa4_compute): x of every block for all atoms (ghosts updated at the input of each block). Reverse: each atom
// recomputes its block and goes through its vector-Jacobian product, adding dE/dxp to its source atoms (row locks); ghost
// rows are folded back to their owners. Gradients on D, radial features and envelope accumulate per edge over the blocks;
// the geometric chain to dE/dd_ij, forces and virials run at the end.

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/particle_neighbors/chunk_neighbors_apply.h>
#include <exanb/mpi/flat_buffer_ghosts.h>
#include <exaStamp/unit_system.h>

#include <algorithm>
#include <limits>

#include "deepmd_dpa3.h"
#include "deepmd_dpa3_driver.h"

namespace exaStamp
{
  using namespace exanb;

  // stage profile (DPA4_PROF=1): wall time per evaluation, printed at exit
  struct DPA4Prof
  {
    bool on = std::getenv("DPA4_PROF") != nullptr;
    int evals = 0;
    double t[6] = {};
    ~DPA4Prof()
    {
      if( !on || evals == 0 ) return;
      static const char* names[6] = { "neighbors, geometry, init", "forward blocks", "read-out, fitting", "reverse blocks", "geometric chain", "forces" };
      std::fprintf( stderr, "[DPA4] profile, ms per evaluation (%d evaluations):\n", evals );
      double s = 0.0;
      for(int k = 0; k < 6; k++) { std::fprintf( stderr, "  %-26s %10.3f\n", names[k], 1e3 * t[k] / evals ); s += t[k]; }
      std::fprintf( stderr, "  %-26s %10.3f\n", "total", 1e3 * s / evals );
    }
  };
  inline DPA4Prof& dpa4_prof() { static DPA4Prof p; return p; }

  struct DPA4Rms { std::vector<double> scale; double eps = 0.0; };

  struct DPA4EqRms { std::vector<double> scale, bias, bw; double eps = 0.0; };   // scale [L][C], bias [C], bw [D]

  // weights of one focus stream
  struct DPA4Focus
  {
    std::vector< std::vector<double> > lin0, linm, gate;    // [(L C)][(L C)], [lmax C][2 lmax C], [C][lmax C] (no gate in the last layer)
    DPA4Rms qk_norm, gate_norm, comp_norm;
    std::vector<double> q, k, logit_w, gate_w, comp_w;      // [C][C], [C][C], [C], [C], [C]
    double log_null = 0.0;
    std::vector<double> ggate, gres;                        // [2C][C], [C]
  };

  struct DPA4Block
  {
    std::vector<double> pre_mix, post_mix;                  // [L][C][W], [L][W][C], W = n_focus C
    std::vector<DPA4Focus> foc;
    bool compete = false;
    std::vector<double> rproj;                              // [C][W], several focus streams only
    std::vector<double> mixw, mixb;                         // [L W][nk], [W]; empty: no mixer
    std::vector<double> gexp, gcon;                         // [L][C][K C], [L][K C][C], shared by the focus streams
    DPA4EqRms pre_norm, post_norm, ffn_norm;                // pre and post norms: empty scale, identity
    std::vector<double> flin1, flin2, fl, fr, fo, fgate;    // [L][C][2 Ce], [L][Ce][C], [Ch][Ch] x3, [2 Ch][Ch]
  };

  struct DPA4Model
  {
    int F_n = 1, W = 0;                                             // focus streams, hidden width F_n C
    double focus_tau = 1.0, focus_ls = 0.0;                         // competition temperature, label smoothing
    int nt = 0, C = 0, lmax = 0, L = 0, D = 0, Dm = 0, nnei = 0, K = 0, Ch = 0, Ce = 0, Ng = 0, nrad = 0, nk = 0;
    int emb = 0, axis = 0, tdim = 0, hrf = 0, fit_h = 0;
    bool store_forward = false;                                     // keep the forward values for the reverse pass
    double gauss = 0.0;                                             // gaussian radial basis: exponent coefficient, centres in freq; 0: bessel
    double rc = 0.0, eps = 0.0, weps = 0.0, deg_floor = 0.0, film_scale = 0.0, film_shift = 0.0, seed_eps = 0.0;
    std::vector<double> env_c, rbf_c, freq, te, tt;                 // series coefficients, frequencies, type embeddings
    std::vector<double> rf1, rf2; DPA4Rms rf_norm;
    std::vector<double> s_r1, s_r2, s_g1, s_g2, s_out;              // env seed layers, s_* activation flags below
    bool s_a[5] = { true, true, true, true, false };
    int s_r = 0, s_g = 0, s_gh = 0;
    DPA4Rms film_s, film_f;
    std::vector<double> C_l2;                                       // [25][35] degree 4 monomials
    std::vector<int> mono;                                          // [35][4]
    std::vector<double> C_l3; std::vector<int> exp_l3;              // lmax 3: [49][84] degree 6 monomials, [84][4] exponents
    int nDq = 35;                                                   // packed Wigner blocks, sum (2 l + 1)^2
    std::vector<double> T, F, rescale, sweight;                     // to grid [Ng][D K], from grid [D K][Ng], [D], [D][K]
    std::vector<int> coeff, deg_m, deg, sign_g, mix_r, mix_s, mix_c;// rows of D per reduced row, their degree, degree per D row, gate group, kernel (r, s, compact)
    std::vector<int> z_row, z_col, z_slot;                          // geometric initial embedding: D[col, row], radial slot
    std::vector<int> rot_r, rot_d, rot_o;                           // reduced row, D column, offset in the packed D
    std::vector<DPA4Block> blocks;
    // read-out
    std::vector<double> rlin1, rl, rr, ro, rgate, rlin2;            // [L][C][2Ce], [2Ch][Hid] x2, [Hid][Ch], [2Ch][Ch], [Ch][C] (l = 0, frame zero)
    int rhid = 0;
    // fitting
    std::vector<double> fw, fb, fout, ebias;                        // [C][2 h], [2 h], [h]
    double fout_b = 0.0;
  };

  // forward values read by the reverse pass, of one focus stream
  struct DPA4FocusState
  {
    std::vector< std::vector<double> > Y, GG, Sd;                   // SO(2) layers but the last one
    std::vector<double> XLf, alpha, q0, sum, gl, gr, cin, zgate, pair, xp;
    double g = 0.0;
  };

  // of one atom in one block
  struct DPA4AtomState
  {
    std::vector<DPA4FocusState> fs;
    std::vector<double> X0, rad, kern, af, alf, XW0, mixed, wn1, x1, wn2, fgl, fgr, fco, fzgate, fpair;
  };

  struct DPA4Buffers
  {
    DPA3Buffers c;
    std::vector<int> ecount;                                        // kept neighbors per atom row
    std::vector<double> envv, rr, rbf, Dq, rf, GD, Grf, Genv, Grbf, Gv, isd, deg, gie, seed;
    std::vector<double> Xp, Kn, GXp, GKn, dX1;
    std::vector< std::vector<double> > X, G;
    std::vector< std::vector<DPA4AtomState> > store;                // [block][atom], dpa4_store_forward only
  };

  // ------------------------------------------------------------------ small dense helpers

  inline double dpa4_sig(double x) { return 1.0 / ( 1.0 + std::exp(-x) ); }

  // array forms, one exponential per value, batched with Accelerate vvexp on macOS
  inline void dpa4_exp_n(double* y, int n)
  {
#   if defined(XSTAMP_DEEPMD_BLAS) && defined(__APPLE__)
    vvexp( y, y, &n );
#   else
    for(int k = 0; k < n; k++) y[k] = std::exp( y[k] );
#   endif
  }
  // s[k] = sigmoid( x[k] ), s may be x
  inline void dpa4_sig_n(const double* x, int n, double* s)
  {
    static thread_local std::vector<double> e;
    if( int(e.size()) < n ) e.resize(n);
    for(int k = 0; k < n; k++) e[k] = -x[k];
    dpa4_exp_n( e.data(), n );
    for(int k = 0; k < n; k++) s[k] = 1.0 / ( 1.0 + e[k] );
  }
  // x[k] = silu( x[k] )
  inline void dpa4_silu_n(double* x, int n)
  {
    static thread_local std::vector<double> e;
    if( int(e.size()) < n ) e.resize(n);
    for(int k = 0; k < n; k++) e[k] = -x[k];
    dpa4_exp_n( e.data(), n );
    for(int k = 0; k < n; k++) x[k] = x[k] / ( 1.0 + e[k] );
  }
  // d[k] *= silu'( z[k] )
  inline void dpa4_dsilu_n(const double* z, int n, double* d)
  {
    static thread_local std::vector<double> e;
    if( int(e.size()) < n ) e.resize(n);
    for(int k = 0; k < n; k++) e[k] = -z[k];
    dpa4_exp_n( e.data(), n );
    for(int k = 0; k < n; k++) { const double s = 1.0 / ( 1.0 + e[k] ); d[k] *= s + z[k] * s * ( 1.0 - s ); }
  }

  // x[D][Cin] -> y[D][Cout], weights per degree [L][Cin][Cout]
  inline void dpa4_so3(int lmax, const double* W, int Cin, int Cout, const double* x, double* y, double beta = 0.0)
  {
    for(int l = 0; l <= lmax; l++)
    {
      const int b = l * l, n = 2 * l + 1;
      dpa3_gemm( false, false, n, Cout, Cin, 1.0, x + size_t(b) * Cin, Cin, W + size_t(l) * Cin * Cout, Cout, beta, y + size_t(b) * Cout, Cout );
    }
  }
  // transpose: dx[D][Cin] (+)= dy[D][Cout] W^T
  inline void dpa4_so3_t(int lmax, const double* W, int Cin, int Cout, const double* dy, double* dx, double beta = 0.0)
  {
    for(int l = 0; l <= lmax; l++)
    {
      const int b = l * l, n = 2 * l + 1;
      dpa3_gemm( false, true, n, Cin, Cout, 1.0, dy + size_t(b) * Cout, Cout, W + size_t(l) * Cin * Cout, Cout, beta, dx + size_t(b) * Cin, Cin );
    }
  }
  // y[C] = x[C] / sqrt( mean(x^2) + eps ) scale[C]
  inline void dpa4_rms(const DPA4Rms& N, int C, const double* x, double* y)
  {
    if( N.scale.empty() ) { std::copy( x, x + C, y ); return; }   // edge_norm off: identity
    double s = 0.0;
#   pragma omp simd reduction(+:s)
    for(int c = 0; c < C; c++) s += x[c] * x[c];
    const double inv = 1.0 / std::sqrt( s / C + N.eps );
    for(int c = 0; c < C; c++) y[c] = x[c] * inv * N.scale[c];
  }
  inline void dpa4_rms_t(const DPA4Rms& N, int C, const double* x, const double* dy, double* dx)
  {
    if( N.scale.empty() ) { for(int c = 0; c < C; c++) dx[c] += dy[c]; return; }
    double s = 0.0, t = 0.0;
#   pragma omp simd reduction(+:s,t)
    for(int c = 0; c < C; c++) { s += x[c] * x[c]; t += dy[c] * x[c] * N.scale[c]; }
    const double inv = 1.0 / std::sqrt( s / C + N.eps );
    const double f = t * inv * inv * inv / C;
    for(int c = 0; c < C; c++) dx[c] += dy[c] * inv * N.scale[c] - f * x[c];
  }
  // equivariant norm, x[D][C] -> y[D][C]; keeps inv and the centered row 0 in w (C + 1 values)
  inline void dpa4_eqrms(const DPA4EqRms& N, const DPA4Model& M, const double* x, double* y, double* w)
  {
    if( N.scale.empty() ) { std::copy( x, x + size_t(M.D) * M.C, y ); return; }
    const int C = M.C, D = M.D;
    double mu = 0.0;
#   pragma omp simd reduction(+:mu)
    for(int c = 0; c < C; c++) mu += x[c];
    mu /= C;
    double mv = 0.0;
    for(int c = 0; c < C; c++) { w[c] = x[c] - mu; mv += w[c] * w[c] * N.bw[0]; }
    for(int d = 1; d < D; d++)
    {
      const double* xd = x + size_t(d) * C;
      double s = 0.0;
#     pragma omp simd reduction(+:s)
      for(int c = 0; c < C; c++) s += xd[c] * xd[c];
      mv += s * N.bw[d];
    }
    const double inv = 1.0 / std::sqrt( mv + N.eps );
    w[C] = inv;
    for(int c = 0; c < C; c++) y[c] = w[c] * inv * N.scale[c] + N.bias[c];
    for(int d = 1; d < D; d++)
    {
      const double* xd = x + size_t(d) * C;
      const double* sc = N.scale.data() + size_t( M.deg[d] ) * C;
      double* yd = y + size_t(d) * C;
      for(int c = 0; c < C; c++) yd[c] = xd[c] * inv * sc[c];
    }
  }
  inline void dpa4_eqrms_t(const DPA4EqRms& N, const DPA4Model& M, const double* x, const double* w, const double* dy, double* dx)
  {
    if( N.scale.empty() ) { for(size_t t = 0; t < size_t(M.D) * M.C; t++) dx[t] += dy[t]; return; }
    const int C = M.C, D = M.D;
    const double inv = w[C];
    double S = 0.0;
    for(int c = 0; c < C; c++) S += dy[c] * N.scale[c] * w[c];
    for(int d = 1; d < D; d++)
    {
      const double* xd = x + size_t(d) * C;
      const double* dd = dy + size_t(d) * C;
      const double* sc = N.scale.data() + size_t( M.deg[d] ) * C;
      double s = 0.0;
#     pragma omp simd reduction(+:s)
      for(int c = 0; c < C; c++) s += dd[c] * sc[c] * xd[c];
      S += s;
    }
    const double f = S * inv * inv * inv;
    double mean = 0.0;
    for(int c = 0; c < C; c++) mean += dy[c] * N.scale[c] * inv - f * N.bw[0] * w[c];
    mean /= C;
    for(int c = 0; c < C; c++) dx[c] += dy[c] * N.scale[c] * inv - f * N.bw[0] * w[c] - mean;
    for(int d = 1; d < D; d++)
    {
      const double* xd = x + size_t(d) * C;
      const double* dd = dy + size_t(d) * C;
      const double* sc = N.scale.data() + size_t( M.deg[d] ) * C;
      double* dxd = dx + size_t(d) * C;
      const double g = f * N.bw[d];
      for(int c = 0; c < C; c++) dxd[c] += dd[c] * sc[c] * inv - g * xd[c];
    }
  }

  // C2 envelope series: u^4 sum_k c_k x^k, x = r / rc, u = 1 - x; value and radial derivative
  inline void dpa4_envelope(const std::vector<double>& c, double r, double rc, double& v, double& dv)
  {
    if( c.empty() ) { v = 1.0; dv = 0.0; return; }
    const double u0 = ( rc - r ) / rc, u = u0 < 0.0 ? 0.0 : ( u0 > 1.0 ? 1.0 : u0 ), x = 1.0 - u;
    double s = c.back(), ds = 0.0;
    for(int k = int(c.size()) - 2; k >= 0; k--) { ds = s + x * ds; s = c[k] + x * s; }
    const double u3 = u * u * u;
    v = u3 * u * s;
    dv = ( u0 < 0.0 || u0 > 1.0 ) ? 0.0 : ( -4.0 * u3 * s + u3 * u * ds ) / rc;
  }

  // quaternion of the edge and its Jacobian J[4][3] dq / dv
  inline void dpa4_quaternion(const double* v, double r, double eps, double weps, double* q, double* J)
  {
    const double rr = std::sqrt( r * r + eps * eps );
    const double u[3] = { v[0] / rr, v[1] / rr, v[2] / rr };
    double du[3][3];
    for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) du[a][b] = ( a == b ? 1.0 / rr : 0.0 ) - v[a] * v[b] / ( rr * rr * rr );
    // normalized charts: P / sqrt( |P|^2 + eps^2 ), P linear in u
    auto chart = [&](const double* P, const int (*A)[3], double* out, double (*dout)[3])
    {
      const double s = P[0]*P[0] + P[1]*P[1] + P[2]*P[2] + P[3]*P[3] + eps*eps;
      const double nn = std::sqrt( s );
      for(int a = 0; a < 4; a++) out[a] = P[a] / nn;
      double dP[4][3];
      for(int a = 0; a < 4; a++) for(int b = 0; b < 3; b++)
      {
        double t = 0.0;
        for(int k = 0; k < 3; k++) t += A[a][k] * du[k][b];
        dP[a][b] = t;
      }
      for(int b = 0; b < 3; b++)
      {
        const double pd = P[0]*dP[0][b] + P[1]*dP[1][b] + P[2]*dP[2][b] + P[3]*dP[3][b];
        for(int a = 0; a < 4; a++) dout[a][b] = dP[a][b] / nn - P[a] * pd / ( nn * s );
      }
    };
    static const int Ap[4][3] = { {0,0,1}, {0,1,0}, {-1,0,0}, {0,0,0} };
    static const int An[4][3] = { {-1,0,0}, {0,0,0}, {0,0,-1}, {0,1,0} };
    const double Pp[4] = { 1.0 + u[2], u[1], -u[0], 0.0 };
    const double Pn[4] = { -u[0], 0.0, 1.0 - u[2], u[1] };
    double qp[4], qn[4], dqp[4][3], dqn[4][3];
    chart( Pp, Ap, qp, dqp );
    chart( Pn, An, qn, dqn );
    const double t = 0.5 * ( u[2] + 1.0 );
    double w = 0.0, dw = 0.0;
    if( t >= 1.0 ) w = 1.0;
    else if( t > 0.0 )
    {
      const double fe = std::numeric_limits<double>::epsilon();
      const double a = std::max( t, fe ), bb = std::max( 1.0 - t, fe );
      const double l = std::exp( -1.0 / a ), rg = std::exp( -1.0 / bb );
      w = l / ( l + rg );
      dw = l * rg * ( 1.0 / ( a * a ) + 1.0 / ( bb * bb ) ) / ( ( l + rg ) * ( l + rg ) );
    }
    const double sg = ( qn[0]*qp[0] + qn[1]*qp[1] + qn[2]*qp[2] + qn[3]*qp[3] ) < 0.0 ? -1.0 : 1.0;
    double B[4], dB[4][3];
    for(int a = 0; a < 4; a++)
    {
      B[a] = ( 1.0 - w ) * qn[a] + w * sg * qp[a];
      for(int b = 0; b < 3; b++) dB[a][b] = ( 1.0 - w ) * dqn[a][b] + w * sg * dqp[a][b] + dw * 0.5 * du[2][b] * ( sg * qp[a] - qn[a] );
    }
    auto normalize = [](const double* P, const double (*dP)[3], double e, double* out, double (*dout)[3])
    {
      const double s = P[0]*P[0] + P[1]*P[1] + P[2]*P[2] + P[3]*P[3] + e*e;
      const double nn = std::sqrt( s );
      for(int a = 0; a < 4; a++) out[a] = P[a] / nn;
      for(int b = 0; b < 3; b++)
      {
        const double pd = P[0]*dP[0][b] + P[1]*dP[1][b] + P[2]*dP[2][b] + P[3]*dP[3][b];
        for(int a = 0; a < 4; a++) dout[a][b] = dP[a][b] / nn - P[a] * pd / ( nn * s );
      }
    };
    double q1[4], dq1[4][3], dq2[4][3];
    normalize( B, dB, eps, q1, dq1 );
    normalize( q1, dq1, weps, q, dq2 );
    for(int a = 0; a < 4; a++) for(int b = 0; b < 3; b++) J[a*3 + b] = dq2[a][b];
  }

  // packed Wigner blocks of one edge: D[0] = 1, then l = 1 (9), l = 2 (25), l = 3 (49)
  inline void dpa4_wigner(const DPA4Model& M, const double* q, double* Dq)
  {
    const double w = q[0], x = q[1], y = q[2], z = q[3];
    Dq[0] = 1.0;
    if( M.lmax < 1 ) return;
    const double R[3][3] = { { 1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y) },
                             { 2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x) },
                             { 2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y) } };
    static const int p[3] = { 1, 2, 0 };
    static const double s[3] = { -1.0, -1.0, 1.0 };
    for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) Dq[1 + a*3 + b] = R[p[a]][p[b]] * s[a] * s[b];
    if( M.lmax < 2 ) return;
    double mono[35];
    for(int m = 0; m < 35; m++) { const int* i = M.mono.data() + m * 4; mono[m] = q[i[0]] * q[i[1]] * q[i[2]] * q[i[3]]; }
    for(int a = 0; a < 25; a++)
    {
      const double* c = M.C_l2.data() + a * 35;
      double t = 0.0;
      for(int m = 0; m < 35; m++) t += c[m] * mono[m];
      Dq[10 + a] = t;
    }
    if( M.lmax < 3 ) return;
    // l = 3: degree 6 monomials, exponents of (w, x, y, z)
    double pw[4][7];
    for(int a = 0; a < 4; a++) { pw[a][0] = 1.0; for(int k = 1; k < 7; k++) pw[a][k] = pw[a][k-1] * q[a]; }
    double mono6[84];
    for(int m = 0; m < 84; m++) { const int* i = M.exp_l3.data() + m * 4; mono6[m] = pw[0][i[0]] * pw[1][i[1]] * pw[2][i[2]] * pw[3][i[3]]; }
    for(int a = 0; a < 49; a++)
    {
      const double* c = M.C_l3.data() + a * 84;
      double t = 0.0;
      for(int m = 0; m < 84; m++) t += c[m] * mono6[m];
      Dq[35 + a] = t;
    }
  }
  // dE/dq from dE/dD (packed)
  inline void dpa4_wigner_t(const DPA4Model& M, const double* q, const double* G, double* gq)
  {
    gq[0] = gq[1] = gq[2] = gq[3] = 0.0;
    if( M.lmax < 1 ) return;
    const double w = q[0], x = q[1], y = q[2], z = q[3];
    static const int p[3] = { 1, 2, 0 };
    static const double s[3] = { -1.0, -1.0, 1.0 };
    double g[3][3];
    for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) g[p[a]][p[b]] = G[1 + a*3 + b] * s[a] * s[b];
    gq[0] += 2.0 * ( - z * g[0][1] + y * g[0][2] + z * g[1][0] - x * g[1][2] - y * g[2][0] + x * g[2][1] );
    gq[1] += -4.0 * x * ( g[1][1] + g[2][2] ) + 2.0 * ( y * g[0][1] + z * g[0][2] + y * g[1][0] - w * g[1][2] + z * g[2][0] + w * g[2][1] );
    gq[2] += -4.0 * y * ( g[0][0] + g[2][2] ) + 2.0 * ( x * g[0][1] + w * g[0][2] + x * g[1][0] + z * g[1][2] - w * g[2][0] + z * g[2][1] );
    gq[3] += -4.0 * z * ( g[0][0] + g[1][1] ) + 2.0 * ( - w * g[0][1] + x * g[0][2] + w * g[1][0] + y * g[1][2] + x * g[2][0] + y * g[2][1] );
    if( M.lmax < 2 ) return;
    double cm[35];
    for(int m = 0; m < 35; m++) cm[m] = 0.0;
    for(int a = 0; a < 25; a++)
    {
      const double ga = G[10 + a];
      if( ga == 0.0 ) continue;
      const double* c = M.C_l2.data() + a * 35;
      for(int m = 0; m < 35; m++) cm[m] += ga * c[m];
    }
    for(int m = 0; m < 35; m++)
    {
      const int* i = M.mono.data() + m * 4;
      const double a = q[i[0]], b = q[i[1]], c = q[i[2]], d = q[i[3]];
      gq[i[0]] += cm[m] * b * c * d;
      gq[i[1]] += cm[m] * a * c * d;
      gq[i[2]] += cm[m] * a * b * d;
      gq[i[3]] += cm[m] * a * b * c;
    }
    if( M.lmax < 3 ) return;
    double pw[4][7];
    for(int a = 0; a < 4; a++) { pw[a][0] = 1.0; for(int k = 1; k < 7; k++) pw[a][k] = pw[a][k-1] * q[a]; }
    double cm6[84];
    for(int m = 0; m < 84; m++) cm6[m] = 0.0;
    for(int a = 0; a < 49; a++)
    {
      const double ga = G[35 + a];
      if( ga == 0.0 ) continue;
      const double* c = M.C_l3.data() + a * 84;
      for(int m = 0; m < 84; m++) cm6[m] += ga * c[m];
    }
    for(int m = 0; m < 84; m++)
    {
      const int* i = M.exp_l3.data() + m * 4;
      for(int a = 0; a < 4; a++)
      {
        if( i[a] == 0 ) continue;
        double t = cm6[m] * i[a] * pw[a][ i[a] - 1 ];
        for(int b = 0; b < 4; b++) if( b != a ) t *= pw[b][ i[b] ];
        gq[a] += t;
      }
    }
  }

  // ------------------------------------------------------------------ scalar path of the grid nets

  // v = c gate(pair), v[0] += SwiGLU(pair); c, v [rows][Cf], pair [2 Cf]
  inline void dpa4_scalar(const double* Gw, int Cf, int rows, const double* c, const double* pair, double* v, double* zgate)
  {
    dpa3_gemm( false, false, 1, Cf, 2 * Cf, 1.0, pair, 2 * Cf, Gw, Cf, 0.0, zgate, Cf );
    for(int k = 0; k < Cf; k++) zgate[k] = dpa4_sig( zgate[k] );
    for(int r = 0; r < rows; r++)
    {
      const double* cr = c + size_t(r) * Cf;
      double* vr = v + size_t(r) * Cf;
      for(int k = 0; k < Cf; k++) vr[k] = cr[k] * zgate[k];
    }
    for(int k = 0; k < Cf; k++) v[k] += pair[k] * dpa4_sig( pair[k] ) * pair[Cf + k];
  }
  // from dv: dc[rows][Cf], dpair[2 Cf] (assigned)
  inline void dpa4_scalar_t(const double* Gw, int Cf, int rows, const double* c, const double* pair, const double* zgate, const double* dv, double* dc, double* dpair, double* tmp)
  {
    for(int k = 0; k < Cf; k++) tmp[k] = 0.0;
    for(int r = 0; r < rows; r++)
    {
      const double* cr = c + size_t(r) * Cf;
      const double* dr = dv + size_t(r) * Cf;
      double* er = dc + size_t(r) * Cf;
      for(int k = 0; k < Cf; k++) { er[k] = dr[k] * zgate[k]; tmp[k] += dr[k] * cr[k]; }
    }
    for(int k = 0; k < Cf; k++) tmp[k] *= zgate[k] * ( 1.0 - zgate[k] );
    dpa3_gemm( false, true, 1, 2 * Cf, Cf, 1.0, tmp, Cf, Gw, Cf, 0.0, dpair, 2 * Cf );
    for(int k = 0; k < Cf; k++)
    {
      const double s = dpa4_sig( pair[k] );
      dpair[k] += dv[k] * ( s + pair[k] * s * ( 1.0 - s ) ) * pair[Cf + k];
      dpair[Cf + k] += dv[k] * pair[k] * s;
    }
  }

  // ------------------------------------------------------------------ one atom in one block

  // temporaries of one thread
  struct DPA4Scratch
  {
    DPA4AtomState own;                                              // state of the current atom when forward values are not stored
    std::vector< std::vector<double> > XL;
    std::vector<double> Ylast, rad, ex, logit, l0, h, vv, res, out;
    std::vector<double> wedge, y, fc, fv;
    std::vector<double> zb, rad0, XW, outw, doutw, dXW, dX0, Grw, dAf, dxf;
    std::vector<double> tmp, tmp2, tmp3, tmp4, dX, dY, dG, dz, dalpha, ddso, dout, dc, dp, dleft, dright, dh, dpairv, df, dxp, drad, dres;
  };

  // forward of block b for atom i; Xp, Kn: nodes (ghosts included); x, y: node row in / out
  inline void dpa4_block_forward(const DPA4Model& M, const DPA4Block& P, const DPA4Buffers& BB, int i, const int64_t e0, const int n,
                                 const int* nbh, const double* x, DPA4AtomState& A, DPA4Scratch& S, double* x2)
  {
    const int C = M.C, D = M.D, Dm = M.Dm, L = M.L, DK = D * M.K, K = M.K, Ch = M.Ch, Ce = M.Ce, F = M.F_n, W = M.W;
    const size_t DC = size_t(D) * C, DmC = size_t(Dm) * C, DW = size_t(D) * W, DmW = size_t(Dm) * W;
    const double* Xp = BB.Xp.data();
    const double* xp = Xp + size_t(i) * DW;
    const int nl = P.foc[0].lin0.size();
    if( int(A.fs.size()) < F ) A.fs.resize( F );
    std::vector<double>& Rd = M.nk > 0 ? S.rad : A.rad;                // the reverse pass reads it without a mixer only
    dpa3_grow( A.X0, size_t(n) * DmW ); dpa3_grow( S.XW, size_t(n) * DmW ); dpa3_grow( S.rad0, size_t(n) * L * C ); dpa3_grow( Rd, size_t(n) * L * W );
    // radial features of the edges, types fused; projected on the hidden width when there are several focus streams
    const int ty = BB.c.type[i];
    const double* tei = M.te.data() + size_t(ty) * C;
    for(int e = 0; e < n; e++)
    {
      const double* rf = BB.rf.data() + size_t(e0 + e) * L * C;
      const double* tej = M.te.data() + size_t( BB.c.type[ nbh[e0 + e] ] ) * C;
      double* r = S.rad0.data() + size_t(e) * L * C;
      for(int l = 0; l < L; l++) for(int c = 0; c < C; c++) r[l*C + c] = rf[l*C + c] + tei[c] + tej[c];
    }
    if( F > 1 ) dpa3_gemm( false, false, n * L, W, C, 1.0, S.rad0.data(), C, P.rproj.data(), W, 0.0, Rd.data(), W );
    else std::copy( S.rad0.begin(), S.rad0.begin() + size_t(n) * L * C, Rd.begin() );
    // rotate the source features to the edge frame
    std::fill( A.X0.begin(), A.X0.begin() + size_t(n) * DmW, 0.0 );
    const int nrot = M.rot_r.size();
    for(int e = 0; e < n; e++)
    {
      const double* xs = Xp + size_t( nbh[e0 + e] ) * DW;
      const double* Dq = BB.Dq.data() + size_t(e0 + e) * M.nDq;
      double* xl = A.X0.data() + size_t(e) * DmW;
      for(int t = 0; t < nrot; t++)
      {
        const double f = Dq[ M.rot_o[t] ];
        const double* sc = xs + size_t( M.rot_d[t] ) * W;
        double* dc = xl + size_t( M.rot_r[t] ) * W;
#       pragma omp simd
        for(int c = 0; c < W; c++) dc[c] += f * sc[c];
      }
    }
    // radial mixer
    double* X = S.XW.data();
    if( M.nk > 0 )
    {
      dpa3_grow( A.kern, size_t(n) * M.nk );
      dpa3_gemm( false, false, n, M.nk, L * W, 1.0, Rd.data(), L * W, P.mixw.data(), M.nk, 0.0, A.kern.data(), M.nk );
      std::fill( X, X + size_t(n) * DmW, 0.0 );
      const int nm = M.mix_r.size();
      for(int e = 0; e < n; e++)
      {
        const double* k = A.kern.data() + size_t(e) * M.nk;
        const double* x0 = A.X0.data() + size_t(e) * DmW;
        double* xo = X + size_t(e) * DmW;
        for(int t = 0; t < nm; t++)
        {
          const double f = k[ M.mix_c[t] ];
          const double* sc = x0 + size_t( M.mix_s[t] ) * W;
          double* dc = xo + size_t( M.mix_r[t] ) * W;
#         pragma omp simd
          for(int c = 0; c < W; c++) dc[c] += f * sc[c] * P.mixb[c];
        }
      }
    }
    else
    {
      for(int e = 0; e < n; e++)
      {
        const double* x0 = A.X0.data() + size_t(e) * DmW;
        const double* rd = Rd.data() + size_t(e) * L * W;
        double* xo = X + size_t(e) * DmW;
        for(int r = 0; r < Dm; r++) { const double* rr = rd + size_t( M.deg_m[r] ) * W; for(int c = 0; c < W; c++) xo[size_t(r)*W + c] = x0[size_t(r)*W + c] * rr[c]; }
      }
    }
    // cross-focus competition on the scalar row after the mixer: af = softmax, alpha = smoothed weights
    dpa3_grow( A.af, size_t(n) * F ); dpa3_grow( A.alf, size_t(n) * F ); dpa3_grow( S.tmp, C );
    if( F > 1 && P.compete )
    {
      dpa3_grow( A.XW0, size_t(n) * W );
      for(int e = 0; e < n; e++)
      {
        const double* x0 = X + size_t(e) * DmW;
        std::copy( x0, x0 + W, A.XW0.begin() + size_t(e) * W );
        double* a = A.af.data() + size_t(e) * F;
        double* al = A.alf.data() + size_t(e) * F;
        double mx = -std::numeric_limits<double>::infinity();
        for(int f = 0; f < F; f++)
        {
          dpa4_rms( P.foc[f].comp_norm, C, x0 + size_t(f) * C, S.tmp.data() );
          double t = 0.0;
          for(int c = 0; c < C; c++) t += S.tmp[c] * P.foc[f].comp_w[c];
          a[f] = t / M.focus_tau;
          mx = std::max( mx, a[f] );
        }
        double den = 0.0;
        for(int f = 0; f < F; f++) { a[f] = std::exp( a[f] - mx ); den += a[f]; }
        for(int f = 0; f < F; f++) { a[f] /= den; al[f] = a[f] * ( 1.0 - M.focus_ls ) + M.focus_ls / F; }
      }
    }
    else std::fill( A.alf.begin(), A.alf.begin() + size_t(n) * F, 1.0 );
    dpa3_grow( S.outw, DW );
    const int nm0 = L * C, nmm = M.lmax * C;
    for(int f = 0; f < F; f++)
    {
      const DPA4Focus& Pf = P.foc[f];
      DPA4FocusState& T = A.fs[f];
      if( int(S.XL.size()) < nl + 1 ) S.XL.resize( nl + 1 );
      if( int(T.Y.size()) < nl ) { T.Y.resize( nl ); T.GG.resize( nl ); T.Sd.resize( nl ); }
      for(auto& v : S.XL) dpa3_grow( v, size_t(n) * DmC );
      dpa3_grow( T.xp, DC );
      for(int d = 0; d < D; d++) for(int c = 0; c < C; c++) T.xp[size_t(d)*C + c] = xp[size_t(d)*W + size_t(f)*C + c];
      for(int e = 0; e < n; e++) for(int r = 0; r < Dm; r++)
      {
        const double* a = X + size_t(e) * DmW + size_t(r) * W + size_t(f) * C;
        double* b = S.XL[0].data() + size_t(e) * DmC + size_t(r) * C;
        for(int c = 0; c < C; c++) b[c] = a[c];
      }
      // SO(2) layers
      for(int k = 0; k < nl; k++)
      {
        const double* xk = S.XL[k].data();
        std::vector<double>& Yk = k < nl - 1 ? T.Y[k] : S.Ylast;           // the reverse pass does not read the last one
        dpa3_grow( Yk, size_t(n) * DmC );
        double* y = Yk.data();
        dpa3_gemm( false, false, n, nm0, nm0, 1.0, xk, int(DmC), Pf.lin0[k].data(), nm0, 0.0, y, int(DmC) );
        if( nmm > 0 )
        {
          const double* wm = Pf.linm[k].data();     // [nmm][2 nmm]: u, v
          const double* neg = xk + nm0;
          const double* pos = xk + nm0 + nmm;
          dpa3_gemm( false, false, n, nmm, nmm, 1.0, neg, int(DmC), wm, 2 * nmm, 0.0, y + nm0, int(DmC) );
          dpa3_gemm( false, false, n, nmm, nmm, -1.0, pos, int(DmC), wm + nmm, 2 * nmm, 1.0, y + nm0, int(DmC) );
          dpa3_gemm( false, false, n, nmm, nmm, 1.0, neg, int(DmC), wm + nmm, 2 * nmm, 0.0, y + nm0 + nmm, int(DmC) );
          dpa3_gemm( false, false, n, nmm, nmm, 1.0, pos, int(DmC), wm, 2 * nmm, 1.0, y + nm0 + nmm, int(DmC) );
        }
        double* xo = S.XL[k+1].data();
        if( k < nl - 1 )
        {
          dpa3_grow( T.GG[k], size_t(n) * nmm ); dpa3_grow( T.Sd[k], size_t(n) * C );
          if( nmm > 0 ) dpa3_gemm( false, false, n, nmm, C, 1.0, y, int(DmC), Pf.gate[k].data(), nmm, 0.0, T.GG[k].data(), nmm );
          dpa4_sig_n( T.GG[k].data(), n * nmm, T.GG[k].data() );
          // sigmoid of the scalar rows
          dpa3_grow( S.zb, size_t(n) * C );
          for(int e = 0; e < n; e++) for(int c = 0; c < C; c++) S.zb[size_t(e)*C + c] = y[size_t(e)*DmC + c];
          dpa4_sig_n( S.zb.data(), n * C, S.zb.data() );
          for(int e = 0; e < n; e++)
          {
            const double* ye = y + size_t(e) * DmC;
            const double* xe = xk + size_t(e) * DmC;
            double* oe = xo + size_t(e) * DmC;
            double* sd = T.Sd[k].data() + size_t(e) * C;
            const double* ge = T.GG[k].data() + size_t(e) * nmm;
            const double* zb = S.zb.data() + size_t(e) * C;
            for(int c = 0; c < C; c++)
            {
              const double s = zb[c];
              oe[c] = xe[c] + ye[c] * s;
              sd[c] = s + ye[c] * s * ( 1.0 - s );
            }
            for(int r = 1; r < Dm; r++)
            {
              const double* gg = ge + size_t( M.sign_g[r] ) * C;
              const double* yr = ye + size_t(r) * C;
              const double* xr = xe + size_t(r) * C;
              double* orr = oe + size_t(r) * C;
#             pragma omp simd
              for(int c = 0; c < C; c++) orr[c] = xr[c] + yr[c] * gg[c];
            }
          }
        }
        else for(size_t t = 0; t < size_t(n) * DmC; t++) xo[t] = xk[t] + y[t];
      }
      // attention
      dpa3_grow( S.l0, C ); dpa3_grow( T.q0, C ); dpa3_grow( S.logit, n ); dpa3_grow( T.alpha, n ); dpa3_grow( S.ex, n );
      dpa4_rms( Pf.qk_norm, C, T.xp.data(), S.l0.data() );
      dpa3_gemm( false, false, 1, C, C, 1.0, S.l0.data(), C, Pf.q.data(), C, 0.0, T.q0.data(), C );
      const double ish = 1.0 / std::sqrt( double(C) );
      double mx = Pf.log_null;
      for(int e = 0; e < n; e++)
      {
        const double* kj = BB.Kn.data() + size_t( nbh[e0 + e] ) * W + size_t(f) * C;
        const double* r0 = Rd.data() + size_t(e) * L * W + size_t(f) * C;
        double s = 0.0, b = 0.0;
#       pragma omp simd reduction(+:s,b)
        for(int c = 0; c < C; c++) { s += T.q0[c] * kj[c]; b += r0[c] * Pf.logit_w[c]; }
        const double env = BB.envv[e0 + e];
        const double l = s * ish + b + 2.0 * std::log( env );
        S.logit[e] = l;
        mx = std::max( mx, l );
      }
      double den = std::exp( Pf.log_null - mx );
      for(int e = 0; e < n; e++) S.ex[e] = S.logit[e] - mx;
      dpa4_exp_n( S.ex.data(), n );
      for(int e = 0; e < n; e++) den += S.ex[e];
      for(int e = 0; e < n; e++) T.alpha[e] = S.ex[e] / den;
      dpa3_grow( T.sum, DC ); dpa3_grow( S.out, DC );
      // rotate back with the competition and attention weights: sum_e alpha_e msg_e
      std::fill( T.sum.begin(), T.sum.begin() + DC, 0.0 );
      {
        dpa3_grow( T.XLf, size_t(n) * DmC );
        std::swap( T.XLf, S.XL[nl] );
        const double* xf = T.XLf.data();
        double* sm = T.sum.data();
        for(int e = 0; e < n; e++)
        {
          const double* Dq = BB.Dq.data() + size_t(e0 + e) * M.nDq;
          const double* xl = xf + size_t(e) * DmC;
          const double al = A.alf[size_t(e) * F + f] * T.alpha[e];
          for(int t = 0; t < nrot; t++)
          {
            const int d = M.rot_d[t];
            const double fa = Dq[ M.rot_o[t] ] * M.rescale[d] * al;
            const double* sc = xl + size_t( M.rot_r[t] ) * C;
            double* dc = sm + size_t(d) * C;
#           pragma omp simd
            for(int c = 0; c < C; c++) dc[c] += fa * sc[c];
          }
        }
      }
      dpa4_rms( Pf.gate_norm, C, T.xp.data(), S.tmp.data() );
      {
        double z = 0.0;
        for(int c = 0; c < C; c++) z += S.tmp[c] * Pf.gate_w[c];
        T.g = dpa4_sig( z );
        for(size_t t = 0; t < DC; t++) S.out[t] = T.sum[t] * T.g;
      }
      // message-node grid product
      const int C3 = K * C;
      const size_t DC3 = size_t(D) * C3;
      dpa3_grow( S.h, 2 * DC3 ); dpa3_grow( T.gl, size_t(M.Ng) * C ); dpa3_grow( T.gr, size_t(M.Ng) * C ); dpa3_grow( T.cin, DC3 ); dpa3_grow( S.vv, DC3 );
      dpa3_grow( T.zgate, C ); dpa3_grow( T.pair, 2 * C ); dpa3_grow( S.res, DC );
      dpa4_so3( M.lmax, P.gexp.data(), C, C3, S.out.data(), S.h.data() );
      dpa4_so3( M.lmax, P.gexp.data(), C, C3, T.xp.data(), S.h.data() + DC3 );
      dpa3_gemm( false, false, M.Ng, C, DK, 1.0, M.T.data(), DK, S.h.data(), C, 0.0, T.gl.data(), C );
      dpa3_gemm( false, false, M.Ng, C, DK, 1.0, M.T.data(), DK, S.h.data() + DC3, C, 0.0, T.gr.data(), C );
      dpa3_grow( S.tmp2, size_t(M.Ng) * C );
      for(size_t t = 0; t < size_t(M.Ng) * C; t++) S.tmp2[t] = T.gl[t] * T.gr[t];
      dpa3_gemm( false, false, DK, C, M.Ng, 1.0, M.F.data(), M.Ng, S.tmp2.data(), C, 0.0, T.cin.data(), C );
      for(int c = 0; c < C; c++) { T.pair[c] = S.out[c]; T.pair[C + c] = T.xp[c]; }
      dpa4_scalar( Pf.ggate.data(), C, DK, T.cin.data(), T.pair.data(), S.vv.data(), T.zgate.data() );
      dpa4_so3( M.lmax, P.gcon.data(), C3, C, S.vv.data(), S.res.data() );
      for(int d = 0; d < D; d++) for(int c = 0; c < C; c++) S.outw[size_t(d)*W + size_t(f)*C + c] = S.out[size_t(d)*C + c] + S.res[size_t(d)*C + c] * Pf.gres[c];
    }
    // final mix, norm, residual
    dpa3_grow( A.mixed, DC ); dpa3_grow( A.wn1, C + 1 ); dpa3_grow( A.x1, DC ); dpa3_grow( S.y, DC ); dpa3_grow( A.wn2, C + 1 );
    dpa4_so3( M.lmax, P.post_mix.data(), W, C, S.outw.data(), A.mixed.data() );
    dpa4_eqrms( P.post_norm, M, A.mixed.data(), S.y.data(), A.wn1.data() );
    for(size_t t = 0; t < DC; t++) A.x1[t] = x[t] + S.y[t];
    // FFN
    dpa4_eqrms( P.ffn_norm, M, A.x1.data(), S.y.data(), A.wn2.data() );
    const size_t DCe = size_t(D) * Ce;
    dpa3_grow( A.fgl, size_t(M.Ng) * Ch ); dpa3_grow( A.fgr, size_t(M.Ng) * Ch ); dpa3_grow( S.fc, DCe ); dpa3_grow( A.fco, DCe ); dpa3_grow( S.fv, DCe );
    dpa3_grow( A.fzgate, Ch ); dpa3_grow( A.fpair, 2 * Ch ); dpa3_grow( S.wedge, 2 * DCe ); dpa3_grow( S.tmp3, 2 * DCe ); dpa3_grow( S.tmp4, size_t(M.Ng) * Ch );
    dpa4_so3( M.lmax, P.flin1.data(), C, 2 * Ce, S.y.data(), S.wedge.data() );
    // left, right: rows [D K][Ch], split in the channels of a row of 2 Ce
    for(int d = 0; d < D; d++)
    {
      const double* h = S.wedge.data() + size_t(d) * 2 * Ce;
      for(int t = 0; t < Ce; t++) { S.tmp3[size_t(d)*Ce + t] = h[t]; S.tmp3[DCe + size_t(d)*Ce + t] = h[Ce + t]; }
    }
    for(int c = 0; c < Ch; c++) { A.fpair[c] = S.tmp3[c]; A.fpair[Ch + c] = S.tmp3[DCe + c]; }
    dpa3_gemm( false, false, DK, Ch, Ch, 1.0, S.tmp3.data(), Ch, P.fl.data(), Ch, 0.0, S.fc.data(), Ch );
    dpa3_gemm( false, false, DK, Ch, Ch, 1.0, S.tmp3.data() + DCe, Ch, P.fr.data(), Ch, 0.0, A.fco.data(), Ch );
    dpa3_gemm( false, false, M.Ng, Ch, DK, 1.0, M.T.data(), DK, S.fc.data(), Ch, 0.0, A.fgl.data(), Ch );
    dpa3_gemm( false, false, M.Ng, Ch, DK, 1.0, M.T.data(), DK, A.fco.data(), Ch, 0.0, A.fgr.data(), Ch );
    for(size_t t = 0; t < size_t(M.Ng) * Ch; t++) S.tmp4[t] = A.fgl[t] * A.fgr[t];
    dpa3_gemm( false, false, DK, Ch, M.Ng, 1.0, M.F.data(), M.Ng, S.tmp4.data(), Ch, 0.0, S.fc.data(), Ch );
    dpa3_gemm( false, false, DK, Ch, Ch, 1.0, S.fc.data(), Ch, P.fo.data(), Ch, 0.0, A.fco.data(), Ch );
    dpa4_scalar( P.fgate.data(), Ch, DK, A.fco.data(), A.fpair.data(), S.fv.data(), A.fzgate.data() );
    dpa4_so3( M.lmax, P.flin2.data(), Ce, C, S.fv.data(), x2 );
    for(size_t t = 0; t < DC; t++) x2[t] += A.x1[t];
  }

  // reverse of dpa4_block_forward: Gout[D][C] dE/dx2 -> dx1 (dE/dx, node term), source atoms and edge accumulators
  inline void dpa4_block_vjp(const DPA4Model& M, const DPA4Block& P, DPA4Buffers& BB, int i, const int64_t e0, const int n,
                             const int* nbh, const double* Gout, DPA4AtomState& A, DPA4Scratch& S, double* dxnode)
  {
    const int C = M.C, D = M.D, Dm = M.Dm, L = M.L, DK = D * M.K, K = M.K, Ch = M.Ch, Ce = M.Ce, F = M.F_n, W = M.W;
    const size_t DC = size_t(D) * C, DmC = size_t(Dm) * C, DCe = size_t(D) * Ce, DW = size_t(D) * W, DmW = size_t(Dm) * W;
    const double* Xp = BB.Xp.data();
    DPA3Buffers& B = BB.c;
    dpa3_grow( S.dX, DC ); dpa3_grow( S.dxp, DW );
    std::fill( S.dxp.begin(), S.dxp.begin() + DW, 0.0 );
    // FFN: x2 = x1 + lin2( scalar( o grid( l, r ) ) )
    dpa3_grow( S.dc, 2 * DCe ); dpa3_grow( S.dp, size_t(M.Ng) * Ch ); dpa3_grow( S.dleft, DCe ); dpa3_grow( S.dright, DCe ); dpa3_grow( S.dh, 2 * DCe );
    dpa3_grow( S.dpairv, 2 * Ch ); dpa3_grow( S.df, Ch ); dpa3_grow( S.tmp2, DCe ); dpa3_grow( S.dG, size_t(M.Ng) * Ch );
    dpa4_so3_t( M.lmax, P.flin2.data(), Ce, C, Gout, S.dc.data() );
    dpa4_scalar_t( P.fgate.data(), Ch, DK, A.fco.data(), A.fpair.data(), A.fzgate.data(), S.dc.data(), S.tmp2.data(), S.dpairv.data(), S.df.data() );
    // fco = fc o: dfc = dfco o^T
    dpa3_gemm( false, true, DK, Ch, Ch, 1.0, S.tmp2.data(), Ch, P.fo.data(), Ch, 0.0, S.dright.data(), Ch );
    // fc = F p: dp = F^T dfc
    dpa3_gemm( true, false, M.Ng, Ch, DK, 1.0, M.F.data(), M.Ng, S.dright.data(), Ch, 0.0, S.dp.data(), Ch );
    for(size_t t = 0; t < size_t(M.Ng) * Ch; t++) S.dG[t] = S.dp[t] * A.fgr[t];
    dpa3_gemm( true, false, DK, Ch, M.Ng, 1.0, M.T.data(), DK, S.dG.data(), Ch, 0.0, S.dright.data(), Ch );
    dpa3_gemm( false, true, DK, Ch, Ch, 1.0, S.dright.data(), Ch, P.fl.data(), Ch, 0.0, S.dleft.data(), Ch );
    for(size_t t = 0; t < size_t(M.Ng) * Ch; t++) S.dG[t] = S.dp[t] * A.fgl[t];
    dpa3_gemm( true, false, DK, Ch, M.Ng, 1.0, M.T.data(), DK, S.dG.data(), Ch, 0.0, S.tmp2.data(), Ch );
    dpa3_gemm( false, true, DK, Ch, Ch, 1.0, S.tmp2.data(), Ch, P.fr.data(), Ch, 0.0, S.dright.data(), Ch );
    for(int c = 0; c < Ch; c++) { S.dleft[c] += S.dpairv[c]; S.dright[c] += S.dpairv[Ch + c]; }
    for(int d = 0; d < D; d++)
    {
      double* h = S.dh.data() + size_t(d) * 2 * Ce;
      for(int t = 0; t < Ce; t++) { h[t] = S.dleft[size_t(d)*Ce + t]; h[Ce + t] = S.dright[size_t(d)*Ce + t]; }
    }
    dpa3_grow( S.dY, DC );
    dpa4_so3_t( M.lmax, P.flin1.data(), C, 2 * Ce, S.dh.data(), S.dY.data() );
    for(size_t t = 0; t < DC; t++) S.dX[t] = Gout[t];
    dpa4_eqrms_t( P.ffn_norm, M, A.x1.data(), A.wn2.data(), S.dY.data(), S.dX.data() );
    // dE/dx1: node term of dE/dx; through the norm and the final mix to dE/dout2
    for(size_t t = 0; t < DC; t++) dxnode[t] = S.dX[t];
    std::fill( S.dY.begin(), S.dY.begin() + DC, 0.0 );
    dpa4_eqrms_t( P.post_norm, M, A.mixed.data(), A.wn1.data(), S.dX.data(), S.dY.data() );
    dpa3_grow( S.doutw, DW );
    dpa4_so3_t( M.lmax, P.post_mix.data(), W, C, S.dY.data(), S.doutw.data() );
    const int nl = P.foc[0].lin0.size(), nrot = M.rot_r.size();
    const int nm0 = L * C, nmm = M.lmax * C;
    dpa3_grow( S.dXW, size_t(n) * DmW ); dpa3_grow( S.Grw, size_t(n) * L * W ); dpa3_grow( S.dAf, size_t(n) * F );
    std::fill( S.Grw.begin(), S.Grw.begin() + size_t(n) * L * W, 0.0 );
    std::fill( S.dAf.begin(), S.dAf.begin() + size_t(n) * F, 0.0 );
    for(int f = 0; f < F; f++)
    {
      const DPA4Focus& Pf = P.foc[f];
      DPA4FocusState& T = A.fs[f];
      dpa3_grow( S.dout, DC ); dpa3_grow( S.dres, DC ); dpa3_grow( S.dxf, DC );
      for(int d = 0; d < D; d++) for(int c = 0; c < C; c++) S.dout[size_t(d)*C + c] = S.doutw[size_t(d)*W + size_t(f)*C + c];
      // grid product: res = scale contract( scalar( F( T L * T R ) ) )
      const int C3 = K * C;
      const size_t DC3 = size_t(D) * C3;
      for(int d = 0; d < D; d++) for(int c = 0; c < C; c++) S.dres[size_t(d)*C + c] = S.dout[size_t(d)*C + c] * Pf.gres[c];
      dpa3_grow( S.dc, 2 * DC3 ); dpa3_grow( S.dz, DC3 ); dpa3_grow( S.drad, size_t(M.Ng) * C );
      dpa3_grow( S.dpairv, 2 * C ); dpa3_grow( S.df, C );
      dpa4_so3_t( M.lmax, P.gcon.data(), C3, C, S.dres.data(), S.dc.data() );
      dpa4_scalar_t( Pf.ggate.data(), C, DK, T.cin.data(), T.pair.data(), T.zgate.data(), S.dc.data(), S.dz.data(), S.dpairv.data(), S.df.data() );
      dpa3_gemm( true, false, M.Ng, C, DK, 1.0, M.F.data(), M.Ng, S.dz.data(), C, 0.0, S.drad.data(), C );
      dpa3_grow( S.dG, size_t(M.Ng) * C );
      for(size_t t = 0; t < size_t(M.Ng) * C; t++) S.dG[t] = S.drad[t] * T.gr[t];
      dpa3_gemm( true, false, DK, C, M.Ng, 1.0, M.T.data(), DK, S.dG.data(), C, 0.0, S.dz.data(), C );
      dpa4_so3_t( M.lmax, P.gexp.data(), C, C3, S.dz.data(), S.dout.data(), 1.0 );    // dE/dout
      for(size_t t = 0; t < size_t(M.Ng) * C; t++) S.dG[t] = S.drad[t] * T.gl[t];
      dpa3_gemm( true, false, DK, C, M.Ng, 1.0, M.T.data(), DK, S.dG.data(), C, 0.0, S.dz.data(), C );
      dpa4_so3_t( M.lmax, P.gexp.data(), C, C3, S.dz.data(), S.dxf.data() );          // dE/dxp_i
      for(int c = 0; c < C; c++) { S.dout[c] += S.dpairv[c]; S.dxf[c] += S.dpairv[C + c]; }
      // out = gate sum_e alpha_e msg_e
      const double g = T.g;
      double dg = 0.0;
      for(size_t t = 0; t < DC; t++) dg += S.dout[t] * T.sum[t];
      {
        // gate = sig( rms(xp0) . w )
        dpa3_grow( S.tmp, C );
        const double fg = dg * g * ( 1.0 - g );
        for(int c = 0; c < C; c++) S.tmp[c] = fg * Pf.gate_w[c];
        dpa4_rms_t( Pf.gate_norm, C, T.xp.data(), S.tmp.data(), S.dxf.data() );
      }
      dpa3_grow( S.dalpha, n );
      // messages -> local frame, packed D gradients, gradient on the competition weight
      dpa3_grow( S.dX, size_t(n) * DmC );
      std::fill( S.dX.begin(), S.dX.begin() + size_t(n) * DmC, 0.0 );
      {
        const double* xf = T.XLf.data();
        for(int e = 0; e < n; e++)
        {
          const double* Dq = BB.Dq.data() + size_t(e0 + e) * M.nDq;
          const double* dm = S.dout.data();                  // dE/dmsg_e = alpha_e g dout
          const double* xl = xf + size_t(e) * DmC;
          double* dl = S.dX.data() + size_t(e) * DmC;
          double* gd = BB.GD.data() + size_t(e0 + e) * M.nDq;
          const double sc = T.alpha[e] * g;
          const double al = A.alf[size_t(e) * F + f];
          double da = 0.0, dae = 0.0;
          for(int t = 0; t < nrot; t++)
          {
            const int d = M.rot_d[t];
            const int r = M.rot_r[t];
            const double f0 = Dq[ M.rot_o[t] ] * M.rescale[d];
            const double fa = f0 * al * sc;
            const double* dmc = dm + size_t(d) * C;
            const double* xr = xl + size_t(r) * C;
            double* dr = dl + size_t(r) * C;
            double ss = 0.0;
#           pragma omp simd reduction(+:ss)
            for(int c = 0; c < C; c++) { dr[c] += fa * dmc[c]; ss += dmc[c] * xr[c]; }
            gd[ M.rot_o[t] ] += ss * M.rescale[d] * al * sc;
            da += ss * f0 * sc;
            dae += ss * f0 * al;
          }
          S.dAf[size_t(e) * F + f] = da;
          S.dalpha[e] = dae * g;
        }
      }
      // alpha: dE/dalpha_e = g dout . msg_e, dE/dlogit_e = alpha_e ( A_e - sum alpha A )
      double mean = 0.0;
      for(int e = 0; e < n; e++) mean += S.dalpha[e] * T.alpha[e];
      const double ish = 1.0 / std::sqrt( double(C) );
      std::fill( S.tmp.begin(), S.tmp.begin() + C, 0.0 );   // dE/dq_i
      for(int e = 0; e < n; e++)
      {
        const double dl = T.alpha[e] * ( S.dalpha[e] - mean );
        const int j = nbh[e0 + e];
        const double* kj = BB.Kn.data() + size_t(j) * W + size_t(f) * C;
        for(int c = 0; c < C; c++) S.tmp[c] += dl * ish * kj[c];
        {
          DPA3RowLock lk( &B.locks[j] );
          double* gk = BB.GKn.data() + size_t(j) * W + size_t(f) * C;
          for(int c = 0; c < C; c++) gk[c] += dl * ish * T.q0[c];
        }
        double* gr = S.Grw.data() + size_t(e) * L * W + size_t(f) * C;
        for(int c = 0; c < C; c++) gr[c] += dl * Pf.logit_w[c];
        BB.Genv[e0 + e] += dl * 2.0 / BB.envv[e0 + e];
      }
      // dE/dq_i -> dE/dxp_i[0]
      dpa3_grow( S.ddso, C );
      dpa3_gemm( false, true, 1, C, C, 1.0, S.tmp.data(), C, Pf.q.data(), C, 0.0, S.ddso.data(), C );
      dpa4_rms_t( Pf.qk_norm, C, T.xp.data(), S.ddso.data(), S.dxf.data() );
      for(int d = 0; d < D; d++) for(int c = 0; c < C; c++) S.dxp[size_t(d)*W + size_t(f)*C + c] += S.dxf[size_t(d)*C + c];
      // SO(2) layers backward
      dpa3_grow( S.dY, size_t(n) * DmC ); dpa3_grow( S.dz, size_t(n) * DmC ); dpa3_grow( S.dG, size_t(n) * ( nmm > 0 ? nmm : 1 ) );
      for(int k = nl - 1; k >= 0; k--)
      {
        const double* y = T.Y[k].data();
        double* dy = S.dY.data();
        if( k < nl - 1 )
        {
          const double* gg = T.GG[k].data();
          std::fill( S.dG.begin(), S.dG.begin() + size_t(n) * nmm, 0.0 );
          for(int e = 0; e < n; e++)
          {
            const double* dxo = S.dX.data() + size_t(e) * DmC;
            const double* ye = y + size_t(e) * DmC;
            const double* ge = gg + size_t(e) * nmm;
            double* de = dy + size_t(e) * DmC;
            double* dg0 = S.dG.data() + size_t(e) * nmm;
            const double* sd = T.Sd[k].data() + size_t(e) * C;
            for(int c = 0; c < C; c++) de[c] = dxo[c] * sd[c];
            for(int r = 1; r < Dm; r++)
            {
              const int gi = M.sign_g[r];
              const double* gr = ge + size_t(gi) * C;
              const double* dr = dxo + size_t(r) * C;
              const double* yr = ye + size_t(r) * C;
              double* er = de + size_t(r) * C;
              double* dgg = dg0 + size_t(gi) * C;
#             pragma omp simd
              for(int c = 0; c < C; c++) { er[c] = dr[c] * gr[c]; dgg[c] += dr[c] * yr[c]; }
            }
          }
          // dE/dz = dg sig'(z), to the scalar row of y
          for(size_t t = 0; t < size_t(n) * nmm; t++) S.dG[t] *= gg[t] * ( 1.0 - gg[t] );
          if( nmm > 0 ) dpa3_gemm( false, true, n, C, nmm, 1.0, S.dG.data(), nmm, Pf.gate[k].data(), nmm, 1.0, dy, int(DmC) );
        }
        else for(size_t t = 0; t < size_t(n) * DmC; t++) dy[t] = S.dX[t];
        // dx_k = dx_{k+1} + dy W^T
        double* dx = S.dz.data();
        for(size_t t = 0; t < size_t(n) * DmC; t++) dx[t] = S.dX[t];
        dpa3_gemm( false, true, n, nm0, nm0, 1.0, dy, int(DmC), Pf.lin0[k].data(), nm0, 1.0, dx, int(DmC) );
        if( nmm > 0 )
        {
          const double* wm = Pf.linm[k].data();
          const double* dneg = dy + nm0;
          const double* dpos = dy + nm0 + nmm;
          // neg' = neg u - pos v, pos' = neg v + pos u, weights [nmm][2 nmm]
          std::vector<double>& U = S.tmp3; dpa3_grow( U, size_t(nmm) * nmm );
          std::vector<double>& V = S.tmp4; dpa3_grow( V, size_t(nmm) * nmm );
          for(int a = 0; a < nmm; a++) for(int b = 0; b < nmm; b++) { U[size_t(a)*nmm + b] = wm[size_t(a)*2*nmm + b]; V[size_t(a)*nmm + b] = wm[size_t(a)*2*nmm + nmm + b]; }
          dpa3_gemm( false, true, n, nmm, nmm, 1.0, dneg, int(DmC), U.data(), nmm, 1.0, dx + nm0, int(DmC) );
          dpa3_gemm( false, true, n, nmm, nmm, 1.0, dpos, int(DmC), V.data(), nmm, 1.0, dx + nm0, int(DmC) );
          dpa3_gemm( false, true, n, nmm, nmm, -1.0, dneg, int(DmC), V.data(), nmm, 1.0, dx + nm0 + nmm, int(DmC) );
          dpa3_gemm( false, true, n, nmm, nmm, 1.0, dpos, int(DmC), U.data(), nmm, 1.0, dx + nm0 + nmm, int(DmC) );
        }
        S.dX.swap( S.dz );
      }
      // dE/dXW, slice of this focus
      for(int e = 0; e < n; e++) for(int r = 0; r < Dm; r++)
      {
        const double* a = S.dX.data() + size_t(e) * DmC + size_t(r) * C;
        double* b = S.dXW.data() + size_t(e) * DmW + size_t(r) * W + size_t(f) * C;
        for(int c = 0; c < C; c++) b[c] = a[c];
      }
    }
    // cross-focus competition backward, to the scalar row of XW
    if( F > 1 && P.compete )
    {
      dpa3_grow( S.tmp, C ); dpa3_grow( S.ddso, C );
      const double* X = A.XW0.data();
      for(int e = 0; e < n; e++)
      {
        const double* a = A.af.data() + size_t(e) * F;
        const double* dA = S.dAf.data() + size_t(e) * F;
        double mean = 0.0;
        for(int f = 0; f < F; f++) mean += a[f] * dA[f];
        for(int f = 0; f < F; f++)
        {
          const double dl = ( 1.0 - M.focus_ls ) * a[f] * ( dA[f] - mean ) / M.focus_tau;
          for(int c = 0; c < C; c++) S.ddso[c] = dl * P.foc[f].comp_w[c];
          dpa4_rms_t( P.foc[f].comp_norm, C, X + size_t(e) * W + size_t(f) * C, S.ddso.data(), S.dXW.data() + size_t(e) * DmW + size_t(f) * C );
        }
      }
    }
    // radial mixer backward: dXW is dE/dXW
    dpa3_grow( S.dX0, size_t(n) * DmW );
    double* dx0 = S.dX0.data();
    if( M.nk > 0 )
    {
      dpa3_grow( S.dG, size_t(n) * M.nk );
      std::fill( dx0, dx0 + size_t(n) * DmW, 0.0 );
      std::fill( S.dG.begin(), S.dG.begin() + size_t(n) * M.nk, 0.0 );
      const int nm = M.mix_r.size();
      for(int e = 0; e < n; e++)
      {
        const double* k = A.kern.data() + size_t(e) * M.nk;
        const double* x0 = A.X0.data() + size_t(e) * DmW;
        const double* dxo = S.dXW.data() + size_t(e) * DmW;
        double* dxi = dx0 + size_t(e) * DmW;
        double* dk = S.dG.data() + size_t(e) * M.nk;
        for(int t = 0; t < nm; t++)
        {
          const int cc = M.mix_c[t];
          const double fk = k[cc];
          const double* dr = dxo + size_t( M.mix_r[t] ) * W;
          const double* xs = x0 + size_t( M.mix_s[t] ) * W;
          double* ds = dxi + size_t( M.mix_s[t] ) * W;
          double ss = 0.0;
#         pragma omp simd reduction(+:ss)
          for(int c = 0; c < W; c++) { const double db = dr[c] * P.mixb[c]; ds[c] += fk * db; ss += db * xs[c]; }
          dk[cc] += ss;
        }
      }
      dpa3_gemm( false, true, n, L * W, M.nk, 1.0, S.dG.data(), M.nk, P.mixw.data(), M.nk, 1.0, S.Grw.data(), L * W );
    }
    else
    {
      for(int e = 0; e < n; e++)
      {
        const double* x0 = A.X0.data() + size_t(e) * DmW;
        const double* rd = A.rad.data() + size_t(e) * L * W;
        const double* dxo = S.dXW.data() + size_t(e) * DmW;
        double* dxi = dx0 + size_t(e) * DmW;
        double* gr = S.Grw.data() + size_t(e) * L * W;
        for(int r = 0; r < Dm; r++)
        {
          const int l = M.deg_m[r];
          for(int c = 0; c < W; c++) { dxi[size_t(r)*W + c] = dxo[size_t(r)*W + c] * rd[size_t(l)*W + c]; gr[size_t(l)*W + c] += dxo[size_t(r)*W + c] * x0[size_t(r)*W + c]; }
        }
      }
    }
    // hidden radial projection backward, to the radial features of the edges
    if( F > 1 ) dpa3_gemm( false, true, n * L, C, W, 1.0, S.Grw.data(), W, P.rproj.data(), W, 1.0, BB.Grf.data() + size_t(e0) * L * C, C );
    else for(size_t t = 0; t < size_t(n) * L * C; t++) BB.Grf[ size_t(e0) * L * C + t ] += S.Grw[t];
    // rotation to the edge frame: dE/dxp of the source atoms, packed D gradients
    dpa3_grow( S.dres, DW );
    for(int e = 0; e < n; e++)
    {
      const int j = nbh[e0 + e];
      const double* xs = Xp + size_t(j) * DW;
      const double* Dq = BB.Dq.data() + size_t(e0 + e) * M.nDq;
      const double* dl = dx0 + size_t(e) * DmW;
      double* gd = BB.GD.data() + size_t(e0 + e) * M.nDq;
      std::fill( S.dres.begin(), S.dres.begin() + DW, 0.0 );
      for(int t = 0; t < nrot; t++)
      {
        const int d = M.rot_d[t];
        const double fq = Dq[ M.rot_o[t] ];
        const double* dr = dl + size_t( M.rot_r[t] ) * W;
        const double* xc = xs + size_t(d) * W;
        double* gc = S.dres.data() + size_t(d) * W;
        double ss = 0.0;
#       pragma omp simd reduction(+:ss)
        for(int c = 0; c < W; c++) { gc[c] += fq * dr[c]; ss += dr[c] * xc[c]; }
        gd[ M.rot_o[t] ] += ss;
      }
      DPA3RowLock lk( &B.locks[j] );
      double* gx = BB.GXp.data() + size_t(j) * DW;
      for(size_t t = 0; t < DW; t++) gx[t] += S.dres[t];
    }
    {
      DPA3RowLock lk( &B.locks[i] );
      double* gx = BB.GXp.data() + size_t(i) * DW;
      for(size_t t = 0; t < DW; t++) gx[t] += S.dxp[t];
    }
  }

  // ------------------------------------------------------------------ whole system

  template<class GridT>
  inline void dpa4_compute( GridT& grid, const GridChunkNeighbors& chunk_neighbors, const Domain& domain, const DPA4Model& M,
                            const DeepmdModelView& fit, const GhostCommunicationScheme* scheme, MPI_Comm comm, DPA4Buffers& BB )
  {
    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;
    static constexpr double conv_energy = EXASTAMP_CONST_QUANTITY( 1. * eV );
    DPA3Buffers& B = BB.c;
    const size_t ncells = grid.number_of_cells();
    const size_t* poff = grid.cell_particle_offset_data();
    const size_t n_nodes = grid.number_of_particles();
    auto* cells = grid.cells();
    const IJK dims = grid.dimension();
    const Mat3d xform = domain.xform();
    const bool identity = domain.xform_is_identity();
    const double rc2 = M.rc * M.rc;
    const int C = M.C, D = M.D, L = M.L, NB = M.blocks.size();
    const size_t DC = size_t(D) * C;
    const int F = M.F_n, W = M.W;
    const size_t DW = size_t(D) * W;
    DPA4Prof& prof = dpa4_prof();
    double tp = prof.on ? omp_get_wtime() : 0.0;
    auto mark = [&](int k) { if( prof.on ) { const double t = omp_get_wtime(); prof.t[k] += t - tp; tp = t; } };
    prof.evals ++;

    // CSR bonds within the cutoff, nearest nnei first
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
    int* nbh = B.nbh.data();
    double* rij = B.rij.data();
    const int64_t* offset = B.offset.data();
    // more than nnei neighbors: keep the nearest
    BB.ecount.resize( n_nodes );
#   pragma omp parallel
    {
      std::vector<int> order;
      std::vector<int> nb2;
      std::vector<double> rr2;
#     pragma omp for schedule(dynamic, 16)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q];
        const int64_t b0 = offset[i];
        const int cnt = offset[i+1] - b0;
        BB.ecount[i] = cnt;
        if( cnt <= M.nnei ) continue;
        order.resize( cnt );
        for(int k = 0; k < cnt; k++) order[k] = k;
        std::stable_sort( order.begin(), order.end(), [&](int x, int y)
        {
          const double* a = rij + 3 * ( b0 + x ); const double* c = rij + 3 * ( b0 + y );
          return a[0]*a[0] + a[1]*a[1] + a[2]*a[2] < c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
        } );
        nb2.resize( cnt ); rr2.resize( 3 * cnt );
        for(int k = 0; k < cnt; k++) { nb2[k] = nbh[b0 + order[k]]; for(int x = 0; x < 3; x++) rr2[3*k + x] = rij[3*( b0 + order[k] ) + x]; }
        for(int k = 0; k < cnt; k++) { nbh[b0 + k] = nb2[k]; for(int x = 0; x < 3; x++) rij[3*( b0 + k ) + x] = rr2[3*k + x]; }
        BB.ecount[i] = M.nnei;
      }
    }
    auto exchange = [&](std::vector<double>& buf, bool to_ghost, size_t row = 0)
    {
      if( row == 0 ) row = DC;
      if( scheme == nullptr ) fatal_error() << "DPA4: no ghost communication scheme (a ghost update must run before deepmd_force)" << std::endl;
      GhostFlatBuffer fb;
      fb.m_data = buf.data(); fb.m_row_size = row; fb.m_row_stride = row; fb.m_nrows = n_nodes; fb.m_type = GhostFlatBuffer::ValueType::Double;
      if( to_ghost ) flat_buffer_ghost_exchange<true>( comm, *scheme, poff, fb, 331 );
      else flat_buffer_ghost_exchange<false>( comm, *scheme, poff, fb, 332 );
    };

    // per edge geometry, radial features, initial node features
    BB.envv.assign( nbonds, 0.0 ); BB.rr.assign( nbonds, 0.0 ); BB.rbf.assign( nbonds * M.nrad, 0.0 ); BB.Dq.assign( nbonds * M.nDq, 0.0 );
    BB.rf.assign( nbonds * L * C, 0.0 );
    BB.GD.assign( nbonds * M.nDq, 0.0 ); BB.Grf.assign( nbonds * L * C, 0.0 ); BB.Genv.assign( nbonds, 0.0 ); BB.Grbf.assign( nbonds * M.nrad, 0.0 ); BB.Gv.assign( nbonds * 3, 0.0 );
    BB.isd.assign( n_nodes, 0.0 ); BB.deg.assign( n_nodes, 0.0 );
    BB.X.resize( NB + 1 ); BB.G.resize( NB + 1 );
    for(auto& v : BB.X) v.assign( n_nodes * DC, 0.0 );
    const int Sc = M.s_r + 2 * M.tdim, Em = M.emb, Ax = M.axis;
    BB.seed.assign( n_nodes * size_t(Em) * 4 , 0.0 );      // env seed aggregate (normalized) [4][emb]
    BB.gie.assign( n_nodes * DC, 0.0 );                    // unnormalized initial l > 0 sums
#   pragma omp parallel
    {
      std::vector<double> Dd( size_t(Em) * Ax ), rt(4), film( 2 * C );
      std::vector<double> filmn( 2 * C ), hb, shb, f1b, f2b, ginb, ghb, gib;
#     pragma omp for schedule(dynamic, 8)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q];
        const int64_t b0 = offset[i];
        const int n = BB.ecount[i];
        double deg = 0.0;
        for(int k = 0; k < n; k++)
        {
          const int64_t e = b0 + k;
          const double* v = rij + 3 * e;
          const double r = std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + M.eps * M.eps );
          BB.rr[e] = r;
          double en, de;
          dpa4_envelope( M.env_c, r, M.rc, en, de );
          BB.envv[e] = en;
          deg += en * en;
          double e7, d7;
          dpa4_envelope( M.rbf_c, r, M.rc, e7, d7 );
          double* rb = BB.rbf.data() + size_t(e) * M.nrad;
          if( M.gauss != 0.0 ) for(int a = 0; a < M.nrad; a++) { const double dr = r - M.freq[a]; rb[a] = std::exp( M.gauss * dr * dr ) * e7; }
          else for(int a = 0; a < M.nrad; a++) { const double x = r * M.freq[a]; rb[a] = M.freq[a] * std::sin(x) / x * e7; }
          double qq[4], J[12];
          dpa4_quaternion( v, r, M.eps, M.weps, qq, J );
          dpa4_wigner( M, qq, BB.Dq.data() + size_t(e) * M.nDq );
        }
        // radial features * envelope, all the edges of the atom
        const int hr = M.hrf;
        const double* Rb = BB.rbf.data() + size_t(b0) * M.nrad;
        dpa3_grow( hb, size_t(n) * hr ); dpa3_grow( shb, size_t(n) * hr );
        dpa3_gemm( false, false, n, hr, M.nrad, 1.0, Rb, M.nrad, M.rf1.data(), hr, 0.0, hb.data(), hr );
        for(int k = 0; k < n; k++) dpa4_rms( M.rf_norm, hr, hb.data() + size_t(k) * hr, shb.data() + size_t(k) * hr );
        dpa4_silu_n( shb.data(), n * hr );
        for(int k = 0; k < n; k++) { const double en = BB.envv[b0 + k]; double* r = shb.data() + size_t(k) * hr; for(int a = 0; a < hr; a++) r[a] *= en; }
        dpa3_gemm( false, false, n, L * C, hr, 1.0, shb.data(), hr, M.rf2.data(), L * C, 0.0, BB.rf.data() + size_t(b0) * L * C, L * C );
        for(int k = 0; k < n; k++)
        {
          const int64_t e = b0 + k;
          // geometric initial embedding
          const int j = nbh[e];
          const double* rfe = BB.rf.data() + size_t(e) * L * C;
          const double* Dqe = BB.Dq.data() + size_t(e) * M.nDq;
          double* gie = BB.gie.data() + size_t(i) * DC;
          for(int t = 0; t < D - 1; t++)
          {
            const double z = Dqe[ M.z_col[t] ];
            const double* rv = rfe + size_t( M.z_slot[t] ) * C;
            double* gg = gie + size_t( t + 1 ) * C;
            for(int c = 0; c < C; c++) gg[c] += z * rv[c];
          }
          // env seed
          (void) j;
        }
        BB.deg[i] = deg;
        const double isd = 1.0 / std::sqrt( deg + M.deg_floor );
        BB.isd[i] = isd;
        double* ga = BB.seed.data() + size_t(i) * Em * 4;
        // rbf proj (activations in M.s_a), type embeddings, all the edges of the atom
        dpa3_grow( f1b, size_t(n) * M.s_r ); dpa3_grow( f2b, size_t(n) * M.s_r ); dpa3_grow( ginb, size_t(n) * Sc ); dpa3_grow( ghb, size_t(n) * M.s_gh ); dpa3_grow( gib, size_t(n) * Em );
        dpa3_gemm( false, false, n, M.s_r, M.nrad, 1.0, Rb, M.nrad, M.s_r1.data(), M.s_r, 0.0, f1b.data(), M.s_r );
        if( M.s_a[0] ) dpa4_silu_n( f1b.data(), n * M.s_r );
        dpa3_gemm( false, false, n, M.s_r, M.s_r, 1.0, f1b.data(), M.s_r, M.s_r2.data(), M.s_r, 0.0, f2b.data(), M.s_r );
        if( M.s_a[1] ) dpa4_silu_n( f2b.data(), n * M.s_r );
        {
          const double* td = M.tt.data() + size_t( B.type[i] ) * M.tdim;
          for(int k = 0; k < n; k++)
          {
            double* g = ginb.data() + size_t(k) * Sc;
            const double* f2 = f2b.data() + size_t(k) * M.s_r;
            const double* ts = M.tt.data() + size_t( B.type[ nbh[b0 + k] ] ) * M.tdim;
            for(int a = 0; a < M.s_r; a++) g[a] = f2[a];
            for(int a = 0; a < M.tdim; a++) { g[M.s_r + a] = ts[a]; g[M.s_r + M.tdim + a] = td[a]; }
          }
        }
        dpa3_gemm( false, false, n, M.s_gh, Sc, 1.0, ginb.data(), Sc, M.s_g1.data(), M.s_gh, 0.0, ghb.data(), M.s_gh );
        if( M.s_a[2] ) dpa4_silu_n( ghb.data(), n * M.s_gh );
        dpa3_gemm( false, false, n, Em, M.s_gh, 1.0, ghb.data(), M.s_gh, M.s_g2.data(), Em, 0.0, gib.data(), Em );
        if( M.s_a[3] ) dpa4_silu_n( gib.data(), n * Em );
        for(int k = 0; k < n; k++)
        {
          const int64_t e = b0 + k;
          const double* v = rij + 3 * e;
          const double r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
          const double inv = 1.0 / std::sqrt( r2 + M.seed_eps * M.seed_eps );
          const double s = BB.envv[e] * inv;
          rt[0] = s; rt[1] = s * v[0] * inv; rt[2] = s * v[1] * inv; rt[3] = s * v[2] * inv;
          const double* gi = gib.data() + size_t(k) * Em;
          for(int x = 0; x < 4; x++) for(int a = 0; a < Em; a++) ga[x*Em + a] += rt[x] * gi[a] * isd;
        }
        // Dd[Em][Ax] = agg^T agg[:, :Ax]
        dpa3_gemm( true, false, Em, Ax, 4, 1.0, ga, Em, ga, Em, 0.0, Dd.data(), Ax );
        dpa3_gemm( false, false, 1, 2 * C, Em * Ax, 1.0, Dd.data(), Em * Ax, M.s_out.data(), 2 * C, 0.0, film.data(), 2 * C );
        if( M.s_a[4] ) dpa4_silu_n( film.data(), 2 * C );
        dpa4_rms( M.film_s, C, film.data(), filmn.data() );
        dpa4_rms( M.film_f, C, film.data() + C, filmn.data() + C );
        double* x = BB.X[0].data() + size_t(i) * DC;
        const double* te = M.te.data() + size_t( B.type[i] ) * C;
        for(int c = 0; c < C; c++) x[c] = te[c] * ( 1.0 + M.film_scale * std::tanh( filmn[c] ) ) + M.film_shift * std::tanh( filmn[C + c] );
        const double* gs = BB.gie.data() + size_t(i) * DC;
        for(size_t t = C; t < DC; t++) x[t] = gs[t] * isd;
      }
    }
    exchange( BB.X[0], true );
    mark(0);

    auto node_prep = [&](int b)
    {
      const DPA4Block& P = M.blocks[b];
      BB.Xp.resize( n_nodes * DW ); BB.Kn.resize( n_nodes * W );
#     pragma omp parallel
      {
        std::vector<double> l0(C), yn(DC), wn(C + 1);
#       pragma omp for schedule(static)
        for(size_t n = 0; n < n_nodes; n++)
        {
          const double* x = BB.X[b].data() + n * DC;
          double* xp = BB.Xp.data() + n * DW;
          dpa4_eqrms( P.pre_norm, M, x, yn.data(), wn.data() );
          dpa4_so3( M.lmax, P.pre_mix.data(), C, W, yn.data(), xp );
          for(int f = 0; f < F; f++)
          {
            dpa4_rms( P.foc[f].qk_norm, C, xp + size_t(f) * C, l0.data() );
            dpa3_gemm( false, false, 1, C, C, 1.0, l0.data(), C, P.foc[f].k.data(), C, 0.0, BB.Kn.data() + n * W + size_t(f) * C, C );
          }
        }
      }
    };

    // forward blocks
    if( M.store_forward ) { BB.store.resize( NB ); for(auto& v : BB.store) v.resize( inum ); }
    for(int b = 0; b < NB; b++)
    {
      node_prep( b );
#     pragma omp parallel
      {
        static thread_local DPA4Scratch S;
#       pragma omp for schedule(dynamic, 1)
        for(int q = 0; q < inum; q++)
        {
          const int i = B.ilist[q];
          DPA4AtomState& A = M.store_forward ? BB.store[b][q] : S.own;
          dpa4_block_forward( M, M.blocks[b], BB, i, offset[i], BB.ecount[i], nbh, BB.X[b].data() + size_t(i) * DC, A, S, BB.X[b+1].data() + size_t(i) * DC );
        }
      }
      if( b + 1 < NB ) exchange( BB.X[b+1], true );
    }

    mark(1);
    // read-out, fitting
    B.eatom.assign( n_nodes, 0.0 );
    BB.G[NB].assign( n_nodes * DC, 0.0 );
    const int Ce = M.Ce, Ch = M.Ch, K = M.K, DK = D * K;
    const size_t DCe = size_t(D) * Ce;
#   pragma omp parallel
    {
      std::vector<double> h( 2 * DCe ), fus( size_t(DK) * 2 * Ch ), ll( size_t(DK) * M.rhid ), rr( size_t(DK) * M.rhid ), s( M.rhid ), co( Ch ), pair( 2 * Ch ), zg( Ch ), z( Ch ), y( C ), desc( C );
      std::vector<double> yh( 2 * M.fit_h ), v( M.fit_h ), dv( M.fit_h ), dy( 2 * M.fit_h ), dd( C ), dz( Ch ), dco( Ch ), ds( M.rhid ), dgt( Ch ), dp( 2 * Ch ), dl( size_t(DK) * M.rhid ), dr( size_t(DK) * M.rhid ), dfus( size_t(DK) * 2 * Ch ), dh( 2 * DCe );
#     pragma omp for schedule(static)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q], ti = B.type[i];
        const double* x = BB.X[NB].data() + size_t(i) * DC;
        dpa4_so3( M.lmax, M.rlin1.data(), C, 2 * Ce, x, h.data() );
        // fused [D K][2 Ch]: left and right of each frame side by side
        for(int d = 0; d < D; d++) for(int kk = 0; kk < K; kk++)
        {
          const double* hl = h.data() + size_t(d) * 2 * Ce + size_t(kk) * Ch;
          const double* hr = h.data() + size_t(d) * 2 * Ce + Ce + size_t(kk) * Ch;
          double* fu = fus.data() + size_t(d*K + kk) * 2 * Ch;
          for(int c = 0; c < Ch; c++) { fu[c] = hl[c]; fu[Ch + c] = hr[c]; }
        }
        for(int c = 0; c < Ch; c++) { pair[c] = h[c]; pair[Ch + c] = h[Ce + c]; }
        dpa3_gemm( false, false, DK, M.rhid, 2 * Ch, 1.0, fus.data(), 2 * Ch, M.rl.data(), M.rhid, 0.0, ll.data(), M.rhid );
        dpa3_gemm( false, false, DK, M.rhid, 2 * Ch, 1.0, fus.data(), 2 * Ch, M.rr.data(), M.rhid, 0.0, rr.data(), M.rhid );
        std::fill( s.begin(), s.end(), 0.0 );
        for(int t = 0; t < DK; t++) { const double w = M.sweight[t]; for(int c = 0; c < M.rhid; c++) s[c] += w * ll[size_t(t)*M.rhid + c] * rr[size_t(t)*M.rhid + c]; }
        dpa3_gemm( false, false, 1, Ch, M.rhid, 1.0, s.data(), M.rhid, M.ro.data(), Ch, 0.0, co.data(), Ch );
        dpa4_scalar( M.rgate.data(), Ch, 1, co.data(), pair.data(), z.data(), zg.data() );
        dpa3_gemm( false, false, 1, C, Ch, 1.0, z.data(), Ch, M.rlin2.data(), C, 0.0, y.data(), C );
        for(int c = 0; c < C; c++) desc[c] = x[c] + y[c];
        // GLU network
        const int H = M.fit_h;
        for(int a = 0; a < 2 * H; a++) yh[a] = M.fb[a];
        dpa3_gemm( false, false, 1, 2 * H, C, 1.0, desc.data(), C, M.fw.data(), 2 * H, 1.0, yh.data(), 2 * H );
        double e = M.fout_b;
        for(int a = 0; a < H; a++) { const double sg = dpa4_sig( yh[H + a] ); v[a] = yh[a] * yh[H + a] * sg; e += v[a] * M.fout[a]; }
        B.eatom[i] = e + M.ebias[ti];
        // backward: dE/dv = out, to yh, to desc
        for(int a = 0; a < H; a++)
        {
          const double sg = dpa4_sig( yh[H + a] );
          dy[a] = M.fout[a] * yh[H + a] * sg;
          dy[H + a] = M.fout[a] * yh[a] * ( sg + yh[H + a] * sg * ( 1.0 - sg ) );
        }
        dpa3_gemm( false, true, 1, C, 2 * H, 1.0, dy.data(), 2 * H, M.fw.data(), 2 * H, 0.0, dd.data(), C );
        // read-out backward
        double* gi = BB.G[NB].data() + size_t(i) * DC;
        for(int c = 0; c < C; c++) gi[c] = dd[c];
        dpa3_gemm( false, true, 1, Ch, C, 1.0, dd.data(), C, M.rlin2.data(), C, 0.0, dz.data(), Ch );
        dpa4_scalar_t( M.rgate.data(), Ch, 1, co.data(), pair.data(), zg.data(), dz.data(), dco.data(), dp.data(), dgt.data() );
        dpa3_gemm( false, true, 1, M.rhid, Ch, 1.0, dco.data(), Ch, M.ro.data(), Ch, 0.0, ds.data(), M.rhid );
        for(int t = 0; t < DK; t++)
        {
          const double w = M.sweight[t];
          for(int c = 0; c < M.rhid; c++) { dl[size_t(t)*M.rhid + c] = w * ds[c] * rr[size_t(t)*M.rhid + c]; dr[size_t(t)*M.rhid + c] = w * ds[c] * ll[size_t(t)*M.rhid + c]; }
        }
        dpa3_gemm( false, true, DK, 2 * Ch, M.rhid, 1.0, dl.data(), M.rhid, M.rl.data(), M.rhid, 0.0, dfus.data(), 2 * Ch );
        dpa3_gemm( false, true, DK, 2 * Ch, M.rhid, 1.0, dr.data(), M.rhid, M.rr.data(), M.rhid, 1.0, dfus.data(), 2 * Ch );
        for(int d = 0; d < D; d++) for(int kk = 0; kk < K; kk++)
        {
          double* hl = dh.data() + size_t(d) * 2 * Ce + size_t(kk) * Ch;
          double* hr = dh.data() + size_t(d) * 2 * Ce + Ce + size_t(kk) * Ch;
          const double* fu = dfus.data() + size_t(d*K + kk) * 2 * Ch;
          for(int c = 0; c < Ch; c++) { hl[c] = fu[c]; hr[c] = fu[Ch + c]; }
        }
        for(int c = 0; c < Ch; c++) { dh[c] += dp[c]; dh[Ce + c] += dp[Ch + c]; }
        dpa4_so3_t( M.lmax, M.rlin1.data(), C, 2 * Ce, dh.data(), gi, 1.0 );
      }
    }

    mark(2);
    // reverse blocks
    if( B.locks.size() != n_nodes ) B.locks.assign( n_nodes, 0 );
    BB.dX1.resize( n_nodes * DC );
    for(int b = NB - 1; b >= 0; b--)
    {
      node_prep( b );
      BB.GXp.assign( n_nodes * DW, 0.0 );
      BB.GKn.assign( n_nodes * W, 0.0 );
#     pragma omp parallel
      {
        static thread_local DPA4Scratch S;
        std::vector<double> x2( DC );
#       pragma omp for schedule(dynamic, 1)
        for(int q = 0; q < inum; q++)
        {
          const int i = B.ilist[q];
          DPA4AtomState& A = M.store_forward ? BB.store[b][q] : S.own;
          if( ! M.store_forward ) dpa4_block_forward( M, M.blocks[b], BB, i, offset[i], BB.ecount[i], nbh, BB.X[b].data() + size_t(i) * DC, A, S, x2.data() );
          dpa4_block_vjp( M, M.blocks[b], BB, i, offset[i], BB.ecount[i], nbh, BB.G[b+1].data() + size_t(i) * DC, A, S, BB.dX1.data() + size_t(i) * DC );
        }
      }
      // dE/dK of all nodes -> dE/dxp[0], ghost rows folded to their owners, then through the pre mix
      const DPA4Block& P = M.blocks[b];
#     pragma omp parallel
      {
        std::vector<double> d(C);
#       pragma omp for schedule(static)
        for(size_t n = 0; n < n_nodes; n++)
        {
          for(int f = 0; f < F; f++)
          {
            dpa3_gemm( false, true, 1, C, C, 1.0, BB.GKn.data() + n * W + size_t(f) * C, C, P.foc[f].k.data(), C, 0.0, d.data(), C );
            dpa4_rms_t( P.foc[f].qk_norm, C, BB.Xp.data() + n * DW + size_t(f) * C, d.data(), BB.GXp.data() + n * DW + size_t(f) * C );
          }
        }
      }
      exchange( BB.GXp, false, DW );
      BB.G[b].assign( n_nodes * DC, 0.0 );
#     pragma omp parallel
      {
      std::vector<double> dn(DC), yn(DC), wn(C + 1);
#     pragma omp for schedule(static)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q];
        double* g = BB.G[b].data() + size_t(i) * DC;
        dpa4_so3_t( M.lmax, P.pre_mix.data(), C, W, BB.GXp.data() + size_t(i) * DW, dn.data() );
        const double* xi = BB.X[b].data() + size_t(i) * DC;
        dpa4_eqrms( P.pre_norm, M, xi, yn.data(), wn.data() );
        dpa4_eqrms_t( P.pre_norm, M, xi, wn.data(), dn.data(), g );
        const double* d = BB.dX1.data() + size_t(i) * DC;
        for(size_t t = 0; t < DC; t++) g[t] += d[t];
      }
      }
    }

    mark(3);
    // initial features backward, geometric chain
#   pragma omp parallel
    {
      std::vector<double> ft( 2 * C ), dfl( 2 * C ), Dg( size_t(Em) * Ax ), Ag( size_t(4) * Em ), Fi( 2 * C ), film( 2 * C ), filmn( 2 * C ), Dd( size_t(Em) * Ax );
      std::vector<double> f1z, f1b, f2z, ginb, ghz, ghb, giz, gib, dgg, dgh, dgi, dg2, dr1, h, sh, dh, dhh;
#     pragma omp for schedule(dynamic, 8)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q];
        const int64_t b0 = offset[i];
        const int n = BB.ecount[i];
        const double isd = BB.isd[i];
        const double* G0 = BB.G[0].data() + size_t(i) * DC;
        const double* te = M.te.data() + size_t( B.type[i] ) * C;
        const double* ga = BB.seed.data() + size_t(i) * Em * 4;
        // film chain, forward values
        dpa3_gemm( true, false, Em, Ax, 4, 1.0, ga, Em, ga, Em, 0.0, Dd.data(), Ax );
        dpa3_gemm( false, false, 1, 2 * C, Em * Ax, 1.0, Dd.data(), Em * Ax, M.s_out.data(), 2 * C, 0.0, film.data(), 2 * C );
        if( M.s_a[4] ) dpa4_silu_n( film.data(), 2 * C );
        dpa4_rms( M.film_s, C, film.data(), filmn.data() );
        dpa4_rms( M.film_f, C, film.data() + C, filmn.data() + C );
        for(int c = 0; c < C; c++)
        {
          const double ts = std::tanh( filmn[c] ), tf = std::tanh( filmn[C + c] );
          ft[c] = G0[c] * te[c] * M.film_scale * ( 1.0 - ts * ts );
          ft[C + c] = G0[c] * M.film_shift * ( 1.0 - tf * tf );
        }
        std::fill( dfl.begin(), dfl.end(), 0.0 );
        dpa4_rms_t( M.film_s, C, film.data(), ft.data(), dfl.data() );
        dpa4_rms_t( M.film_f, C, film.data() + C, ft.data() + C, dfl.data() + C );
        dpa3_gemm( false, true, 1, Em * Ax, 2 * C, 1.0, dfl.data(), 2 * C, M.s_out.data(), 2 * C, 0.0, Dg.data(), Em * Ax );
        // Dd[a][b] = sum_x ga[x][a] ga[x][b]: dga[x][a] = sum_b Dg[a][b] ga[x][b] + (a < Ax) sum_e Dg[e][a] ga[x][e]
        for(int x = 0; x < 4; x++)
        {
          double* ax = Ag.data() + size_t(x) * Em;
          const double* gx = ga + size_t(x) * Em;
          for(int a = 0; a < Em; a++)
          {
            double t = 0.0;
            for(int b = 0; b < Ax; b++) t += Dg[size_t(a)*Ax + b] * gx[b];
            ax[a] = t;
          }
          for(int a = 0; a < Ax; a++)
          {
            double t = 0.0;
            for(int e = 0; e < Em; e++) t += Dg[size_t(e)*Ax + a] * gx[e];
            ax[a] += t;
          }
        }
        // isd gradient: ga = isd sum(raw), gie = isd sum(raw)
        double gsd = 0.0;
        for(size_t t = 0; t < size_t(4) * Em; t++) gsd += Ag[t] * ga[t];         // Ag dga / isd d(isd)  (ga already has isd)
        gsd /= isd;
        const double* gs = BB.gie.data() + size_t(i) * DC;
        for(size_t t = C; t < DC; t++) gsd += G0[t] * gs[t];
        const double ddeg = gsd * ( -0.5 ) * isd * isd * isd;
        const int hr = M.hrf;
        const double* Rb = BB.rbf.data() + size_t(b0) * M.nrad;
        // forward of g, all the edges of the atom
        dpa3_grow( f1z, size_t(n) * M.s_r ); dpa3_grow( f1b, size_t(n) * M.s_r ); dpa3_grow( f2z, size_t(n) * M.s_r ); dpa3_grow( ginb, size_t(n) * Sc );
        dpa3_grow( ghz, size_t(n) * M.s_gh ); dpa3_grow( ghb, size_t(n) * M.s_gh ); dpa3_grow( giz, size_t(n) * Em ); dpa3_grow( gib, size_t(n) * Em );
        dpa3_grow( dgg, size_t(n) * Em ); dpa3_grow( dgh, size_t(n) * M.s_gh ); dpa3_grow( dgi, size_t(n) * Sc ); dpa3_grow( dg2, size_t(n) * M.s_r ); dpa3_grow( dr1, size_t(n) * M.s_r );
        dpa3_gemm( false, false, n, M.s_r, M.nrad, 1.0, Rb, M.nrad, M.s_r1.data(), M.s_r, 0.0, f1z.data(), M.s_r );
        std::copy( f1z.begin(), f1z.begin() + size_t(n) * M.s_r, f1b.begin() );
        if( M.s_a[0] ) dpa4_silu_n( f1b.data(), n * M.s_r );
        dpa3_gemm( false, false, n, M.s_r, M.s_r, 1.0, f1b.data(), M.s_r, M.s_r2.data(), M.s_r, 0.0, f2z.data(), M.s_r );
        std::copy( f2z.begin(), f2z.begin() + size_t(n) * M.s_r, f1b.begin() );
        if( M.s_a[1] ) dpa4_silu_n( f1b.data(), n * M.s_r );
        {
          const double* td = M.tt.data() + size_t( B.type[i] ) * M.tdim;
          for(int k = 0; k < n; k++)
          {
            double* g = ginb.data() + size_t(k) * Sc;
            const double* f2 = f1b.data() + size_t(k) * M.s_r;
            const double* ts = M.tt.data() + size_t( B.type[ nbh[b0 + k] ] ) * M.tdim;
            for(int a = 0; a < M.s_r; a++) g[a] = f2[a];
            for(int a = 0; a < M.tdim; a++) { g[M.s_r + a] = ts[a]; g[M.s_r + M.tdim + a] = td[a]; }
          }
        }
        dpa3_gemm( false, false, n, M.s_gh, Sc, 1.0, ginb.data(), Sc, M.s_g1.data(), M.s_gh, 0.0, ghz.data(), M.s_gh );
        std::copy( ghz.begin(), ghz.begin() + size_t(n) * M.s_gh, ghb.begin() );
        if( M.s_a[2] ) dpa4_silu_n( ghb.data(), n * M.s_gh );
        dpa3_gemm( false, false, n, Em, M.s_gh, 1.0, ghb.data(), M.s_gh, M.s_g2.data(), Em, 0.0, giz.data(), Em );
        std::copy( giz.begin(), giz.begin() + size_t(n) * Em, gib.begin() );
        if( M.s_a[3] ) dpa4_silu_n( gib.data(), n * Em );
        for(int k = 0; k < n; k++)
        {
          const int64_t e = b0 + k;
          const double* v = rij + 3 * e;
          const double r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
          const double inv = 1.0 / std::sqrt( r2 + M.seed_eps * M.seed_eps );
          const double en = BB.envv[e];
          BB.Genv[e] += ddeg * 2.0 * en;
          const double s = en * inv;
          const double rt[4] = { s, s * v[0] * inv, s * v[1] * inv, s * v[2] * inv };
          const double* gi = gib.data() + size_t(k) * Em;
          // dE/drt, dE/dg (outer product rt g isd)
          double rtg[4] = { 0.0, 0.0, 0.0, 0.0 };
          double* dg = dgg.data() + size_t(k) * Em;
          std::fill( dg, dg + Em, 0.0 );
          for(int x = 0; x < 4; x++)
          {
            const double* ax = Ag.data() + size_t(x) * Em;
            for(int a = 0; a < Em; a++) { rtg[x] += ax[a] * gi[a] * isd; dg[a] += ax[a] * rt[x] * isd; }
          }
          // rt: env part into Genv, vector part into Gv
          BB.Genv[e] += rtg[0] * inv + inv * inv * ( rtg[1]*v[0] + rtg[2]*v[1] + rtg[3]*v[2] );
          const double i3 = en * inv * inv * inv;
          const double i4 = en * inv * inv;
          double* gv = BB.Gv.data() + 3 * e;
          const double rv = rtg[1]*v[0] + rtg[2]*v[1] + rtg[3]*v[2];
          for(int x = 0; x < 3; x++) gv[x] += -rtg[0] * i3 * v[x] + rtg[1 + x] * i4 - 2.0 * rv * v[x] * i4 * inv * inv;
          // geometric initial embedding
          double* gd = BB.GD.data() + size_t(e) * M.nDq;
          double* gr = BB.Grf.data() + size_t(e) * L * C;
          const double* rfe = BB.rf.data() + size_t(e) * L * C;
          const double* Dqe = BB.Dq.data() + size_t(e) * M.nDq;
          for(int t = 0; t < D - 1; t++)
          {
            const double z = Dqe[ M.z_col[t] ];
            const double* gg = G0 + size_t( t + 1 ) * C;
            const double* rv2 = rfe + size_t( M.z_slot[t] ) * C;
            double* gs2 = gr + size_t( M.z_slot[t] ) * C;
            double sum = 0.0;
            for(int c = 0; c < C; c++) { gs2[c] += gg[c] * z * isd; sum += gg[c] * rv2[c]; }
            gd[ M.z_col[t] ] += sum * isd;
          }
          // env derivative: rf / env
          double gen = 0.0;
          for(int t = 0; t < L * C; t++) gen += gr[t] * ( en != 0.0 ? rfe[t] / en : 0.0 );
          BB.Genv[e] += gen;
        }
        // g layers backward, all the edges of the atom
        double* Grb = BB.Grbf.data() + size_t(b0) * M.nrad;
        if( M.s_a[3] ) dpa4_dsilu_n( giz.data(), n * Em, dgg.data() );
        dpa3_gemm( false, true, n, M.s_gh, Em, 1.0, dgg.data(), Em, M.s_g2.data(), Em, 0.0, dgh.data(), M.s_gh );
        if( M.s_a[2] ) dpa4_dsilu_n( ghz.data(), n * M.s_gh, dgh.data() );
        dpa3_gemm( false, true, n, Sc, M.s_gh, 1.0, dgh.data(), M.s_gh, M.s_g1.data(), M.s_gh, 0.0, dgi.data(), Sc );
        for(int k = 0; k < n; k++) for(int a = 0; a < M.s_r; a++) dg2[size_t(k)*M.s_r + a] = dgi[size_t(k)*Sc + a];
        if( M.s_a[1] ) dpa4_dsilu_n( f2z.data(), n * M.s_r, dg2.data() );
        dpa3_gemm( false, true, n, M.s_r, M.s_r, 1.0, dg2.data(), M.s_r, M.s_r2.data(), M.s_r, 0.0, dr1.data(), M.s_r );
        if( M.s_a[0] ) dpa4_dsilu_n( f1z.data(), n * M.s_r, dr1.data() );
        dpa3_gemm( false, true, n, M.nrad, M.s_r, 1.0, dr1.data(), M.s_r, M.s_r1.data(), M.s_r, 1.0, Grb, M.nrad );
        // radial features = env act( rms( rbf W1 ) ) W2: dE/dfeatures
        dpa3_grow( h, size_t(n) * hr ); dpa3_grow( sh, size_t(n) * hr ); dpa3_grow( dh, size_t(n) * hr ); dpa3_grow( dhh, size_t(n) * hr );
        dpa3_gemm( false, false, n, hr, M.nrad, 1.0, Rb, M.nrad, M.rf1.data(), hr, 0.0, h.data(), hr );
        for(int k = 0; k < n; k++) dpa4_rms( M.rf_norm, hr, h.data() + size_t(k) * hr, sh.data() + size_t(k) * hr );
        dpa3_gemm( false, true, n, hr, L * C, 1.0, BB.Grf.data() + size_t(b0) * L * C, L * C, M.rf2.data(), L * C, 0.0, dh.data(), hr );
        for(int k = 0; k < n; k++) { const double en = BB.envv[b0 + k]; double* r = dh.data() + size_t(k) * hr; for(int a = 0; a < hr; a++) r[a] *= en; }
        dpa4_dsilu_n( sh.data(), n * hr, dh.data() );
        std::fill( dhh.begin(), dhh.begin() + size_t(n) * hr, 0.0 );
        for(int k = 0; k < n; k++) dpa4_rms_t( M.rf_norm, hr, h.data() + size_t(k) * hr, dh.data() + size_t(k) * hr, dhh.data() + size_t(k) * hr );
        dpa3_gemm( false, true, n, M.nrad, hr, 1.0, dhh.data(), hr, M.rf1.data(), hr, 1.0, Grb, M.nrad );
        for(int k = 0; k < n; k++)
        {
          const int64_t e = b0 + k;
          const double* v = rij + 3 * e;
          double* gd = BB.GD.data() + size_t(e) * M.nDq;
          double* gv = BB.Gv.data() + 3 * e;
          // radial: dE/dr from env and rbf, Wigner chain
          const double r = BB.rr[e];
          double e5, d5, e7, d7;
          dpa4_envelope( M.env_c, r, M.rc, e5, d5 );
          dpa4_envelope( M.rbf_c, r, M.rc, e7, d7 );
          double gr0 = BB.Genv[e] * d5;
          const double* grb = BB.Grbf.data() + size_t(e) * M.nrad;
          if( M.gauss != 0.0 ) for(int a = 0; a < M.nrad; a++)
          {
            const double dr = r - M.freq[a], ra = std::exp( M.gauss * dr * dr );
            gr0 += grb[a] * ra * ( d7 + 2.0 * M.gauss * dr * e7 );
          }
          else for(int a = 0; a < M.nrad; a++)
          {
            const double f = M.freq[a], x = r * f, sx = std::sin(x), cx = std::cos(x);
            gr0 += grb[a] * ( f * sx / x * d7 + f * f * ( cx / x - sx / ( x * x ) ) * e7 );
          }
          double qq[4], J[12], gq[4];
          dpa4_quaternion( v, r, M.eps, M.weps, qq, J );
          dpa4_wigner_t( M, qq, gd, gq );
          for(int x = 0; x < 3; x++)
          {
            double t = gr0 * v[x] / r;
            for(int a = 0; a < 4; a++) t += gq[a] * J[a*3 + x];
            gv[x] += t;
          }
        }
      }
    }

    mark(4);
    // forces, per-atom virials, energies into the grid
    B.f.assign( 3 * n_nodes, 0.0 );
    B.vatom.assign( 9 * n_nodes, 0.0 );
#   pragma omp parallel for schedule(dynamic, 16)
    for(int q = 0; q < inum; q++)
    {
      const int i = B.ilist[q];
      double fi[3] = { 0, 0, 0 }, vi[9] = { 0,0,0, 0,0,0, 0,0,0 };
      const int n = BB.ecount[i];
      for(int64_t b = offset[i]; b < offset[i] + n; b++)
      {
        const double* g = BB.Gv.data() + 3 * b;
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
    mark(5);
  }

}
