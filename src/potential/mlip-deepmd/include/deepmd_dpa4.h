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

// DPA4 (SeZM) descriptor, host only: SO(3) node features x[D][C] (D = (lmax+1)^2), l <= 2, m <= 1, one focus, one
// attention head, as the presets of deepmd (exaStamp/scripts/python/deepmd_dpa4.py is the numpy specification).
//
//   per edge j -> i : quaternion q(d_ij), Wigner blocks D(q), envelope, Bessel radial basis, radial features
//   x0 = FiLM(type embedding; env seed) on l = 0, geometric initial embedding on l > 0
//   block: xp = pre_mix(x);  edge message: rotate xp_j to the edge frame, radial mixer, SO(2) layers, rotate back
//          out = gate_i sum_j alpha_ij msg_ij;  out += grid(out, xp_i);  x1 = x + norm(post_mix(out));
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

  struct DPA4Block
  {
    std::vector<double> pre_mix, post_mix;                  // [L][C][C]
    std::vector< std::vector<double> > lin0, linm, gate;    // [(L C)][(L C)], [lmax C][2 lmax C], [C][lmax C] (no gate in the last layer)
    std::vector<double> mixw, mixb;                         // [L C][nk], [C]; empty: no mixer
    DPA4Rms qk_norm, gate_norm;
    std::vector<double> q, k, logit_w, gate_w;              // [C][C], [C][C], [C], [C]
    double log_null = 0.0;
    std::vector<double> gexp, gcon, ggate, gres;            // [L][C][K C], [L][K C][C], [2C][C], [C]
    DPA4EqRms post_norm, ffn_norm;
    std::vector<double> flin1, flin2, fl, fr, fo, fgate;    // [L][C][2 Ce], [L][Ce][C], [Ch][Ch] x3, [2 Ch][Ch]
  };

  struct DPA4Model
  {
    int nt = 0, C = 0, lmax = 0, L = 0, D = 0, Dm = 0, nnei = 0, K = 0, Ch = 0, Ce = 0, Ng = 0, nrad = 0, nk = 0;
    int emb = 0, axis = 0, tdim = 0, hrf = 0, fit_h = 0;
    double rc = 0.0, eps = 0.0, weps = 0.0, deg_floor = 0.0, film_scale = 0.0, film_shift = 0.0, seed_eps = 0.0;
    std::vector<double> env_c, rbf_c, freq, te, tt;                 // series coefficients, frequencies, type embeddings
    std::vector<double> rf1, rf2; DPA4Rms rf_norm;
    std::vector<double> s_r1, s_r2, s_g1, s_g2, s_out;              // env seed layers, s_* activation flags below
    bool s_a[5] = { true, true, true, true, false };
    int s_r = 0, s_g = 0, s_gh = 0;
    DPA4Rms film_s, film_f;
    std::vector<double> C_l2;                                       // [25][35] degree 4 monomials
    std::vector<int> mono;                                          // [35][4]
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

  struct DPA4Buffers
  {
    DPA3Buffers c;
    std::vector<int> ecount;                                        // kept neighbors per atom row
    std::vector<double> envv, rr, rbf, Dq, rf, GD, Grf, Genv, Grbf, Gv, isd, deg, gie, seed;
    std::vector<double> Xp, Kn, GXp, GKn, dX1;
    std::vector< std::vector<double> > X, G;
  };

  // ------------------------------------------------------------------ small dense helpers

  inline double dpa4_sig(double x) { return 1.0 / ( 1.0 + std::exp(-x) ); }

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
    double s = 0.0;
#   pragma omp simd reduction(+:s)
    for(int c = 0; c < C; c++) s += x[c] * x[c];
    const double inv = 1.0 / std::sqrt( s / C + N.eps );
    for(int c = 0; c < C; c++) y[c] = x[c] * inv * N.scale[c];
  }
  inline void dpa4_rms_t(const DPA4Rms& N, int C, const double* x, const double* dy, double* dx)
  {
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

  // packed Wigner blocks of one edge: D[0] = 1, then l = 1 (9), l = 2 (25)
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

  struct DPA4Scratch
  {
    std::vector<double> X0, rad, kern, zg, msg, wedge, alpha, ex, logit, q0, l0, sum, out, gl, gr, cin, vv, zgate, pair, res, out2, mixed;
    std::vector<double> so, wn1, x1, y, wn2, h, fgl, fgr, fc, fco, fv, fzgate, fpair;
    std::vector< std::vector<double> > XL, Y, GG, Sd;
    std::vector<double> tmp, tmp2, tmp3, tmp4, dmsg, dX, dY, dG, dz, dalpha, ddso, dout, dc, dp, dleft, dright, dh, dpairv, df, dxp, drad, dres;
    int n = 0;
  };

  // forward of block b for atom i; Xp, Kn: nodes (ghosts included); x, y: node row in / out
  inline void dpa4_block_forward(const DPA4Model& M, const DPA4Block& P, const DPA4Buffers& BB, int i, const int64_t e0, const int n,
                                 const int* nbh, const double* x, DPA4Scratch& S, double* x2)
  {
    const int C = M.C, D = M.D, Dm = M.Dm, L = M.L, DK = D * M.K, K = M.K, Ch = M.Ch, Ce = M.Ce;
    const size_t DC = size_t(D) * C, DmC = size_t(Dm) * C;
    const double* Xp = BB.Xp.data();
    const double* xp = Xp + size_t(i) * DC;
    S.n = n;
    const int nl = P.lin0.size();
    if( int(S.XL.size()) < nl + 1 ) { S.XL.resize( nl + 1 ); S.Y.resize( nl ); S.GG.resize( nl ); S.Sd.resize( nl ); }
    for(auto& v : S.XL) dpa3_grow( v, size_t(n) * DmC );
    dpa3_grow( S.X0, size_t(n) * DmC ); dpa3_grow( S.rad, size_t(n) * L * C ); dpa3_grow( S.msg, size_t(n) * DC );
    // radial features of the edges, types fused
    const int ty = BB.c.type[i];
    const double* tei = M.te.data() + size_t(ty) * C;
    for(int e = 0; e < n; e++)
    {
      const double* rf = BB.rf.data() + size_t(e0 + e) * L * C;
      const double* tej = M.te.data() + size_t( BB.c.type[ nbh[e0 + e] ] ) * C;
      double* r = S.rad.data() + size_t(e) * L * C;
      for(int l = 0; l < L; l++) for(int c = 0; c < C; c++) r[l*C + c] = rf[l*C + c] + tei[c] + tej[c];
    }
    // rotate the source features to the edge frame
    std::fill( S.X0.begin(), S.X0.begin() + size_t(n) * DmC, 0.0 );
    const int nrot = M.rot_r.size();
    for(int e = 0; e < n; e++)
    {
      const double* xs = Xp + size_t( nbh[e0 + e] ) * DC;
      const double* Dq = BB.Dq.data() + size_t(e0 + e) * 35;
      double* xl = S.X0.data() + size_t(e) * DmC;
      for(int t = 0; t < nrot; t++)
      {
        const double f = Dq[ M.rot_o[t] ];
        const double* sc = xs + size_t( M.rot_d[t] ) * C;
        double* dc = xl + size_t( M.rot_r[t] ) * C;
#       pragma omp simd
        for(int c = 0; c < C; c++) dc[c] += f * sc[c];
      }
    }
    // radial mixer
    double* X = S.XL[0].data();
    if( M.nk > 0 )
    {
      dpa3_grow( S.kern, size_t(n) * M.nk );
      dpa3_gemm( false, false, n, M.nk, L * C, 1.0, S.rad.data(), L * C, P.mixw.data(), M.nk, 0.0, S.kern.data(), M.nk );
      std::fill( X, X + size_t(n) * DmC, 0.0 );
      const int nm = M.mix_r.size();
      for(int e = 0; e < n; e++)
      {
        const double* k = S.kern.data() + size_t(e) * M.nk;
        const double* x0 = S.X0.data() + size_t(e) * DmC;
        double* xo = X + size_t(e) * DmC;
        for(int t = 0; t < nm; t++)
        {
          const double f = k[ M.mix_c[t] ];
          const double* sc = x0 + size_t( M.mix_s[t] ) * C;
          double* dc = xo + size_t( M.mix_r[t] ) * C;
#         pragma omp simd
          for(int c = 0; c < C; c++) dc[c] += f * sc[c] * P.mixb[c];
        }
      }
    }
    else
    {
      for(int e = 0; e < n; e++)
      {
        const double* x0 = S.X0.data() + size_t(e) * DmC;
        const double* rd = S.rad.data() + size_t(e) * L * C;
        double* xo = X + size_t(e) * DmC;
        for(int r = 0; r < Dm; r++) { const double* rr = rd + size_t( M.deg_m[r] ) * C; for(int c = 0; c < C; c++) xo[size_t(r)*C + c] = x0[size_t(r)*C + c] * rr[c]; }
      }
    }
    // SO(2) layers
    const int nm0 = L * C, nmm = M.lmax * C;
    for(int k = 0; k < nl; k++)
    {
      const double* xk = S.XL[k].data();
      dpa3_grow( S.Y[k], size_t(n) * DmC );
      double* y = S.Y[k].data();
      dpa3_gemm( false, false, n, nm0, nm0, 1.0, xk, int(DmC), P.lin0[k].data(), nm0, 0.0, y, int(DmC) );
      if( nmm > 0 )
      {
        const double* wm = P.linm[k].data();     // [nmm][2 nmm]: u, v
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
        dpa3_grow( S.GG[k], size_t(n) * nmm ); dpa3_grow( S.Sd[k], size_t(n) * C );
        const int gs = M.lmax;
        if( nmm > 0 ) dpa3_gemm( false, false, n, nmm, C, 1.0, y, int(DmC), P.gate[k].data(), nmm, 0.0, S.GG[k].data(), nmm );
        for(size_t t = 0; t < size_t(n) * nmm; t++) S.GG[k][t] = dpa4_sig( S.GG[k][t] );
        for(int e = 0; e < n; e++)
        {
          const double* ye = y + size_t(e) * DmC;
          const double* xe = xk + size_t(e) * DmC;
          double* oe = xo + size_t(e) * DmC;
          double* sd = S.Sd[k].data() + size_t(e) * C;
          const double* ge = S.GG[k].data() + size_t(e) * nmm;
          for(int c = 0; c < C; c++)
          {
            const double s = dpa4_sig( ye[c] );
            oe[c] = xe[c] + ye[c] * s;
            sd[c] = s + ye[c] * s * ( 1.0 - s );
          }
          for(int r = 1; r < Dm; r++)
          {
            const double* gg = ge + size_t( M.sign_g[r] ) * C;
            const double* yr = ye + size_t(r) * C;
            const double* xr = xe + size_t(r) * C;
            double* orr = oe + size_t(r) * C;
#           pragma omp simd
            for(int c = 0; c < C; c++) orr[c] = xr[c] + yr[c] * gg[c];
          }
          (void) gs;
        }
      }
      else for(size_t t = 0; t < size_t(n) * DmC; t++) xo[t] = xk[t] + y[t];
    }
    // rotate back
    std::fill( S.msg.begin(), S.msg.begin() + size_t(n) * DC, 0.0 );
    {
      const double* xf = S.XL[nl].data();
      for(int e = 0; e < n; e++)
      {
        const double* Dq = BB.Dq.data() + size_t(e0 + e) * 35;
        const double* xl = xf + size_t(e) * DmC;
        double* m = S.msg.data() + size_t(e) * DC;
        for(int t = 0; t < nrot; t++)
        {
          const int d = M.rot_d[t];
          const double f = Dq[ M.rot_o[t] ] * M.rescale[d];
          const double* sc = xl + size_t( M.rot_r[t] ) * C;
          double* dc = m + size_t(d) * C;
#         pragma omp simd
          for(int c = 0; c < C; c++) dc[c] += f * sc[c];
        }
      }
    }
    // attention
    dpa3_grow( S.l0, C ); dpa3_grow( S.q0, C ); dpa3_grow( S.logit, n ); dpa3_grow( S.alpha, n ); dpa3_grow( S.ex, n );
    dpa4_rms( P.qk_norm, C, xp, S.l0.data() );
    dpa3_gemm( false, false, 1, C, C, 1.0, S.l0.data(), C, P.q.data(), C, 0.0, S.q0.data(), C );
    const double ish = 1.0 / std::sqrt( double(C) );
    double mx = P.log_null;
    for(int e = 0; e < n; e++)
    {
      const double* kj = BB.Kn.data() + size_t( nbh[e0 + e] ) * C;
      const double* r0 = S.rad.data() + size_t(e) * L * C;
      double s = 0.0, b = 0.0;
#     pragma omp simd reduction(+:s,b)
      for(int c = 0; c < C; c++) { s += S.q0[c] * kj[c]; b += r0[c] * P.logit_w[c]; }
      const double env = BB.envv[e0 + e];
      const double l = s * ish + b + 2.0 * std::log( env );
      S.logit[e] = l;
      mx = std::max( mx, l );
    }
    double den = std::exp( P.log_null - mx );
    for(int e = 0; e < n; e++) { S.ex[e] = std::exp( S.logit[e] - mx ); den += S.ex[e]; }
    for(int e = 0; e < n; e++) S.alpha[e] = S.ex[e] / den;
    dpa3_grow( S.sum, DC ); dpa3_grow( S.out, DC );
    std::fill( S.sum.begin(), S.sum.begin() + DC, 0.0 );
    for(int e = 0; e < n; e++)
    {
      const double a = S.alpha[e];
      const double* m = S.msg.data() + size_t(e) * DC;
      double* s = S.sum.data();
#     pragma omp simd
      for(size_t t = 0; t < DC; t++) s[t] += a * m[t];
    }
    dpa3_grow( S.tmp, C );
    dpa4_rms( P.gate_norm, C, xp, S.tmp.data() );
    {
      double z = 0.0;
      for(int c = 0; c < C; c++) z += S.tmp[c] * P.gate_w[c];
      const double g = dpa4_sig( z );
      S.so.resize(1); S.so[0] = g;
      for(size_t t = 0; t < DC; t++) S.out[t] = S.sum[t] * g;
    }
    // message-node grid product
    const int C3 = K * C;
    const size_t DC3 = size_t(D) * C3;
    dpa3_grow( S.h, 2 * DC3 ); dpa3_grow( S.gl, size_t(M.Ng) * C ); dpa3_grow( S.gr, size_t(M.Ng) * C ); dpa3_grow( S.cin, DC3 ); dpa3_grow( S.vv, DC3 );
    dpa3_grow( S.zgate, C ); dpa3_grow( S.pair, 2 * C ); dpa3_grow( S.res, DC ); dpa3_grow( S.out2, DC );
    dpa4_so3( M.lmax, P.gexp.data(), C, C3, S.out.data(), S.h.data() );
    dpa4_so3( M.lmax, P.gexp.data(), C, C3, xp, S.h.data() + DC3 );
    dpa3_gemm( false, false, M.Ng, C, DK, 1.0, M.T.data(), DK, S.h.data(), C, 0.0, S.gl.data(), C );
    dpa3_gemm( false, false, M.Ng, C, DK, 1.0, M.T.data(), DK, S.h.data() + DC3, C, 0.0, S.gr.data(), C );
    dpa3_grow( S.tmp2, size_t(M.Ng) * C );
    for(size_t t = 0; t < size_t(M.Ng) * C; t++) S.tmp2[t] = S.gl[t] * S.gr[t];
    dpa3_gemm( false, false, DK, C, M.Ng, 1.0, M.F.data(), M.Ng, S.tmp2.data(), C, 0.0, S.cin.data(), C );
    for(int c = 0; c < C; c++) { S.pair[c] = S.out[c]; S.pair[C + c] = xp[c]; }
    dpa4_scalar( P.ggate.data(), C, DK, S.cin.data(), S.pair.data(), S.vv.data(), S.zgate.data() );
    dpa4_so3( M.lmax, P.gcon.data(), C3, C, S.vv.data(), S.res.data() );
    for(int d = 0; d < D; d++) for(int c = 0; c < C; c++) S.out2[size_t(d)*C + c] = S.out[size_t(d)*C + c] + S.res[size_t(d)*C + c] * P.gres[c];
    // final mix, norm, residual
    dpa3_grow( S.mixed, DC ); dpa3_grow( S.wn1, C + 1 ); dpa3_grow( S.x1, DC ); dpa3_grow( S.y, DC ); dpa3_grow( S.wn2, C + 1 );
    dpa4_so3( M.lmax, P.post_mix.data(), C, C, S.out2.data(), S.mixed.data() );
    dpa4_eqrms( P.post_norm, M, S.mixed.data(), S.y.data(), S.wn1.data() );
    for(size_t t = 0; t < DC; t++) S.x1[t] = x[t] + S.y[t];
    // FFN
    dpa4_eqrms( P.ffn_norm, M, S.x1.data(), S.y.data(), S.wn2.data() );
    const size_t DCe = size_t(D) * Ce;
    dpa3_grow( S.fgl, size_t(M.Ng) * Ch ); dpa3_grow( S.fgr, size_t(M.Ng) * Ch ); dpa3_grow( S.fc, DCe ); dpa3_grow( S.fco, DCe ); dpa3_grow( S.fv, DCe );
    dpa3_grow( S.fzgate, Ch ); dpa3_grow( S.fpair, 2 * Ch ); dpa3_grow( S.wedge, 2 * DCe ); dpa3_grow( S.tmp3, 2 * DCe ); dpa3_grow( S.tmp4, size_t(M.Ng) * Ch );
    dpa4_so3( M.lmax, P.flin1.data(), C, 2 * Ce, S.y.data(), S.wedge.data() );
    // left, right: rows [D K][Ch], split in the channels of a row of 2 Ce
    for(int d = 0; d < D; d++)
    {
      const double* h = S.wedge.data() + size_t(d) * 2 * Ce;
      for(int t = 0; t < Ce; t++) { S.tmp3[size_t(d)*Ce + t] = h[t]; S.tmp3[DCe + size_t(d)*Ce + t] = h[Ce + t]; }
    }
    for(int c = 0; c < Ch; c++) { S.fpair[c] = S.tmp3[c]; S.fpair[Ch + c] = S.tmp3[DCe + c]; }
    dpa3_gemm( false, false, DK, Ch, Ch, 1.0, S.tmp3.data(), Ch, P.fl.data(), Ch, 0.0, S.fc.data(), Ch );
    dpa3_gemm( false, false, DK, Ch, Ch, 1.0, S.tmp3.data() + DCe, Ch, P.fr.data(), Ch, 0.0, S.fco.data(), Ch );
    dpa3_gemm( false, false, M.Ng, Ch, DK, 1.0, M.T.data(), DK, S.fc.data(), Ch, 0.0, S.fgl.data(), Ch );
    dpa3_gemm( false, false, M.Ng, Ch, DK, 1.0, M.T.data(), DK, S.fco.data(), Ch, 0.0, S.fgr.data(), Ch );
    for(size_t t = 0; t < size_t(M.Ng) * Ch; t++) S.tmp4[t] = S.fgl[t] * S.fgr[t];
    dpa3_gemm( false, false, DK, Ch, M.Ng, 1.0, M.F.data(), M.Ng, S.tmp4.data(), Ch, 0.0, S.fc.data(), Ch );
    dpa3_gemm( false, false, DK, Ch, Ch, 1.0, S.fc.data(), Ch, P.fo.data(), Ch, 0.0, S.fco.data(), Ch );
    dpa4_scalar( P.fgate.data(), Ch, DK, S.fco.data(), S.fpair.data(), S.fv.data(), S.fzgate.data() );
    dpa4_so3( M.lmax, P.flin2.data(), Ce, C, S.fv.data(), x2 );
    for(size_t t = 0; t < DC; t++) x2[t] += S.x1[t];
  }

  // reverse of dpa4_block_forward: Gout[D][C] dE/dx2 -> dx1 (dE/dx, node term), source atoms and edge accumulators
  inline void dpa4_block_vjp(const DPA4Model& M, const DPA4Block& P, DPA4Buffers& BB, int i, const int64_t e0, const int n,
                             const int* nbh, const double* Gout, DPA4Scratch& S, double* dxnode)
  {
    const int C = M.C, D = M.D, Dm = M.Dm, L = M.L, DK = D * M.K, K = M.K, Ch = M.Ch, Ce = M.Ce;
    const size_t DC = size_t(D) * C, DmC = size_t(Dm) * C, DCe = size_t(D) * Ce;
    const double* Xp = BB.Xp.data();
    const double* xp = Xp + size_t(i) * DC;
    DPA3Buffers& B = BB.c;
    dpa3_grow( S.dX, DC ); dpa3_grow( S.dxp, DC );
    std::fill( S.dxp.begin(), S.dxp.begin() + DC, 0.0 );
    // FFN: x2 = x1 + lin2( scalar( o grid( l, r ) ) )
    dpa3_grow( S.dc, 2 * DCe ); dpa3_grow( S.dp, size_t(M.Ng) * Ch ); dpa3_grow( S.dleft, DCe ); dpa3_grow( S.dright, DCe ); dpa3_grow( S.dh, 2 * DCe );
    dpa3_grow( S.dpairv, 2 * Ch ); dpa3_grow( S.df, Ch ); dpa3_grow( S.tmp2, DCe ); dpa3_grow( S.dG, size_t(M.Ng) * Ch );
    dpa4_so3_t( M.lmax, P.flin2.data(), Ce, C, Gout, S.dc.data() );
    dpa4_scalar_t( P.fgate.data(), Ch, DK, S.fco.data(), S.fpair.data(), S.fzgate.data(), S.dc.data(), S.tmp2.data(), S.dpairv.data(), S.df.data() );
    // fco = fc o: dfc = dfco o^T
    dpa3_gemm( false, true, DK, Ch, Ch, 1.0, S.tmp2.data(), Ch, P.fo.data(), Ch, 0.0, S.dright.data(), Ch );
    // fc = F p: dp = F^T dfc
    dpa3_gemm( true, false, M.Ng, Ch, DK, 1.0, M.F.data(), M.Ng, S.dright.data(), Ch, 0.0, S.dp.data(), Ch );
    for(size_t t = 0; t < size_t(M.Ng) * Ch; t++) S.dG[t] = S.dp[t] * S.fgr[t];
    dpa3_gemm( true, false, DK, Ch, M.Ng, 1.0, M.T.data(), DK, S.dG.data(), Ch, 0.0, S.dright.data(), Ch );
    dpa3_gemm( false, true, DK, Ch, Ch, 1.0, S.dright.data(), Ch, P.fl.data(), Ch, 0.0, S.dleft.data(), Ch );
    for(size_t t = 0; t < size_t(M.Ng) * Ch; t++) S.dG[t] = S.dp[t] * S.fgl[t];
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
    dpa4_eqrms_t( P.ffn_norm, M, S.x1.data(), S.wn2.data(), S.dY.data(), S.dX.data() );
    // dE/dx1: node term of dE/dx; through the norm and the final mix to dE/dout2
    for(size_t t = 0; t < DC; t++) dxnode[t] = S.dX[t];
    std::fill( S.dY.begin(), S.dY.begin() + DC, 0.0 );
    dpa4_eqrms_t( P.post_norm, M, S.mixed.data(), S.wn1.data(), S.dX.data(), S.dY.data() );
    dpa3_grow( S.dout, DC ); dpa3_grow( S.dres, DC );
    dpa4_so3_t( M.lmax, P.post_mix.data(), C, C, S.dY.data(), S.dout.data() );
    // grid product: res = scale contract( scalar( F( T L * T R ) ) )
    const int C3 = K * C;
    const size_t DC3 = size_t(D) * C3;
    for(int d = 0; d < D; d++) for(int c = 0; c < C; c++) S.dres[size_t(d)*C + c] = S.dout[size_t(d)*C + c] * P.gres[c];
    dpa3_grow( S.dc, 2 * DC3 ); dpa3_grow( S.dz, DC3 ); dpa3_grow( S.drad, size_t(M.Ng) * C );
    dpa3_grow( S.dpairv, 2 * C ); dpa3_grow( S.df, C );
    dpa4_so3_t( M.lmax, P.gcon.data(), C3, C, S.dres.data(), S.dc.data() );
    dpa4_scalar_t( P.ggate.data(), C, DK, S.cin.data(), S.pair.data(), S.zgate.data(), S.dc.data(), S.dz.data(), S.dpairv.data(), S.df.data() );
    dpa3_gemm( true, false, M.Ng, C, DK, 1.0, M.F.data(), M.Ng, S.dz.data(), C, 0.0, S.drad.data(), C );
    dpa3_grow( S.dG, size_t(M.Ng) * C );
    for(size_t t = 0; t < size_t(M.Ng) * C; t++) S.dG[t] = S.drad[t] * S.gr[t];
    dpa3_gemm( true, false, DK, C, M.Ng, 1.0, M.T.data(), DK, S.dG.data(), C, 0.0, S.dz.data(), C );
    dpa4_so3_t( M.lmax, P.gexp.data(), C, C3, S.dz.data(), S.dout.data(), 1.0 );    // dE/dout
    for(size_t t = 0; t < size_t(M.Ng) * C; t++) S.dG[t] = S.drad[t] * S.gl[t];
    dpa3_gemm( true, false, DK, C, M.Ng, 1.0, M.T.data(), DK, S.dG.data(), C, 0.0, S.dz.data(), C );
    dpa4_so3_t( M.lmax, P.gexp.data(), C, C3, S.dz.data(), S.dxp.data() );          // dE/dxp_i
    for(int c = 0; c < C; c++) { S.dout[c] += S.dpairv[c]; S.dxp[c] += S.dpairv[C + c]; }
    // out = gate sum_e alpha_e msg_e
    const double g = S.so[0];
    double dg = 0.0;
    for(size_t t = 0; t < DC; t++) dg += S.dout[t] * S.sum[t];
    {
      // gate = sig( rms(xp0) . w )
      dpa3_grow( S.tmp, C );
      const double f = dg * g * ( 1.0 - g );
      for(int c = 0; c < C; c++) S.tmp[c] = f * P.gate_w[c];
      dpa4_rms_t( P.gate_norm, C, xp, S.tmp.data(), S.dxp.data() );
    }
    // alpha: dE/dalpha_e = g dout . msg_e, dE/dlogit_e = alpha_e ( A_e - sum alpha A )
    dpa3_grow( S.dalpha, n );
    double mean = 0.0;
    for(int e = 0; e < n; e++)
    {
      const double* m = S.msg.data() + size_t(e) * DC;
      const double* dd = S.dout.data();
      double a = 0.0;
#     pragma omp simd reduction(+:a)
      for(size_t t = 0; t < DC; t++) a += dd[t] * m[t];
      S.dalpha[e] = a * g;
      mean += S.dalpha[e] * S.alpha[e];
    }
    dpa3_grow( S.dmsg, size_t(n) * DC );
    const double ish = 1.0 / std::sqrt( double(C) );
    dpa3_grow( S.tmp, C );
    std::fill( S.tmp.begin(), S.tmp.begin() + C, 0.0 );   // dE/dq_i
    for(int e = 0; e < n; e++)
    {
      const double dl = S.alpha[e] * ( S.dalpha[e] - mean );
      const int j = nbh[e0 + e];
      const double* kj = BB.Kn.data() + size_t(j) * C;
      for(int c = 0; c < C; c++) S.tmp[c] += dl * ish * kj[c];
      {
        DPA3RowLock lk( &B.locks[j] );
        double* gk = BB.GKn.data() + size_t(j) * C;
        for(int c = 0; c < C; c++) gk[c] += dl * ish * S.q0[c];
      }
      double* gr = BB.Grf.data() + size_t(e0 + e) * L * C;
      for(int c = 0; c < C; c++) gr[c] += dl * P.logit_w[c];
      BB.Genv[e0 + e] += dl * 2.0 / BB.envv[e0 + e];
      const double f = S.alpha[e] * g;
      const double* dd = S.dout.data();
      double* dm = S.dmsg.data() + size_t(e) * DC;
      for(size_t t = 0; t < DC; t++) dm[t] = f * dd[t];
    }
    // dE/dq_i -> dE/dxp_i[0]
    dpa3_grow( S.ddso, C );
    dpa3_gemm( false, true, 1, C, C, 1.0, S.tmp.data(), C, P.q.data(), C, 0.0, S.ddso.data(), C );
    dpa4_rms_t( P.qk_norm, C, xp, S.ddso.data(), S.dxp.data() );
    // messages -> local frame, packed D gradients
    dpa3_grow( S.dX, size_t(n) * DmC );
    std::fill( S.dX.begin(), S.dX.begin() + size_t(n) * DmC, 0.0 );
    const int nl = P.lin0.size(), nrot = M.rot_r.size();
    {
      const double* xf = S.XL[nl].data();
      for(int e = 0; e < n; e++)
      {
        const double* Dq = BB.Dq.data() + size_t(e0 + e) * 35;
        const double* dm = S.dmsg.data() + size_t(e) * DC;
        const double* xl = xf + size_t(e) * DmC;
        double* dl = S.dX.data() + size_t(e) * DmC;
        double* gd = BB.GD.data() + size_t(e0 + e) * 35;
        for(int t = 0; t < nrot; t++)
        {
          const int d = M.rot_d[t];
          const int r = M.rot_r[t];
          const double f = Dq[ M.rot_o[t] ] * M.rescale[d];
          const double* dmc = dm + size_t(d) * C;
          const double* xr = xl + size_t(r) * C;
          double* dr = dl + size_t(r) * C;
          double s = 0.0;
#         pragma omp simd reduction(+:s)
          for(int c = 0; c < C; c++) { dr[c] += f * dmc[c]; s += dmc[c] * xr[c]; }
          gd[ M.rot_o[t] ] += s * M.rescale[d];
        }
      }
    }
    // SO(2) layers backward
    const int nm0 = L * C, nmm = M.lmax * C;
    dpa3_grow( S.dY, size_t(n) * DmC ); dpa3_grow( S.dz, size_t(n) * DmC ); dpa3_grow( S.dG, size_t(n) * ( nmm > 0 ? nmm : 1 ) );
    for(int k = nl - 1; k >= 0; k--)
    {
      const double* y = S.Y[k].data();
      double* dy = S.dY.data();
      if( k < nl - 1 )
      {
        const double* gg = S.GG[k].data();
        std::fill( S.dG.begin(), S.dG.begin() + size_t(n) * nmm, 0.0 );
        for(int e = 0; e < n; e++)
        {
          const double* dxo = S.dX.data() + size_t(e) * DmC;
          const double* ye = y + size_t(e) * DmC;
          const double* ge = gg + size_t(e) * nmm;
          double* de = dy + size_t(e) * DmC;
          double* dg0 = S.dG.data() + size_t(e) * nmm;
          const double* sd = S.Sd[k].data() + size_t(e) * C;
          for(int c = 0; c < C; c++) de[c] = dxo[c] * sd[c];
          for(int r = 1; r < Dm; r++)
          {
            const int gi = M.sign_g[r];
            const double* gr = ge + size_t(gi) * C;
            const double* dr = dxo + size_t(r) * C;
            const double* yr = ye + size_t(r) * C;
            double* er = de + size_t(r) * C;
            double* dgg = dg0 + size_t(gi) * C;
#           pragma omp simd
            for(int c = 0; c < C; c++) { er[c] = dr[c] * gr[c]; dgg[c] += dr[c] * yr[c]; }
          }
        }
        // dE/dz = dg sig'(z), to the scalar row of y
        for(size_t t = 0; t < size_t(n) * nmm; t++) S.dG[t] *= gg[t] * ( 1.0 - gg[t] );
        if( nmm > 0 ) dpa3_gemm( false, true, n, C, nmm, 1.0, S.dG.data(), nmm, P.gate[k].data(), nmm, 1.0, dy, int(DmC) );
      }
      else for(size_t t = 0; t < size_t(n) * DmC; t++) dy[t] = S.dX[t];
      // dx_k = dx_{k+1} + dy W^T
      double* dx = S.dz.data();
      for(size_t t = 0; t < size_t(n) * DmC; t++) dx[t] = S.dX[t];
      dpa3_gemm( false, true, n, nm0, nm0, 1.0, dy, int(DmC), P.lin0[k].data(), nm0, 1.0, dx, int(DmC) );
      if( nmm > 0 )
      {
        const double* wm = P.linm[k].data();
        const double* dneg = dy + nm0;
        const double* dpos = dy + nm0 + nmm;
        // neg' = neg u - pos v, pos' = neg v + pos u, weights [nmm][2 nmm]
        dpa3_grow( S.tmp2, size_t(nmm) * 2 * nmm );
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
    // radial mixer backward: dX is dE/dXL[0]
    dpa3_grow( S.dY, size_t(n) * DmC );
    double* dx0 = S.dY.data();
    if( M.nk > 0 )
    {
      dpa3_grow( S.dG, size_t(n) * M.nk );
      std::fill( dx0, dx0 + size_t(n) * DmC, 0.0 );
      std::fill( S.dG.begin(), S.dG.begin() + size_t(n) * M.nk, 0.0 );
      const int nm = M.mix_r.size();
      for(int e = 0; e < n; e++)
      {
        const double* k = S.kern.data() + size_t(e) * M.nk;
        const double* x0 = S.X0.data() + size_t(e) * DmC;
        const double* dxo = S.dX.data() + size_t(e) * DmC;
        double* dxi = dx0 + size_t(e) * DmC;
        double* dk = S.dG.data() + size_t(e) * M.nk;
        for(int t = 0; t < nm; t++)
        {
          const int cc = M.mix_c[t];
          const double f = k[cc];
          const double* dr = dxo + size_t( M.mix_r[t] ) * C;
          const double* xs = x0 + size_t( M.mix_s[t] ) * C;
          double* ds = dxi + size_t( M.mix_s[t] ) * C;
          double s = 0.0;
#         pragma omp simd reduction(+:s)
          for(int c = 0; c < C; c++) { const double db = dr[c] * P.mixb[c]; ds[c] += f * db; s += db * xs[c]; }
          dk[cc] += s;
        }
      }
      dpa3_grow( S.tmp3, size_t(n) * L * C );
      dpa3_gemm( false, true, n, L * C, M.nk, 1.0, S.dG.data(), M.nk, P.mixw.data(), M.nk, 0.0, S.tmp3.data(), L * C );
      for(size_t t = 0; t < size_t(n) * L * C; t++) BB.Grf[ size_t(e0) * L * C + t ] += S.tmp3[t];
    }
    else
    {
      for(int e = 0; e < n; e++)
      {
        const double* x0 = S.X0.data() + size_t(e) * DmC;
        const double* rd = S.rad.data() + size_t(e) * L * C;
        const double* dxo = S.dX.data() + size_t(e) * DmC;
        double* dxi = dx0 + size_t(e) * DmC;
        double* gr = BB.Grf.data() + size_t(e0 + e) * L * C;
        for(int r = 0; r < Dm; r++)
        {
          const int l = M.deg_m[r];
          for(int c = 0; c < C; c++) { dxi[size_t(r)*C + c] = dxo[size_t(r)*C + c] * rd[size_t(l)*C + c]; gr[size_t(l)*C + c] += dxo[size_t(r)*C + c] * x0[size_t(r)*C + c]; }
        }
      }
    }
    // rotation to the edge frame: dE/dxp of the source atoms, packed D gradients
    for(int e = 0; e < n; e++)
    {
      const int j = nbh[e0 + e];
      const double* xs = Xp + size_t(j) * DC;
      const double* Dq = BB.Dq.data() + size_t(e0 + e) * 35;
      const double* dl = dx0 + size_t(e) * DmC;
      double* gd = BB.GD.data() + size_t(e0 + e) * 35;
      dpa3_grow( S.dres, DC );
      std::fill( S.dres.begin(), S.dres.begin() + DC, 0.0 );
      for(int t = 0; t < nrot; t++)
      {
        const int d = M.rot_d[t];
        const double f = Dq[ M.rot_o[t] ];
        const double* dr = dl + size_t( M.rot_r[t] ) * C;
        const double* xc = xs + size_t(d) * C;
        double* gc = S.dres.data() + size_t(d) * C;
        double s = 0.0;
#       pragma omp simd reduction(+:s)
        for(int c = 0; c < C; c++) { gc[c] += f * dr[c]; s += dr[c] * xc[c]; }
        gd[ M.rot_o[t] ] += s;
      }
      DPA3RowLock lk( &B.locks[j] );
      double* gx = BB.GXp.data() + size_t(j) * DC;
      for(size_t t = 0; t < DC; t++) gx[t] += S.dres[t];
    }
    {
      DPA3RowLock lk( &B.locks[i] );
      double* gx = BB.GXp.data() + size_t(i) * DC;
      for(size_t t = 0; t < DC; t++) gx[t] += S.dxp[t];
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
    auto exchange = [&](std::vector<double>& buf, bool to_ghost)
    {
      if( scheme == nullptr ) fatal_error() << "DPA4: no ghost communication scheme (a ghost update must run before deepmd_force)" << std::endl;
      GhostFlatBuffer fb;
      fb.m_data = buf.data(); fb.m_row_size = DC; fb.m_row_stride = DC; fb.m_nrows = n_nodes; fb.m_type = GhostFlatBuffer::ValueType::Double;
      if( to_ghost ) flat_buffer_ghost_exchange<true>( comm, *scheme, poff, fb, 331 );
      else flat_buffer_ghost_exchange<false>( comm, *scheme, poff, fb, 332 );
    };

    // per edge geometry, radial features, initial node features
    BB.envv.assign( nbonds, 0.0 ); BB.rr.assign( nbonds, 0.0 ); BB.rbf.assign( nbonds * M.nrad, 0.0 ); BB.Dq.assign( nbonds * 35, 0.0 );
    BB.rf.assign( nbonds * L * C, 0.0 );
    BB.GD.assign( nbonds * 35, 0.0 ); BB.Grf.assign( nbonds * L * C, 0.0 ); BB.Genv.assign( nbonds, 0.0 ); BB.Grbf.assign( nbonds * M.nrad, 0.0 ); BB.Gv.assign( nbonds * 3, 0.0 );
    BB.isd.assign( n_nodes, 0.0 ); BB.deg.assign( n_nodes, 0.0 );
    BB.X.resize( NB + 1 ); BB.G.resize( NB + 1 );
    for(auto& v : BB.X) v.assign( n_nodes * DC, 0.0 );
    const int Sc = M.s_r + 2 * M.tdim, Em = M.emb, Ax = M.axis;
    BB.seed.assign( n_nodes * size_t(Em) * 4 , 0.0 );      // env seed aggregate (normalized) [4][emb]
    BB.gie.assign( n_nodes * DC, 0.0 );                    // unnormalized initial l > 0 sums
#   pragma omp parallel
    {
      std::vector<double> h( M.rf_norm.scale.size() ), sh( M.rf_norm.scale.size() ), ag( size_t(4) * Em ), Dd( size_t(Em) * Ax ), rt(4), gi( std::max( Em, std::max( M.s_gh, std::max( M.s_r, 1 ) ) ) ), f1( std::max( M.s_r, 1 ) * 4 ), film( 2 * C ), ft( 2 * C );
      std::vector<double> ginp( std::max( Sc, 1 ) ), gh( std::max( M.s_gh, 1 ) ), filmn( 2 * C );
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
          for(int a = 0; a < M.nrad; a++) { const double x = r * M.freq[a]; rb[a] = M.freq[a] * std::sin(x) / x * e7; }
          double qq[4], J[12];
          dpa4_quaternion( v, r, M.eps, M.weps, qq, J );
          dpa4_wigner( M, qq, BB.Dq.data() + size_t(e) * 35 );
          // radial features * envelope
          const int hr = M.hrf;
          dpa3_gemm( false, false, 1, hr, M.nrad, 1.0, rb, M.nrad, M.rf1.data(), hr, 0.0, h.data(), hr );
          dpa4_rms( M.rf_norm, hr, h.data(), sh.data() );
          for(int a = 0; a < hr; a++) sh[a] = sh[a] * dpa4_sig( sh[a] );
          dpa3_gemm( false, false, 1, L * C, hr, en, sh.data(), hr, M.rf2.data(), L * C, 0.0, BB.rf.data() + size_t(e) * L * C, L * C );
          // geometric initial embedding
          const int j = nbh[e];
          const double* rfe = BB.rf.data() + size_t(e) * L * C;
          const double* Dqe = BB.Dq.data() + size_t(e) * 35;
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
        for(int k = 0; k < n; k++)
        {
          const int64_t e = b0 + k;
          const double* v = rij + 3 * e;
          const double r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
          const double inv = 1.0 / std::sqrt( r2 + M.seed_eps * M.seed_eps );
          const double s = BB.envv[e] * inv;
          rt[0] = s; rt[1] = s * v[0] * inv; rt[2] = s * v[1] * inv; rt[3] = s * v[2] * inv;
          const double* rb = BB.rbf.data() + size_t(e) * M.nrad;
          // rbf proj (activations in M.s_a), type embeddings
          dpa3_gemm( false, false, 1, M.s_r, M.nrad, 1.0, rb, M.nrad, M.s_r1.data(), M.s_r, 0.0, f1.data(), M.s_r );
          if( M.s_a[0] ) for(int a = 0; a < M.s_r; a++) f1[a] = f1[a] * dpa4_sig( f1[a] );
          dpa3_gemm( false, false, 1, M.s_r, M.s_r, 1.0, f1.data(), M.s_r, M.s_r2.data(), M.s_r, 0.0, ginp.data(), M.s_r );
          if( M.s_a[1] ) for(int a = 0; a < M.s_r; a++) ginp[a] = ginp[a] * dpa4_sig( ginp[a] );
          const double* ts = M.tt.data() + size_t( B.type[ nbh[e] ] ) * M.tdim;
          const double* td = M.tt.data() + size_t( B.type[i] ) * M.tdim;
          for(int a = 0; a < M.tdim; a++) { ginp[M.s_r + a] = ts[a]; ginp[M.s_r + M.tdim + a] = td[a]; }
          dpa3_gemm( false, false, 1, M.s_gh, Sc, 1.0, ginp.data(), Sc, M.s_g1.data(), M.s_gh, 0.0, gh.data(), M.s_gh );
          if( M.s_a[2] ) for(int a = 0; a < M.s_gh; a++) gh[a] = gh[a] * dpa4_sig( gh[a] );
          dpa3_gemm( false, false, 1, Em, M.s_gh, 1.0, gh.data(), M.s_gh, M.s_g2.data(), Em, 0.0, gi.data(), Em );
          if( M.s_a[3] ) for(int a = 0; a < Em; a++) gi[a] = gi[a] * dpa4_sig( gi[a] );
          for(int x = 0; x < 4; x++) for(int a = 0; a < Em; a++) ga[x*Em + a] += rt[x] * gi[a] * isd;
        }
        // Dd[Em][Ax] = agg^T agg[:, :Ax]
        dpa3_gemm( true, false, Em, Ax, 4, 1.0, ga, Em, ga, Em, 0.0, Dd.data(), Ax );
        dpa3_gemm( false, false, 1, 2 * C, Em * Ax, 1.0, Dd.data(), Em * Ax, M.s_out.data(), 2 * C, 0.0, film.data(), 2 * C );
        if( M.s_a[4] ) for(int a = 0; a < 2 * C; a++) film[a] = film[a] * dpa4_sig( film[a] );
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
      BB.Xp.resize( n_nodes * DC ); BB.Kn.resize( n_nodes * C );
#     pragma omp parallel
      {
        std::vector<double> l0(C);
#       pragma omp for schedule(static)
        for(size_t n = 0; n < n_nodes; n++)
        {
          const double* x = BB.X[b].data() + n * DC;
          double* xp = BB.Xp.data() + n * DC;
          dpa4_so3( M.lmax, P.pre_mix.data(), C, C, x, xp );
          dpa4_rms( P.qk_norm, C, xp, l0.data() );
          dpa3_gemm( false, false, 1, C, C, 1.0, l0.data(), C, P.k.data(), C, 0.0, BB.Kn.data() + n * C, C );
        }
      }
    };

    // forward blocks
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
          dpa4_block_forward( M, M.blocks[b], BB, i, offset[i], BB.ecount[i], nbh, BB.X[b].data() + size_t(i) * DC, S, BB.X[b+1].data() + size_t(i) * DC );
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
      BB.GXp.assign( n_nodes * DC, 0.0 );
      BB.GKn.assign( n_nodes * C, 0.0 );
#     pragma omp parallel
      {
        static thread_local DPA4Scratch S;
        std::vector<double> x2( DC );
#       pragma omp for schedule(dynamic, 1)
        for(int q = 0; q < inum; q++)
        {
          const int i = B.ilist[q];
          dpa4_block_forward( M, M.blocks[b], BB, i, offset[i], BB.ecount[i], nbh, BB.X[b].data() + size_t(i) * DC, S, x2.data() );
          dpa4_block_vjp( M, M.blocks[b], BB, i, offset[i], BB.ecount[i], nbh, BB.G[b+1].data() + size_t(i) * DC, S, BB.dX1.data() + size_t(i) * DC );
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
          dpa3_gemm( false, true, 1, C, C, 1.0, BB.GKn.data() + n * C, C, P.k.data(), C, 0.0, d.data(), C );
          dpa4_rms_t( P.qk_norm, C, BB.Xp.data() + n * DC, d.data(), BB.GXp.data() + n * DC );
        }
      }
      exchange( BB.GXp, false );
      BB.G[b].assign( n_nodes * DC, 0.0 );
#     pragma omp parallel for schedule(static)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q];
        double* g = BB.G[b].data() + size_t(i) * DC;
        dpa4_so3_t( M.lmax, P.pre_mix.data(), C, C, BB.GXp.data() + size_t(i) * DC, g );
        const double* d = BB.dX1.data() + size_t(i) * DC;
        for(size_t t = 0; t < DC; t++) g[t] += d[t];
      }
    }

    mark(3);
    // initial features backward, geometric chain
#   pragma omp parallel
    {
      std::vector<double> ft( 2 * C ), dfl( 2 * C ), Dg( size_t(Em) * Ax ), Ag( size_t(4) * Em ), dgg( Em ), dgh( std::max( M.s_gh, 1 ) ), dgi( std::max( Sc, 1 ) ), dr1( std::max( M.s_r, 1 ) ), dh( M.rf_norm.scale.size() ), h( M.rf_norm.scale.size() ), sh( M.rf_norm.scale.size() ), rtg(4);
      std::vector<double> Fi( 2 * C ), gi( Em ), gh( std::max( M.s_gh, 1 ) ), ginp( std::max( Sc, 1 ) ), f1( std::max( M.s_r, 1 ) ), film( 2 * C ), filmn( 2 * C ), Dd( size_t(Em) * Ax );
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
        if( M.s_a[4] ) for(int a = 0; a < 2 * C; a++) film[a] = film[a] * dpa4_sig( film[a] );
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
          const double* rb = BB.rbf.data() + size_t(e) * M.nrad;
          // forward of g
          dpa3_gemm( false, false, 1, M.s_r, M.nrad, 1.0, rb, M.nrad, M.s_r1.data(), M.s_r, 0.0, f1.data(), M.s_r );
          std::vector<double> f1z( f1 );
          if( M.s_a[0] ) for(int a = 0; a < M.s_r; a++) f1[a] = f1[a] * dpa4_sig( f1[a] );
          dpa3_gemm( false, false, 1, M.s_r, M.s_r, 1.0, f1.data(), M.s_r, M.s_r2.data(), M.s_r, 0.0, ginp.data(), M.s_r );
          std::vector<double> f2z( ginp.begin(), ginp.begin() + M.s_r );
          if( M.s_a[1] ) for(int a = 0; a < M.s_r; a++) ginp[a] = ginp[a] * dpa4_sig( ginp[a] );
          const double* ts = M.tt.data() + size_t( B.type[ nbh[e] ] ) * M.tdim;
          const double* td = M.tt.data() + size_t( B.type[i] ) * M.tdim;
          for(int a = 0; a < M.tdim; a++) { ginp[M.s_r + a] = ts[a]; ginp[M.s_r + M.tdim + a] = td[a]; }
          dpa3_gemm( false, false, 1, M.s_gh, Sc, 1.0, ginp.data(), Sc, M.s_g1.data(), M.s_gh, 0.0, gh.data(), M.s_gh );
          std::vector<double> ghz( gh );
          if( M.s_a[2] ) for(int a = 0; a < M.s_gh; a++) gh[a] = gh[a] * dpa4_sig( gh[a] );
          dpa3_gemm( false, false, 1, Em, M.s_gh, 1.0, gh.data(), M.s_gh, M.s_g2.data(), Em, 0.0, gi.data(), Em );
          std::vector<double> giz( gi );
          if( M.s_a[3] ) for(int a = 0; a < Em; a++) gi[a] = gi[a] * dpa4_sig( gi[a] );
          // dE/drt, dE/dg (outer product rt g isd)
          std::fill( rtg.begin(), rtg.end(), 0.0 );
          std::fill( dgg.begin(), dgg.end(), 0.0 );
          for(int x = 0; x < 4; x++)
          {
            const double* ax = Ag.data() + size_t(x) * Em;
            for(int a = 0; a < Em; a++) { rtg[x] += ax[a] * gi[a] * isd; dgg[a] += ax[a] * rt[x] * isd; }
          }
          // g layers backward
          auto dact = [](double z) { const double s0 = dpa4_sig( z ); return s0 + z * s0 * ( 1.0 - s0 ); };
          if( M.s_a[3] ) for(int a = 0; a < Em; a++) dgg[a] *= dact( giz[a] );
          dpa3_gemm( false, true, 1, M.s_gh, Em, 1.0, dgg.data(), Em, M.s_g2.data(), Em, 0.0, dgh.data(), M.s_gh );
          if( M.s_a[2] ) for(int a = 0; a < M.s_gh; a++) dgh[a] *= dact( ghz[a] );
          dpa3_gemm( false, true, 1, Sc, M.s_gh, 1.0, dgh.data(), M.s_gh, M.s_g1.data(), M.s_gh, 0.0, dgi.data(), Sc );
          if( M.s_a[1] ) for(int a = 0; a < M.s_r; a++) dgi[a] *= dact( f2z[a] );
          dpa3_gemm( false, true, 1, M.s_r, M.s_r, 1.0, dgi.data(), M.s_r, M.s_r2.data(), M.s_r, 0.0, dr1.data(), M.s_r );
          if( M.s_a[0] ) for(int a = 0; a < M.s_r; a++) dr1[a] *= dact( f1z[a] );
          dpa3_gemm( false, true, 1, M.nrad, M.s_r, 1.0, dr1.data(), M.s_r, M.s_r1.data(), M.s_r, 1.0, BB.Grbf.data() + size_t(e) * M.nrad, M.nrad );
          // rt: env part into Genv, vector part into Gv
          BB.Genv[e] += rtg[0] * inv + inv * inv * ( rtg[1]*v[0] + rtg[2]*v[1] + rtg[3]*v[2] );
          const double i3 = en * inv * inv * inv;
          const double i4 = en * inv * inv;
          double* gv = BB.Gv.data() + 3 * e;
          const double rv = rtg[1]*v[0] + rtg[2]*v[1] + rtg[3]*v[2];
          for(int x = 0; x < 3; x++) gv[x] += -rtg[0] * i3 * v[x] + rtg[1 + x] * i4 - 2.0 * rv * v[x] * i4 * inv * inv;
          // geometric initial embedding
          double* gd = BB.GD.data() + size_t(e) * 35;
          double* gr = BB.Grf.data() + size_t(e) * L * C;
          const double* rfe = BB.rf.data() + size_t(e) * L * C;
          const double* Dqe = BB.Dq.data() + size_t(e) * 35;
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
          // radial features = env act( rms( rbf W1 ) ) W2: dE/dfeatures
          const int hr = M.hrf;
          dpa3_gemm( false, false, 1, hr, M.nrad, 1.0, rb, M.nrad, M.rf1.data(), hr, 0.0, h.data(), hr );
          dpa4_rms( M.rf_norm, hr, h.data(), sh.data() );
          // env derivative: rf / env
          double gen = 0.0;
          for(int t = 0; t < L * C; t++) gen += gr[t] * ( en != 0.0 ? rfe[t] / en : 0.0 );
          BB.Genv[e] += gen;
          dpa3_gemm( false, true, 1, hr, L * C, en, gr, L * C, M.rf2.data(), L * C, 0.0, dh.data(), hr );
          for(int a = 0; a < hr; a++) dh[a] *= dact( sh[a] );
          std::fill( f1.begin(), f1.end(), 0.0 );
          std::vector<double> dhh( hr, 0.0 );
          dpa4_rms_t( M.rf_norm, hr, h.data(), dh.data(), dhh.data() );
          dpa3_gemm( false, true, 1, M.nrad, hr, 1.0, dhh.data(), hr, M.rf1.data(), hr, 1.0, BB.Grbf.data() + size_t(e) * M.nrad, M.nrad );
          // radial: dE/dr from env and rbf, Wigner chain
          const double r = BB.rr[e];
          double e5, d5, e7, d7;
          dpa4_envelope( M.env_c, r, M.rc, e5, d5 );
          dpa4_envelope( M.rbf_c, r, M.rc, e7, d7 );
          double gr0 = BB.Genv[e] * d5;
          const double* grb = BB.Grbf.data() + size_t(e) * M.nrad;
          for(int a = 0; a < M.nrad; a++)
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
