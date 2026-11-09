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

// DeePMD force functor (descriptors se_e2_a, se_e2_r, se_e3), one central atom i per call, full
// neighbor list within rcut. se_e2_r and se_e3 are described where they differ, below.
//
// Forward, se_e2_a (deepmd/dpmodel/descriptor/se_e2_a.py):
//   R_j    = sw(r) ( 1/r, d/r^2 ),  d = r_j - r_i            environment matrix row
//   Rh_j   = ( R_j - davg[ti][tj] ) / dstd[ti][tj]
//   G_j    = embedding[tj or ti,tj]( Rh_j[0] )               ng values
//   gr     = ( sum_j G_j x Rh_j + empty slots ) / nnei      ng x 4
//   D      = gr gr[:axis_neuron]^T                            ng x axis_neuron
//   E_i    = fitting[ti](D) + bias[ti]
// The neighbor list of type t holds sel[t] slots: the nearest sel[t] neighbors are kept, empty
// slots contribute G(pad) x Rh(pad) with Rh(pad) = -davg/dstd (a constant, no force).
//
// Backward: hand-written vector-Jacobian products, dE/dD -> dE/dgr -> per neighbor dE/dRh, dE/dG
// -> dE/ds (embedding backward, recomputed) -> dE/dd. The force dE/dd_j acts on j (-) and i (+);
// forces on neighbors (ghosts included) are atomically scattered and folded back by
// update_force_energy_from_ghost. Per-atom virial: sum_j tensor(f_j, d_j), on the central atom.
//
// Cooperative execution (one GPU block per central atom, BlockSharedComputeBuffer): the per-atom
// scratch is in the compute buffer extension, i.e. block shared memory; neighbors are distributed
// over the threads of the block for the embedding networks (thread-local work, block atomics into
// gr), output neurons for the fitting network. On the host the same code runs with one thread
// (ONIKA_CU_BLOCK_SIZE = 1, block synchronizations and block atomics are no-ops).

#include <cmath>
#include <type_traits>
#include <exanb/compute/compute_pair_buffer.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <exanb/compute/compute_pair_traits.h>
#include <exanb/core/concurent_add_contributions.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/cuda/cuda.h>

#include "deepmd_model.h"

#ifndef XSTAMP_DEEPMD_MAX_NEIGHBORS
#define XSTAMP_DEEPMD_MAX_NEIGHBORS 256
#endif

namespace exaStamp
{
  using namespace exanb;

  static constexpr size_t DP_MAX_NEIGHBORS = XSTAMP_DEEPMD_MAX_NEIGHBORS;

  // fitting network activations of the forward pass (z per layer, ping-pong outputs x0 / x1)
  struct DeepmdFitScratch
  {
    double z[DP_MAX_FIT_LAYERS][DP_MAX_FIT_WIDTH];
    double x0[DP_MAX_FIT_WIDTH];
    double x1[DP_MAX_FIT_WIDTH];
  };

  // dE/dgr is needed only once the fitting backward pass is done: it reuses the fitting scratch
  union DeepmdFitOrGrad
  {
    DeepmdFitScratch net;
    double dgr[DP_MAX_EMB_WIDTH * 4];   // two-body descriptors
    double dv[DP_MAX_NEIGHBORS * 3];    // se_e3 : dE/dRh[1:4] per neighbor
  };

  // se_e3 / se_t_tebd cache the normalized rows of the kept neighbors at the end of the descriptor array,
  // whose head holds D and dE/dD (ng values, + tebd_dim for se_t_tebd); se_t_tebd also keeps the switch values and
  // their dE/dsw before the rows
  static constexpr int DP_E3_ROWS_OFFSET = DP_MAX_DESCRIPTOR - 4 * int(DP_MAX_NEIGHBORS);
  static constexpr int DP_T_SW_OFFSET = DP_E3_ROWS_OFFSET - 2 * int(DP_MAX_NEIGHBORS);
  static_assert( DP_T_SW_OFFSET >= DP_MAX_EMB_WIDTH + 64 , "XSTAMP_DEEPMD_MAX_DESCRIPTOR must hold 6 x XSTAMP_DEEPMD_MAX_NEIGHBORS + XSTAMP_DEEPMD_MAX_EMB_WIDTH + 64 values" );

  // Neighbor types and per-atom scratch, in block shared memory on GPU
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) DeepmdComputeBuffer
  {
    alignas(onika::memory::DEFAULT_ALIGNMENT) int type[DP_MAX_NEIGHBORS];  // model type of the kept neighbors
    alignas(onika::memory::DEFAULT_ALIGNMENT) double gr[DP_MAX_EMB_WIDTH * 4];
    alignas(onika::memory::DEFAULT_ALIGNMENT) double desc[DP_MAX_DESCRIPTOR]; // D, then dE/dD
    alignas(onika::memory::DEFAULT_ALIGNMENT) DeepmdFitOrGrad fit;
    double fi[3];
    double vir[9];
    int count[DP_MAX_TYPES];
    int ns;
  };

  // Plain pointers to the model, built from a DeepmdContext
  struct DeepmdModelView
  {
    const double * weights = nullptr;
    const DPLayer * layers = nullptr;
    const DPNet * embedding = nullptr;
    const DPNet * fitting = nullptr;
    const double * davg = nullptr;
    const double * dstd = nullptr;
    const double * pad = nullptr;
    const double * energy_bias = nullptr;
    const int * sel = nullptr;
    const int * species_type = nullptr;
    int descriptor = DP_SE_E2_A;
    int ntypes = 0;
    int ng = 0;
    int axis_neuron = 0;
    int nnei = 0;
    bool type_one_side = true;
    bool fitting_mixed_types = false;
    double rcut = 0.0;
    double rcut_smth = 0.0;
    double protection = 0.0;
    double conv_energy = 1.0;  // eV -> internal energy unit
    // mixed neighbor list, type embeddings (se_t_tebd)
    bool mixed = false;
    bool strip = false;
    bool smooth = true;
    bool concat_output_tebd = false;
    int tebd_dim = 0;
    const double * tebd = nullptr;
    const double * strip_table = nullptr;

    inline DeepmdModelView() = default;
    inline DeepmdModelView(const DeepmdContext& c, double conv)
      : weights(c.weights.data()), layers(c.layers.data()), embedding(c.embedding.data()), fitting(c.fitting.data())
      , davg(c.davg.data()), dstd(c.dstd.data()), pad(c.pad.data()), energy_bias(c.energy_bias.data())
      , sel(c.sel.data()), species_type(c.species_type.data())
      , descriptor(c.descriptor), ntypes(c.ntypes), ng(c.ng), axis_neuron(c.axis_neuron), nnei(c.nnei)
      , type_one_side(c.type_one_side), fitting_mixed_types(c.fitting_mixed_types)
      , rcut(c.rcut), rcut_smth(c.rcut_smth), protection(c.protection), conv_energy(conv)
      , mixed(c.mixed), strip(c.strip), smooth(c.attn_smooth), concat_output_tebd(c.concat_output_tebd), tebd_dim(c.tebd_dim)
      , tebd(c.tebd.data()), strip_table(c.strip_table.data()) {}
  };

  // smooth switch w(r) = 1 - 10u^3 + 15u^4 - 6u^5, u = (r - rs) / (rc - rs) clipped to [0,1], and dw/dr
  ONIKA_HOST_DEVICE_FUNC inline void dp_switch(double r, double rs, double rc, double& w, double& dw)
  {
    if( r <= rs ) { w = 1.0; dw = 0.0; return; }
    if( r >= rc ) { w = 0.0; dw = 0.0; return; }
    const double u = ( r - rs ) / ( rc - rs );
    const double u2 = u * u;
    w = u2 * u * ( -6.0 * u2 + 15.0 * u - 10.0 ) + 1.0;
    dw = u2 * ( -30.0 * u2 + 60.0 * u - 30.0 ) / ( rc - rs );
  }

  struct DeepmdForceOp
  {
    DeepmdModelView m;

    // normalized environment matrix row of neighbor d, and the derivative factors of the raw row
    // R0 = w/q, Rk = w d_k / q^2, q = r + protection
    struct EnvRow
    {
      double rh[4];
      double r, w, dw, q;
    };

    ONIKA_HOST_DEVICE_FUNC inline void env_row(const Vec3d& d, double d2, int ti, int tj, EnvRow& e) const
    {
      e.r = std::sqrt(d2);
      dp_switch(e.r, m.rcut_smth, m.rcut, e.w, e.dw);
      e.q = e.r + m.protection;
      const double raw[4] = { e.w / e.q, e.w * d.x / ( e.q * e.q ), e.w * d.y / ( e.q * e.q ), e.w * d.z / ( e.q * e.q ) };
      const double* avg = m.davg + ( ti * m.ntypes + tj ) * 4;
      const double* std = m.dstd + ( ti * m.ntypes + tj ) * 4;
      for(int c = 0; c < 4; c++) e.rh[c] = ( raw[c] - avg[c] ) / std[c];
    }

    // dE/dd from dE/dRh (normalized row) for neighbor d
    ONIKA_HOST_DEVICE_FUNC inline Vec3d env_row_grad(const Vec3d& d, const EnvRow& e, const double* drh, int ti, int tj) const
    {
      const double* std = m.dstd + ( ti * m.ntypes + tj ) * 4;
      const double g0 = drh[0] / std[0];
      const Vec3d gk = { drh[1] / std[1], drh[2] / std[2], drh[3] / std[3] };
      const double q2 = e.q * e.q;
      // dR0/dd = ( w'/q - w/q^2 ) d/r ;  dRk/dd_m = w delta_km / q^2 + d_k ( w'/q^2 - 2w/q^3 ) d_m / r
      const double radial = g0 * ( e.dw / e.q - e.w / q2 ) + dot(gk, d) * ( e.dw / q2 - 2.0 * e.w / ( q2 * e.q ) );
      return gk * ( e.w / q2 ) + d * ( radial / e.r );
    }

    // two-body descriptors : by neighbor type, or by (center, neighbor) type pair
    ONIKA_HOST_DEVICE_FUNC inline const DPNet& embedding_net(int ti, int tj) const
    {
      return m.embedding[ m.type_one_side ? tj : ti * m.ntypes + tj ];
    }

    template <class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3d &virial, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    // fitting network forward pass, threads over output neurons; returns the output (block shared)
    ONIKA_HOST_DEVICE_FUNC inline const double* fit_forward(const DPNet& net, const double* in, DeepmdFitScratch& s, int tid, int bs) const
    {
      const double* xin = in;
      for(int l = 0; l < net.nlayers; l++)
      {
        const DPLayer& L = m.layers[net.first_layer + l];
        const double* w = m.weights + L.w;
        const double* b = m.weights + L.b;
        double* y = ( l % 2 == 0 ) ? s.x0 : s.x1;
        for(int o = tid; o < L.nout; o += bs)
        {
          double z = b[o];
          for(int i = 0; i < L.nin; i++) z += xin[i] * w[ size_t(i) * L.nout + o ];
          s.z[l][o] = z;
          double v = dp_activation(L.activation, z);
          if( L.idt != size_t(-1) ) v *= m.weights[L.idt + o];
          if( L.resnet == 1 ) v += xin[o];
          else if( L.resnet == 2 ) v += xin[ o % L.nin ];
          y[o] = v;
        }
        ONIKA_CU_BLOCK_SYNC();
        xin = y;
      }
      return xin;
    }

    // fitting network backward pass from dE/dout = 1: dE/din in din. dE/dz overwrites z, layer by layer
    ONIKA_HOST_DEVICE_FUNC inline void fit_backward(const DPNet& net, DeepmdFitScratch& s, double* din, int tid, int bs) const
    {
      double* dy = s.x0;
      if( tid == 0 ) dy[0] = 1.0;
      ONIKA_CU_BLOCK_SYNC();
      for(int l = net.nlayers - 1; l >= 0; l--)
      {
        const DPLayer& L = m.layers[net.first_layer + l];
        const double* w = m.weights + L.w;
        double* dz = s.z[l];
        for(int o = tid; o < L.nout; o += bs)
        {
          double g = dy[o] * dp_activation_grad(L.activation, dz[o]);
          if( L.idt != size_t(-1) ) g *= m.weights[L.idt + o];
          dz[o] = g;
        }
        ONIKA_CU_BLOCK_SYNC();
        double* dx = ( l == 0 ) ? din : ( ( dy == s.x0 ) ? s.x1 : s.x0 );
        for(int i = tid; i < L.nin; i += bs)
        {
          const double* wi = w + size_t(i) * L.nout;
          double v = 0.0;
          for(int o = 0; o < L.nout; o++) v += wi[o] * dz[o];
          if( L.resnet == 1 ) v += dy[i];
          else if( L.resnet == 2 ) v += dy[i] + dy[i + L.nin];
          dx[i] = v;
        }
        ONIKA_CU_BLOCK_SYNC();
        dy = dx;
      }
    }

    // forward part: neighbor selection, environment rows, descriptor D in buf.ext.desc (+ type embedding);
    // the state the backward part needs stays in buf.ext
    template <class ComputeBufferT>
    ONIKA_HOST_DEVICE_FUNC inline void forward_desc(size_t n, ComputeBufferT &buf, int type) const
    {
      const int tid = ONIKA_CU_THREAD_IDX;
      const int bs = ONIKA_CU_BLOCK_SIZE;
      auto& x = buf.ext;
      if( n > ComputeBufferT::MaxNeighbors )
      {
        if( tid == 0 ) printf("deepmd_force: %d neighbors within the cutoff exceed buffer capacity (%d), rebuild with larger XSTAMP_DEEPMD_MAX_NEIGHBORS\n", int(n), int(ComputeBufferT::MaxNeighbors));
        ONIKA_CU_ABORT();
      }

      const int ti = m.species_type[type];
      const int nt = m.ntypes;
      const int ng = m.ng;
      const int m2 = m.axis_neuron;
      const double inv_nnei = 1.0 / m.nnei;

      // ---- neighbors inside rcut, mapped to model types; per type, keep the sel[t] nearest (thread 0)
      if( tid == 0 )
      {
        const double rc2 = m.rcut * m.rcut;
        for(int t = 0; t < nt; t++) x.count[t] = 0;
        int ns = 0;
        for(size_t jj = 0; jj < n; ++jj)
        {
          if( buf.d2[jj] < rc2 )
          {
            const int tj = m.species_type[ buf.nbh_pt[jj][field::type] ];
            if( size_t(ns) != jj ) { buf.copy(jj, ns); }
            x.type[ns] = tj;
            ++ x.count[tj];
            ++ ns;
          }
        }
        while( m.mixed && ns > m.nnei ) // mixed list: drop the farthest neighbor, any type
        {
          int far = 0;
          for(int jj = 1; jj < ns; jj++) if( buf.d2[jj] > buf.d2[far] ) far = jj;
          -- x.count[ x.type[far] ];
          -- ns;
          if( far != ns ) { buf.copy(ns, far); x.type[far] = x.type[ns]; }
        }
        for(int t = 0; t < nt && !m.mixed; t++)
        {
          while( x.count[t] > m.sel[t] ) // drop the farthest neighbor of type t
          {
            int far = -1;
            for(int jj = 0; jj < ns; jj++) if( x.type[jj] == t && ( far < 0 || buf.d2[jj] > buf.d2[far] ) ) far = jj;
            -- ns;
            if( far != ns ) { buf.copy(ns, far); x.type[far] = x.type[ns]; }
            -- x.count[t];
          }
        }
        x.ns = ns;
        for(int k = 0; k < 3; k++) x.fi[k] = 0.0;
        for(int k = 0; k < 9; k++) x.vir[k] = 0.0;
      }
      ONIKA_CU_BLOCK_SYNC();
      const int ns = x.ns;

      // ---- forward
      //   se_e2_a : gr = ( sum_j G_j x Rh_j + empty slots ) / nnei,  D = gr gr[:m2]^T
      //   se_e2_r : gr[:,0] = sum_j G_j + empty slots,  D = gr[:,0] / ( 5 nnei )
      //   se_e3   : D = sum_{j,k, tj<=tk} th_jk G_tjtk(th_jk) / ( sel_tj sel_tk ),  th_jk = Rh_j[1:4] . Rh_k[1:4]
      const int desc_kind = m.descriptor;
      double* rows = x.desc + DP_E3_ROWS_OFFSET;
      for(int k = tid; k < ng * 4; k += bs)
      {
        double s = 0.0;
        if( !m.mixed ) for(int t = 0; t < nt; t++) s += ( m.sel[t] - x.count[t] ) * m.pad[ size_t( ti * nt + t ) * ng * 4 + k ];
        x.gr[k] = s;
      }
      ONIKA_CU_BLOCK_SYNC();

      DPEmbeddingWork wk; // thread-local
      if( desc_kind != DP_SE_E3 && desc_kind != DP_SE_T_TEBD )
      {
        for(int jj = tid; jj < ns; jj += bs)
        {
          const int tj = x.type[jj];
          const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] }; // r_j - r_i
          EnvRow e;
          env_row(d, buf.d2[jj], ti, tj, e);
          const double* G = dp_mlp_forward(embedding_net(ti, tj), m.layers, m.weights, e.rh, wk);
          const double one[4] = { 1.0, 0.0, 0.0, 0.0 };
          const double* prod = ( desc_kind == DP_SE_E2_R ) ? one : e.rh;
          for(int a = 0; a < ng; a++) for(int c = 0; c < 4; c++) ONIKA_CU_BLOCK_ATOMIC_ADD( x.gr[a*4 + c] , G[a] * prod[c] );
        }
      }
      else if( desc_kind == DP_SE_T_TEBD )
      {
        // D = sum_{j,k} th_jk G_jk(th_jk) / nnei^2 over all ordered pairs, G_jk = net_{tj,tk}(th) (concat) or
        //     g(th) ( 1 + T_{tj,tk} [sw_j sw_k] ) (strip)
        double* swv = x.desc + DP_T_SW_OFFSET;
        for(int jj = tid; jj < ns; jj += bs)
        {
          const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
          EnvRow e;
          env_row(d, buf.d2[jj], ti, x.type[jj], e);
          for(int c = 0; c < 4; c++) rows[jj*4 + c] = e.rh[c];
          swv[jj] = e.w;
          swv[DP_MAX_NEIGHBORS + jj] = 0.0;   // dE/dsw
        }
        ONIKA_CU_BLOCK_SYNC();
        const double w = 1.0 / ( double(m.nnei) * m.nnei );
        for(int p = tid; p < ns * ns; p += bs)
        {
          const int j = p / ns, k = p % ns, tj = x.type[j], tk = x.type[k];
          const double th = rows[j*4+1] * rows[k*4+1] + rows[j*4+2] * rows[k*4+2] + rows[j*4+3] * rows[k*4+3];
          const double* G = dp_mlp_forward(m.embedding[ m.strip ? 0 : tj * nt + tk ], m.layers, m.weights, &th, wk);
          if( !m.strip ) { for(int a = 0; a < ng; a++) ONIKA_CU_BLOCK_ATOMIC_ADD( x.gr[a*4] , G[a] * th * w ); }
          else
          {
            const double* T = m.strip_table + size_t( tj * nt + tk ) * ng;
            const double f = m.smooth ? swv[j] * swv[k] : 1.0;
            for(int a = 0; a < ng; a++) ONIKA_CU_BLOCK_ATOMIC_ADD( x.gr[a*4] , G[a] * ( 1.0 + T[a] * f ) * th * w );
          }
        }
      }
      else
      {
        for(int jj = tid; jj < ns; jj += bs)
        {
          const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
          EnvRow e;
          env_row(d, buf.d2[jj], ti, x.type[jj], e);
          for(int c = 0; c < 4; c++) rows[jj*4 + c] = e.rh[c];
        }
        ONIKA_CU_BLOCK_SYNC();
        for(int p = tid; p < ns * ns; p += bs)
        {
          const int j = p / ns, k = p % ns, tj = x.type[j], tk = x.type[k];
          if( tj > tk ) continue;
          const double th = rows[j*4+1] * rows[k*4+1] + rows[j*4+2] * rows[k*4+2] + rows[j*4+3] * rows[k*4+3];
          const double wp = th / ( double(m.sel[tj]) * m.sel[tk] );
          const double* G = dp_mlp_forward(m.embedding[ tj * nt + tk ], m.layers, m.weights, &th, wk);
          for(int a = 0; a < ng; a++) ONIKA_CU_BLOCK_ATOMIC_ADD( x.gr[a*4] , G[a] * wp );
        }
      }
      ONIKA_CU_BLOCK_SYNC();

      if( desc_kind == DP_SE_E2_A )
      {
        for(int k = tid; k < ng * 4; k += bs) x.gr[k] *= inv_nnei;
        ONIKA_CU_BLOCK_SYNC();
        for(int k = tid; k < ng * m2; k += bs)
        {
          const double* ga = x.gr + ( k / m2 ) * 4;
          const double* gb = x.gr + ( k % m2 ) * 4;
          x.desc[k] = ga[0] * gb[0] + ga[1] * gb[1] + ga[2] * gb[2] + ga[3] * gb[3];
        }
      }
      else
      {
        const double scale = ( desc_kind == DP_SE_E2_R ) ? inv_nnei / 5.0 : 1.0;
        for(int a = tid; a < ng; a += bs) x.desc[a] = x.gr[a*4] * scale;
        if( desc_kind == DP_SE_T_TEBD && m.concat_output_tebd ) for(int c = tid; c < m.tebd_dim; c += bs) x.desc[ng + c] = m.tebd[ ti * m.tebd_dim + c ];
      }
      ONIKA_CU_BLOCK_SYNC();
    }

    // backward part, from dE/dD in buf.ext.desc: forces on the neighbors (atomic), central atom force and virial
    // block-reduced into buf.ext.fi / buf.ext.vir
    template <bool compute_virial, class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void backward_desc(ComputeBufferT &buf, int type, CellParticlesT cells) const
    {
      const int tid = ONIKA_CU_THREAD_IDX;
      const int bs = ONIKA_CU_BLOCK_SIZE;
      auto& x = buf.ext;
      const int ti = m.species_type[type];
      const int nt = m.ntypes;
      const int ng = m.ng;
      const int m2 = m.axis_neuron;
      const double inv_nnei = 1.0 / m.nnei;
      const int ns = x.ns;
      const int desc_kind = m.descriptor;
      double* rows = x.desc + DP_E3_ROWS_OFFSET;
      DPEmbeddingWork wk; // thread-local
      Vec3d fi = { 0.0, 0.0, 0.0 };
      Mat3d vir = { 0.,0.,0., 0.,0.,0., 0.,0.,0. };
      const double conv = m.conv_energy;
      double g0[DP_MAX_EMB_WIDTH], g1[DP_MAX_EMB_WIDTH];

      // force on neighbor jj from dE/dRh, its scatter, and the central atom / virial accumulation
      auto neighbor_force = [&](int jj, const Vec3d& d, const EnvRow& e, const double* drh, double dsw_j = 0.0)
      {
        const Vec3d g = ( env_row_grad(d, e, drh, ti, x.type[jj]) + d * ( dsw_j * e.dw / e.r ) ) * conv; // dE/dd
        const Vec3d fj = g * -1.0;
        fi = fi + g;
        if constexpr ( compute_virial ) { vir = vir + tensor(fj, d); }
        size_t cell_j = 0, p_j = 0;
        buf.nbh.get(jj, cell_j, p_j);
        atomic_add_contribution(cells[cell_j][field::fx][p_j], fj.x);
        atomic_add_contribution(cells[cell_j][field::fy][p_j], fj.y);
        atomic_add_contribution(cells[cell_j][field::fz][p_j], fj.z);
      };

      if( desc_kind != DP_SE_E3 && desc_kind != DP_SE_T_TEBD )
      {
        // dE/dgr, scaled to dE/d of the unnormalized sum
        for(int k = tid; k < ng * 4; k += bs)
        {
          const int a = k / 4, c = k % 4;
          double s = 0.0;
          if( desc_kind == DP_SE_E2_A )
          {
            for(int b = 0; b < m2; b++) s += x.desc[a*m2 + b] * x.gr[b*4 + c];
            if( a < m2 ) for(int a2 = 0; a2 < ng; a2++) s += x.desc[a2*m2 + a] * x.gr[a2*4 + c];
            s *= inv_nnei;
          }
          else if( c == 0 ) s = x.desc[a] * inv_nnei / 5.0;
          x.fit.dgr[k] = s;
        }
        ONIKA_CU_BLOCK_SYNC();
        const double* dgr = x.fit.dgr;

        // per neighbor : recompute the embedding, dE/dRh and dE/dG, embedding backward, dE/dd
        for(int jj = tid; jj < ns; jj += bs)
        {
          const int tj = x.type[jj];
          const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
          EnvRow e;
          env_row(d, buf.d2[jj], ti, tj, e);
          const DPNet& emb = embedding_net(ti, tj);
          const double* G = dp_mlp_forward(emb, m.layers, m.weights, e.rh, wk);
          const double one[4] = { 1.0, 0.0, 0.0, 0.0 };
          const double* prod = ( desc_kind == DP_SE_E2_R ) ? one : e.rh;

          double drh[4] = { 0.0, 0.0, 0.0, 0.0 };
          for(int a = 0; a < ng; a++)
          {
            double s = 0.0;
            for(int c = 0; c < 4; c++) { s += dgr[a*4 + c] * prod[c]; if( desc_kind == DP_SE_E2_A ) drh[c] += dgr[a*4 + c] * G[a]; }
            g0[a] = s; // dE/dG
          }
          double ds = 0.0;
          dp_mlp_backward(emb, m.layers, m.weights, e.rh, wk, g0, &ds, g0, g1);
          drh[0] += ds;
          neighbor_force(jj, d, e, drh);
        }
      }
      else if( desc_kind == DP_SE_T_TEBD )
      {
        // pairs : G = g f, f = 1 + T [sw_j sw_k] ; dE/dth = w ( dD . G + th dD . f dg/dth ),
        //         dE/dsw_j += w th dD . g T sw_k (strip, smooth), symmetrically for k
        double* swv = x.desc + DP_T_SW_OFFSET;
        double* dsw = swv + DP_MAX_NEIGHBORS;
        for(int k = tid; k < ns * 3; k += bs) x.fit.dv[k] = 0.0;
        ONIKA_CU_BLOCK_SYNC();
        const double w = 1.0 / ( double(m.nnei) * m.nnei );
        for(int p = tid; p < ns * ns; p += bs)
        {
          const int j = p / ns, k = p % ns, tj = x.type[j], tk = x.type[k];
          const double th = rows[j*4+1] * rows[k*4+1] + rows[j*4+2] * rows[k*4+2] + rows[j*4+3] * rows[k*4+3];
          const DPNet& emb = m.embedding[ m.strip ? 0 : tj * nt + tk ];
          const double* G = dp_mlp_forward(emb, m.layers, m.weights, &th, wk);
          const double* T = m.strip ? m.strip_table + size_t( tj * nt + tk ) * ng : nullptr;
          const double f2 = ( m.strip && m.smooth ) ? swv[j] * swv[k] : 1.0;
          double dth = 0.0, dsf = 0.0;
          for(int a = 0; a < ng; a++)
          {
            const double f = m.strip ? 1.0 + T[a] * f2 : 1.0;
            dth += x.desc[a] * G[a] * f * w;
            g0[a] = x.desc[a] * th * w * f;
            if( m.strip && m.smooth ) dsf += x.desc[a] * G[a] * T[a] * th * w;
          }
          double dg = 0.0;
          dp_mlp_backward(emb, m.layers, m.weights, &th, wk, g0, &dg, g0, g1);
          dth += dg;
          for(int c = 0; c < 3; c++)
          {
            ONIKA_CU_BLOCK_ATOMIC_ADD( x.fit.dv[j*3 + c] , dth * rows[k*4 + 1 + c] );
            ONIKA_CU_BLOCK_ATOMIC_ADD( x.fit.dv[k*3 + c] , dth * rows[j*4 + 1 + c] );
          }
          if( dsf != 0.0 )
          {
            ONIKA_CU_BLOCK_ATOMIC_ADD( dsw[j] , dsf * swv[k] );
            ONIKA_CU_BLOCK_ATOMIC_ADD( dsw[k] , dsf * swv[j] );
          }
        }
        ONIKA_CU_BLOCK_SYNC();
        for(int jj = tid; jj < ns; jj += bs)
        {
          const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
          EnvRow e;
          env_row(d, buf.d2[jj], ti, x.type[jj], e);
          const double drh[4] = { 0.0, x.fit.dv[jj*3], x.fit.dv[jj*3 + 1], x.fit.dv[jj*3 + 2] };
          neighbor_force(jj, d, e, drh, dsw[jj]);
        }
      }
      else
      {
        // pairs : dE/dth = w ( dD . G + dD . dG/dth th ), accumulated into dE/dRh[1:4] of both neighbors
        for(int k = tid; k < ns * 3; k += bs) x.fit.dv[k] = 0.0;
        ONIKA_CU_BLOCK_SYNC();
        for(int p = tid; p < ns * ns; p += bs)
        {
          const int j = p / ns, k = p % ns, tj = x.type[j], tk = x.type[k];
          if( tj > tk ) continue;
          const double th = rows[j*4+1] * rows[k*4+1] + rows[j*4+2] * rows[k*4+2] + rows[j*4+3] * rows[k*4+3];
          const double w = 1.0 / ( double(m.sel[tj]) * m.sel[tk] );
          const DPNet& emb = m.embedding[ tj * nt + tk ];
          const double* G = dp_mlp_forward(emb, m.layers, m.weights, &th, wk);
          double dth = 0.0;
          for(int a = 0; a < ng; a++) { dth += x.desc[a] * G[a] * w; g0[a] = x.desc[a] * th * w; }
          double dg = 0.0;
          dp_mlp_backward(emb, m.layers, m.weights, &th, wk, g0, &dg, g0, g1);
          dth += dg;
          for(int c = 0; c < 3; c++)
          {
            ONIKA_CU_BLOCK_ATOMIC_ADD( x.fit.dv[j*3 + c] , dth * rows[k*4 + 1 + c] );
            ONIKA_CU_BLOCK_ATOMIC_ADD( x.fit.dv[k*3 + c] , dth * rows[j*4 + 1 + c] );
          }
        }
        ONIKA_CU_BLOCK_SYNC();
        for(int jj = tid; jj < ns; jj += bs)
        {
          const Vec3d d = { buf.drx[jj], buf.dry[jj], buf.drz[jj] };
          EnvRow e;
          env_row(d, buf.d2[jj], ti, x.type[jj], e);
          const double drh[4] = { 0.0, x.fit.dv[jj*3], x.fit.dv[jj*3 + 1], x.fit.dv[jj*3 + 2] };
          neighbor_force(jj, d, e, drh);
        }
      }

      // ---- block reduction of the central atom force and virial, written by thread 0
      ONIKA_CU_BLOCK_ATOMIC_ADD( x.fi[0] , fi.x );
      ONIKA_CU_BLOCK_ATOMIC_ADD( x.fi[1] , fi.y );
      ONIKA_CU_BLOCK_ATOMIC_ADD( x.fi[2] , fi.z );
      if constexpr ( compute_virial )
      {
        const double v[9] = { vir.m11, vir.m12, vir.m13, vir.m21, vir.m22, vir.m23, vir.m31, vir.m32, vir.m33 };
        for(int k = 0; k < 9; k++) ONIKA_CU_BLOCK_ATOMIC_ADD( x.vir[k] , v[k] );
      }
      ONIKA_CU_BLOCK_SYNC();
    }

    // All force updates are atomics, locks are not needed.
    template <class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT>
    ONIKA_HOST_DEVICE_FUNC inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3dT &virial, CellParticlesT cells, GridCellLocksT, ParticleLockT &) const
    {
      static constexpr bool compute_virial = std::is_same_v<Mat3dT, Mat3d>;
      const int tid = ONIKA_CU_THREAD_IDX;
      const int bs = ONIKA_CU_BLOCK_SIZE;
      auto& x = buf.ext;
      const int ti = m.species_type[type];
      forward_desc(n, buf, type);
      const DPNet& fit = m.fitting[ m.fitting_mixed_types ? 0 : ti ];
      const double e_atom = fit_forward(fit, x.desc, x.fit.net, tid, bs)[0] + m.energy_bias[ti];
      ONIKA_CU_BLOCK_SYNC();
      // ---- backward : dE/dD (in desc), then dE/d(sum over neighbors) per descriptor
      fit_backward(fit, x.fit.net, x.desc, tid, bs);
      backward_desc<compute_virial>(buf, type, cells);
      if( tid == 0 )
      {
        // central atom force may be concurrently updated by other atoms
        atomic_add_contribution(fx, x.fi[0]);
        atomic_add_contribution(fy, x.fi[1]);
        atomic_add_contribution(fz, x.fi[2]);
        en += e_atom * m.conv_energy;
        if constexpr ( compute_virial )
        {
          virial += Mat3d{ x.vir[0], x.vir[1], x.vir[2], x.vir[3], x.vir[4], x.vir[5], x.vir[6], x.vir[7], x.vir[8] };
        }
      }
    }
  };

}

namespace exanb
{
  template <> struct ComputePairTraits<exaStamp::DeepmdForceOp>
  {
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool BlockSharedComputeBuffer = true; // one block per central atom
    static inline constexpr bool CudaCompatible = true;
    static inline constexpr bool HasParticleContextStart = false;
    static inline constexpr bool HasParticleContextStop = false;
    static inline constexpr bool HasParticleContext = false;
    static inline constexpr bool RequiresNbhOptionalData = true;
    static inline constexpr bool ComputeEmptyNeighborhood = true; // an isolated atom has an energy (bias and empty slots)
  };
}
