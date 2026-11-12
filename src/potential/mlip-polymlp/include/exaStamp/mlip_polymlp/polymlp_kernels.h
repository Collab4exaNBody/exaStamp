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

// Polynomial MLP onika kernels, CPU and GPU, fp64. Per chunk of central atoms ordered by type:
//   1. basis    (a block per atom): its threads over a batch of neighbours (gaussians, Y_lm into the block scratch),
//                then over the sources (sums over the batch, in a fixed order) -> A [source][atom];
//   2. readout  (a thread per atom): features D, E_i from the polynomial, dE/dD into GD, dE/dA into W;
//   3. forces   (a block per atom): its threads over the neighbours, radial and Y_lm derivatives contracted with W.
// Input: the CSR bonds (mlip_flat_bonds.h, r_ij = r_j - r_i). Output: forces on all atoms (ghosts included), per-atom
// energies, per-atom virials (the whole bond virial on the central atom).

#include <exaStamp/mlip_polymlp/polymlp_tables.h>
#include <exaStamp/mlip_utils/mlip_flat_bonds.h>
#include <exaStamp/mlip_utils/mlip_ylm.h>

#include <onika/cuda/cuda.h>
#include <onika/parallel/parallel_for.h>
#include <onika/parallel/block_parallel_for.h>

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace exaStamp
{
  namespace mlip_polymlp
  {
    namespace polymlp_onika
    {
      using mlip_utils::BondsView;
      using mlip_utils::OutputsView;

      // per-chunk arrays [entry][atom], blocked by groups of V atoms: element (q, a) at (a >> vs) * group + q * V + (a & vm).
      // GPU: one group of the whole chunk (atom fastest); host: groups of kPolymlpLanes atoms (SIMD lanes of the readout)
      struct ChunkView
      {
        unsigned vs = 62;
        size_t vm = ~size_t(0), V = 1;
        size_t sg = 0, fg = 0, cg = 0;   // group strides: A / W, D / GD, C / G
        ONIKA_HOST_DEVICE_FUNC inline size_t ia(int q, size_t a) const { return ( a >> vs ) * sg + q * V + ( a & vm ); }
        ONIKA_HOST_DEVICE_FUNC inline size_t id(int k, size_t a) const { return ( a >> vs ) * fg + k * V + ( a & vm ); }
        ONIKA_HOST_DEVICE_FUNC inline size_t ic(int c, size_t a) const { return ( a >> vs ) * cg + c * V + ( a & vm ); }
        const int* order = nullptr;
        double *Are = nullptr, *Aim = nullptr, *Wre = nullptr, *Wim = nullptr, *D = nullptr, *GD = nullptr;
        double *Cre = nullptr, *Cim = nullptr, *Gre = nullptr, *Gim = nullptr;
        double* scratch = nullptr;
        size_t nscratch = 0;
      };

      ONIKA_HOST_DEVICE_FUNC inline size_t block_slot(size_t a)
      {
        return ONIKA_CU_VALUE_IF_CUDA( a , size_t( ONIKA_CU_BLOCK_IDX ) );
      }

      ONIKA_HOST_DEVICE_FUNC inline size_t basis_record(const PolymlpView& t) { return size_t(t.ng) + 2 * size_t(t.nlh) + 1; }
      ONIKA_HOST_DEVICE_FUNC inline size_t forces_record(const PolymlpView& t) { return 2 * size_t(t.ng) + 8 * size_t(t.nlh); }

      // gaussians f_n(r) (and derivatives)
      template<bool DER>
      ONIKA_HOST_DEVICE_FUNC inline void gaussians( const PolymlpView& t, double r, double* f, double* df )
      {
        const double pi = M_PI;
        const double fc = 0.5 * ( 1.0 + cos( pi * r / t.rc ) );
        const double dfc = -0.5 * pi * sin( pi * r / t.rc ) / t.rc;
        for(int n=0;n<t.ng;n++)
        {
          const double x = r - t.mu[n];
          const double e = exp( -t.beta[n] * x * x );
          f[n] = e * fc;
          if( DER ) df[n] = e * ( dfc - 2.0 * t.beta[n] * x * fc );
        }
      }

      struct BasisFunctor
      {
        PolymlpView t;
        BondsView bv;
        ChunkView ch;
        int c0 = 0;

        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t a ) const
        {
          const int B = ONIKA_CU_BLOCK_SIZE, tid = ONIKA_CU_THREAD_IDX;
          const int i = ch.order[c0+a], ti = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[ti], ns = t.src0[ti+1] - s0;
          const size_t PT = basis_record( t );
          const double rc2 = t.rc * t.rc;
          double* scr = ch.scratch + block_slot(a) * ch.nscratch;

          for(int q=tid;q<ns;q+=B) { ch.Are[ch.ia(q,a)] = 0.0; ch.Aim[ch.ia(q,a)] = 0.0; }

          for(int k0=0;k0<jnum;k0+=B)
          {
            ONIKA_CU_BLOCK_SYNC();
            {
              double* f = scr + size_t(tid) * PT;
              double* yr = f + t.ng;
              double* yi = yr + t.nlh;
              double& flag = yi[t.nlh];
              flag = -1.0;
              const int jj = k0 + tid;
              if( jj < jnum )
              {
                const int j = bv.j[b0+jj], nu = bv.type[j];
                const double* d = bv.rij + 3 * ( b0 + jj );
                const double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                if( nu >= 0 && r2 < rc2 && r2 > 0.0 )
                {
                  const double r = sqrt( r2 ), inv_r = 1.0 / r;
                  gaussians<false>( t, r, f, nullptr );
                  mlip_utils::ylm<false>( t.lmax, t.ylm_a, t.ylm_b, d[0]*inv_r, d[1]*inv_r, d[2]*inv_r, inv_r, yr, yi, nullptr, nullptr );
                  flag = nu;
                }
              }
            }
            ONIKA_CU_BLOCK_SYNC();
            const int nb = jnum - k0 < B ? jnum - k0 : B;
            if( B == 1 )
            {
              const double* f = scr;
              const double* yr = f + t.ng;
              const double* yi = yr + t.nlh;
              const int nu = int( yi[t.nlh] );
              if( nu < 0 ) continue;
              const int qa = t.srcb[ ti * ( t.nel + 1 ) + nu ], qb = t.srcb[ ti * ( t.nel + 1 ) + nu + 1 ];
              for(int q=qa;q<qb;q++)
              {
                const int sq = s0 + q, n = t.src_n[sq], lm = t.src_lm[sq];
                ch.Are[ch.ia(q,a)] += f[n] * yr[lm];
                ch.Aim[ch.ia(q,a)] += f[n] * yi[lm];
              }
              continue;
            }
            for(int q=tid;q<ns;q+=B)
            {
              const int sq = s0 + q, nu = t.src_nu[sq], n = t.src_n[sq], lm = t.src_lm[sq];
              double sr = 0.0, si = 0.0;
              for(int tt=0;tt<nb;tt++)
              {
                const double* f = scr + size_t(tt) * PT;
                const double* yr = f + t.ng;
                const double* yi = yr + t.nlh;
                if( int( yi[t.nlh] ) != nu ) continue;
                sr += f[n] * yr[lm];
                si += f[n] * yi[lm];
              }
              ch.Are[ch.ia(q,a)] += sr;
              ch.Aim[ch.ia(q,a)] += si;
            }
          }
        }
      };

      struct ReadoutFunctor
      {
        PolymlpView t;
        BondsView bv;
        ChunkView ch;
        OutputsView out;
        int c0 = 0;

        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t a ) const
        {
          const int i = ch.order[c0+a], ti = bv.type[i];
          const int ns = t.src0[ti+1] - t.src0[ti], h0 = t.ch0[ti], nch = t.ch0[ti+1] - h0, nf = t.ft0[ti+1] - t.ft0[ti];
          const int n0 = t.no0[ti], n1 = t.no0[ti+1];

          // channel values
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c];
            const double sc = t.ch_scale[h0+c];
            ch.Cre[ch.ic(c,a)] = sc * ch.Are[ch.ia(q,a)];
            ch.Cim[ch.ic(c,a)] = t.ch_conj[h0+c] ? -sc * ch.Aim[ch.ia(q,a)] : sc * ch.Aim[ch.ia(q,a)];
            ch.Gre[ch.ic(c,a)] = 0.0;
            ch.Gim[ch.ic(c,a)] = 0.0;
          }
          for(int k=0;k<nf;k++) { ch.D[ch.id(k,a)] = 0.0; ch.GD[ch.id(k,a)] = 0.0; }

          // features: depth-first walk of the product tree, value of the path at each depth
          double pr[kPolymlpMaxFactors+1], pi[kPolymlpMaxFactors+1];
          pr[0] = 1.0; pi[0] = 0.0;
          for(int n=n0;n<n1;n++)
          {
            const int d = t.node_depth[n], c = t.node_ch[n];
            const double cr = ch.Cre[ch.ic(c,a)], ci = ch.Cim[ch.ic(c,a)];
            pr[d] = pr[d-1] * cr - pi[d-1] * ci;
            pi[d] = pr[d-1] * ci + pi[d-1] * cr;
            for(int64_t s=t.nto[n];s<t.nto[n+1];s++) ch.D[ch.id(t.nt_feat[s],a)] += t.nt_coef[s] * pr[d];
          }

          // polynomial: energy, dE/dD
          double e = 0.0;
          for(int u=t.po0[ti];u<t.po0[ti+1];u++)
          {
            const int64_t f0 = t.pfo[u];
            const int np = int( t.pfo[u+1] - f0 );
            const double w = t.pcoef[u];
            double v[kPolymlpMaxFactors], suf[kPolymlpMaxFactors+1];
            for(int f=0;f<np;f++) v[f] = ch.D[ch.id(t.pfac[f0+f],a)];
            suf[np] = 1.0;
            for(int f=np-1;f>=0;f--) suf[f] = suf[f+1] * v[f];
            e += w * suf[0];
            double pre = w;
            for(int f=0;f<np;f++) { ch.GD[ch.id(t.pfac[f0+f],a)] += pre * suf[f+1]; pre *= v[f]; }
          }
          out.eatom[i] += e;

          // reverse mode along the tree: adjoint of each node of the open path (d(E)/d(path value), Re linear), closed
          // when the walk leaves it: G(channel) += adj * parent value, adj(parent) += adj * channel
          double ar[kPolymlpMaxFactors+1], ai[kPolymlpMaxFactors+1];
          int pc[kPolymlpMaxFactors+1];
          int top = 0;
          auto close = [&](int k)
          {
            const int c = pc[k];
            const double cr = ch.Cre[ch.ic(c,a)], ci = ch.Cim[ch.ic(c,a)];
            ch.Gre[ch.ic(c,a)] += ar[k] * pr[k-1] - ai[k] * pi[k-1];
            ch.Gim[ch.ic(c,a)] += ar[k] * pi[k-1] + ai[k] * pr[k-1];
            ar[k-1] += ar[k] * cr - ai[k] * ci;
            ai[k-1] += ar[k] * ci + ai[k] * cr;
          };
          ar[0] = 0.0; ai[0] = 0.0;
          for(int n=n0;n<n1;n++)
          {
            const int d = t.node_depth[n], c = t.node_ch[n];
            while( top >= d ) { close( top ); -- top; }
            const double cr = ch.Cre[ch.ic(c,a)], ci = ch.Cim[ch.ic(c,a)];
            pr[d] = pr[d-1] * cr - pi[d-1] * ci;
            pi[d] = pr[d-1] * ci + pi[d-1] * cr;
            double g = 0.0;
            for(int64_t s=t.nto[n];s<t.nto[n+1];s++) g += ch.GD[ch.id(t.nt_feat[s],a)] * t.nt_coef[s];
            ar[d] = g; ai[d] = 0.0;
            pc[d] = c;
            top = d;
          }
          while( top >= 1 ) { close( top ); -- top; }

          // channel gradients into dE/dA: channel = s A (Re: dA_re s G_re, dA_im -s G_im) or s conj(A) (dA_im +s G_im)
          for(int q=0;q<ns;q++) { ch.Wre[ch.ia(q,a)] = 0.0; ch.Wim[ch.ia(q,a)] = 0.0; }
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c];
            const double sc = t.ch_scale[h0+c], Gr = ch.Gre[ch.ic(c,a)], Gi = ch.Gim[ch.ic(c,a)];
            ch.Wre[ch.ia(q,a)] += sc * Gr;
            ch.Wim[ch.ia(q,a)] += t.ch_conj[h0+c] ? sc * Gi : -sc * Gi;
          }
        }
      };

      struct ForcesFunctor
      {
        PolymlpView t;
        BondsView bv;
        ChunkView ch;
        OutputsView out;
        int c0 = 0;

        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t a ) const
        {
          const int B = ONIKA_CU_BLOCK_SIZE, tid = ONIKA_CU_THREAD_IDX;
          const int i = ch.order[c0+a], ti = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[ti];
          const double rc2 = t.rc * t.rc;
          double* rec = ch.scratch + block_slot(a) * ch.nscratch + size_t(tid) * forces_record( t );
          double* f = rec;
          double* df = f + t.ng;
          double* yr = df + t.ng;
          double* yi = yr + t.nlh;
          double* gr = yi + t.nlh;
          double* gi = gr + 3 * t.nlh;
          double fi[3] = { 0.0, 0.0, 0.0 }, vi[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          for(int jj=tid;jj<jnum;jj+=B)
          {
            const int j = bv.j[b0+jj], nu = bv.type[j];
            if( nu < 0 ) continue;
            const double* d = bv.rij + 3 * ( b0 + jj );
            const double dx = d[0], dy = d[1], dz = d[2], r2 = dx*dx + dy*dy + dz*dz;
            if( r2 >= rc2 || r2 == 0.0 ) continue;
            const double r = sqrt( r2 ), inv_r = 1.0 / r, ux = dx * inv_r, uy = dy * inv_r, uz = dz * inv_r;
            gaussians<true>( t, r, f, df );
            mlip_utils::ylm<true>( t.lmax, t.ylm_a, t.ylm_b, ux, uy, uz, inv_r, yr, yi, gr, gi );
            double rad = 0.0, ang[3] = { 0.0, 0.0, 0.0 };
            const int qa = t.srcb[ ti * ( t.nel + 1 ) + nu ], qb = t.srcb[ ti * ( t.nel + 1 ) + nu + 1 ];
            for(int q=qa;q<qb;q++)
            {
              const int sq = s0 + q, n = t.src_n[sq], lm = t.src_lm[sq];
              const double wr = ch.Wre[ch.ia(q,a)], wi = ch.Wim[ch.ia(q,a)];
              rad += df[n] * ( wr * yr[lm] + wi * yi[lm] );
              if( lm == 0 ) continue;
              const double fv = f[n];
              for(int k=0;k<3;k++) ang[k] += fv * ( wr * gr[3*lm+k] + wi * gi[3*lm+k] );
            }
            const double fv[3] = { rad * ux + ang[0], rad * uy + ang[1], rad * uz + ang[2] };    // dE_i / d r_ij
            for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(j) + k ] , -fv[k] );
            for(int k=0;k<3;k++) fi[k] += fv[k];
            vi[0] -= dx * fv[0]; vi[1] -= dy * fv[1]; vi[2] -= dz * fv[2];
            vi[3] -= dx * fv[1]; vi[4] -= dx * fv[2]; vi[5] -= dy * fv[2];
          }
          for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(i) + k ] , fi[k] );
          if( out.vatom != nullptr ) for(int k=0;k<6;k++) ONIKA_CU_ATOMIC_ADD( out.vatom[ 6 * size_t(i) + k ] , vi[k] );
        }
      };

      static constexpr int kPolymlpLanes = 8;

      // host readout of a group of kPolymlpLanes atoms of one type, vectorised over the atoms (scalar readout for a
      // group mixing types or a partial group)
      struct ReadoutBatchFunctor
      {
        ReadoutFunctor one;
        int nc = 0;

        inline void operator () ( size_t g ) const
        {
          constexpr int L = kPolymlpLanes;
          const PolymlpView& t = one.t;
          const ChunkView& ch = one.ch;
          const size_t a0 = g * L;
          const int nb = int( nc - a0 < size_t(L) ? nc - a0 : L );
          const int c0 = one.c0;
          const int ti = one.bv.type[ ch.order[c0+a0] ];
          bool same = nb == L;
          for(int b=1;b<nb && same;b++) same = one.bv.type[ ch.order[c0+a0+b] ] == ti;
          if( !same ) { for(int b=0;b<nb;b++) one( a0 + b ); return; }

          const int ns = t.src0[ti+1] - t.src0[ti], h0 = t.ch0[ti], nch = t.ch0[ti+1] - h0, nf = t.ft0[ti+1] - t.ft0[ti];
          const int n0 = t.no0[ti], n1 = t.no0[ti+1];
          double* __restrict__ Cre = ch.Cre + ch.ic(0,a0);
          double* __restrict__ Cim = ch.Cim + ch.ic(0,a0);
          double* __restrict__ Gre = ch.Gre + ch.ic(0,a0);
          double* __restrict__ Gim = ch.Gim + ch.ic(0,a0);
          const double* __restrict__ Are = ch.Are + ch.ia(0,a0);
          const double* __restrict__ Aim = ch.Aim + ch.ia(0,a0);
          double* __restrict__ Wre = ch.Wre + ch.ia(0,a0);
          double* __restrict__ Wim = ch.Wim + ch.ia(0,a0);
          double* __restrict__ D = ch.D + ch.id(0,a0);
          double* __restrict__ GD = ch.GD + ch.id(0,a0);
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c];
            const double sr = t.ch_scale[h0+c], si = t.ch_conj[h0+c] ? -sr : sr;
#           pragma omp simd
            for(int b=0;b<L;b++) { Cre[c*L+b] = sr * Are[q*L+b]; Cim[c*L+b] = si * Aim[q*L+b]; Gre[c*L+b] = 0.0; Gim[c*L+b] = 0.0; }
          }
          for(int k=0;k<nf*L;k++) { D[k] = 0.0; GD[k] = 0.0; }
          double pr[kPolymlpMaxFactors+1][L], pi[kPolymlpMaxFactors+1][L], ar[kPolymlpMaxFactors+1][L], ai[kPolymlpMaxFactors+1][L];
          for(int b=0;b<L;b++) { pr[0][b] = 1.0; pi[0][b] = 0.0; ar[0][b] = 0.0; ai[0][b] = 0.0; }
          for(int n=n0;n<n1;n++)
          {
            const int d = t.node_depth[n], c = t.node_ch[n];
#           pragma omp simd
            for(int b=0;b<L;b++)
            {
              const double cr = Cre[c*L+b], ci = Cim[c*L+b];
              pr[d][b] = pr[d-1][b] * cr - pi[d-1][b] * ci;
              pi[d][b] = pr[d-1][b] * ci + pi[d-1][b] * cr;
            }
            for(int64_t s=t.nto[n];s<t.nto[n+1];s++)
            {
              const int k = t.nt_feat[s];
              const double cf = t.nt_coef[s];
#             pragma omp simd
              for(int b=0;b<L;b++) D[k*L+b] += cf * pr[d][b];
            }
          }
          double e[L];
          for(int b=0;b<L;b++) e[b] = 0.0;
          for(int u=t.po0[ti];u<t.po0[ti+1];u++)
          {
            const int64_t f0 = t.pfo[u];
            const int np = int( t.pfo[u+1] - f0 );
            const double w = t.pcoef[u];
            if( np == 1 )
            {
              const int k0 = t.pfac[f0];
#             pragma omp simd
              for(int b=0;b<L;b++) { e[b] += w * D[k0*L+b]; GD[k0*L+b] += w; }
            }
            else if( np == 2 )
            {
              const int k0 = t.pfac[f0], k1 = t.pfac[f0+1];
#             pragma omp simd
              for(int b=0;b<L;b++)
              {
                const double v0 = D[k0*L+b], v1 = D[k1*L+b];
                e[b] += w * v0 * v1; GD[k0*L+b] += w * v1; GD[k1*L+b] += w * v0;
              }
            }
            else
            {
              for(int b=0;b<L;b++)
              {
                double v[kPolymlpMaxFactors], suf[kPolymlpMaxFactors+1];
                for(int f=0;f<np;f++) v[f] = D[t.pfac[f0+f]*L+b];
                suf[np] = 1.0;
                for(int f=np-1;f>=0;f--) suf[f] = suf[f+1] * v[f];
                e[b] += w * suf[0];
                double pre = w;
                for(int f=0;f<np;f++) { GD[t.pfac[f0+f]*L+b] += pre * suf[f+1]; pre *= v[f]; }
              }
            }
          }
          int pc[kPolymlpMaxFactors+1];
          int top = 0;
          auto close = [&](int k)
          {
            const int c = pc[k];
#           pragma omp simd
            for(int b=0;b<L;b++)
            {
              const double cr = Cre[c*L+b], ci = Cim[c*L+b];
              Gre[c*L+b] += ar[k][b] * pr[k-1][b] - ai[k][b] * pi[k-1][b];
              Gim[c*L+b] += ar[k][b] * pi[k-1][b] + ai[k][b] * pr[k-1][b];
              ar[k-1][b] += ar[k][b] * cr - ai[k][b] * ci;
              ai[k-1][b] += ar[k][b] * ci + ai[k][b] * cr;
            }
          };
          for(int n=n0;n<n1;n++)
          {
            const int d = t.node_depth[n], c = t.node_ch[n];
            while( top >= d ) { close( top ); -- top; }
            double g[L];
            for(int b=0;b<L;b++) g[b] = 0.0;
            for(int64_t s=t.nto[n];s<t.nto[n+1];s++)
            {
              const int k = t.nt_feat[s];
              const double cf = t.nt_coef[s];
#             pragma omp simd
              for(int b=0;b<L;b++) g[b] += GD[k*L+b] * cf;
            }
#           pragma omp simd
            for(int b=0;b<L;b++)
            {
              const double cr = Cre[c*L+b], ci = Cim[c*L+b];
              pr[d][b] = pr[d-1][b] * cr - pi[d-1][b] * ci;
              pi[d][b] = pr[d-1][b] * ci + pi[d-1][b] * cr;
              ar[d][b] = g[b]; ai[d][b] = 0.0;
            }
            pc[d] = c;
            top = d;
          }
          while( top >= 1 ) { close( top ); -- top; }
          for(int k=0;k<ns*L;k++) { Wre[k] = 0.0; Wim[k] = 0.0; }
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c];
            const double sc = t.ch_scale[h0+c], si = t.ch_conj[h0+c] ? sc : -sc;
#           pragma omp simd
            for(int b=0;b<L;b++) { Wre[q*L+b] += sc * Gre[c*L+b]; Wim[q*L+b] += si * Gim[c*L+b]; }
          }
          for(int b=0;b<L;b++) one.out.eatom[ ch.order[c0+a0+b] ] += e[b];
        }
      };

      // host basis and forces of one atom: neighbours of each type in groups of kPolymlpLanes, gaussians and Y_lm
      // vectorised over the lanes. Thread scratch: neighbour list, lane arrays.
      inline size_t host_lane_scratch(const PolymlpView& t, size_t jmax)
      {
        return ( jmax + 1 ) / 2 + 1 + size_t(kPolymlpLanes) * ( 2 * size_t(t.ng) + 8 * size_t(t.nlh) + 8 );
      }

      template<bool DER>
      struct HostLanes
      {
        static constexpr int L = kPolymlpLanes;
        int* jl;
        double *f, *df, *yr, *yi, *gr, *gi, *ux, *uy, *uz, *ir, *rr, *dx, *dy, *dz;
        HostLanes(const PolymlpView& t, double* scr, size_t jmax)
        {
          jl = reinterpret_cast<int*>( scr );
          double* p = scr + ( jmax + 1 ) / 2 + 1;
          auto take = [&](size_t n) { double* x = p; p += n; return x; };
          f = take( L * t.ng ); df = take( L * t.ng );
          yr = take( L * t.nlh ); yi = take( L * t.nlh ); gr = take( 3 * L * t.nlh ); gi = take( 3 * L * t.nlh );
          ux = take( L ); uy = take( L ); uz = take( L ); ir = take( L ); rr = take( L ); dx = take( L ); dy = take( L ); dz = take( L );
        }
        static int list( const PolymlpView& t, const BondsView& bv, int64_t b0, int jnum, int nu, int* jl )
        {
          int n = 0;
          const double rc2 = t.rc * t.rc;
          for(int jj=0;jj<jnum;jj++)
          {
            if( bv.type[ bv.j[b0+jj] ] != nu ) continue;
            const double* d = bv.rij + 3 * ( b0 + jj );
            const double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            if( r2 < rc2 && r2 > 0.0 ) jl[n++] = jj;
          }
          return n;
        }
        void fill( const PolymlpView& t, const BondsView& bv, int64_t b0, int k0, int nb )
        {
          for(int x=0;x<L;x++)
          {
            if( x >= nb ) { ux[x] = 0.0; uy[x] = 0.0; uz[x] = 1.0; ir[x] = 1.0; rr[x] = 2.0 * t.rc; dx[x] = dy[x] = dz[x] = 0.0; continue; }
            const double* d = bv.rij + 3 * ( b0 + jl[k0+x] );
            const double r = sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
            rr[x] = r; ir[x] = 1.0 / r; dx[x] = d[0]; dy[x] = d[1]; dz[x] = d[2];
            ux[x] = d[0] * ir[x]; uy[x] = d[1] * ir[x]; uz[x] = d[2] * ir[x];
          }
          // gaussians times the cosine cutoff, zero for unused lanes (r = 2 rc)
          const double pi = M_PI;
          double fc[L], dfc[L];
#         pragma omp simd
          for(int x=0;x<L;x++)
          {
            const bool in = rr[x] < t.rc;
            fc[x] = in ? 0.5 * ( 1.0 + cos( pi * rr[x] / t.rc ) ) : 0.0;
            dfc[x] = in ? -0.5 * pi * sin( pi * rr[x] / t.rc ) / t.rc : 0.0;
          }
          for(int n=0;n<t.ng;n++)
          {
            const double be = t.beta[n], m = t.mu[n];
#           pragma omp simd
            for(int x=0;x<L;x++)
            {
              const double y = rr[x] - m;
              const double e = exp( -be * y * y );
              f[n*L+x] = e * fc[x];
              if( DER ) df[n*L+x] = e * ( dfc[x] - 2.0 * be * y * fc[x] );
            }
          }
          mlip_utils::ylm_lanes<L,DER>( t.lmax, t.ylm_a, t.ylm_b, ux, uy, uz, ir, yr, yi, gr, gi );
        }
      };

      struct BasisHostFunctor
      {
        PolymlpView t;
        BondsView bv;
        ChunkView ch;
        int c0 = 0;
        size_t jmax = 0;

        inline void operator () ( size_t a ) const
        {
          constexpr int L = kPolymlpLanes;
          const int i = ch.order[c0+a], ti = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[ti], ns = t.src0[ti+1] - s0;
          HostLanes<false> h( t, ch.scratch + block_slot(a) * ch.nscratch, jmax );
          for(int q=0;q<ns;q++) { ch.Are[ch.ia(q,a)] = 0.0; ch.Aim[ch.ia(q,a)] = 0.0; }
          for(int nu=0;nu<t.nel;nu++)
          {
            const int qa = t.srcb[ ti * ( t.nel + 1 ) + nu ], qb = t.srcb[ ti * ( t.nel + 1 ) + nu + 1 ];
            if( qa == qb ) continue;
            const int cnt = HostLanes<false>::list( t, bv, b0, jnum, nu, h.jl );
            for(int k0=0;k0<cnt;k0+=L)
            {
              h.fill( t, bv, b0, k0, cnt - k0 < L ? cnt - k0 : L );
              for(int q=qa;q<qb;q++)
              {
                const int sq = s0 + q, n = t.src_n[sq], lm = t.src_lm[sq];
                double sr = 0.0, si = 0.0;
#               pragma omp simd reduction(+:sr,si)
                for(int x=0;x<L;x++) { sr += h.f[n*L+x] * h.yr[lm*L+x]; si += h.f[n*L+x] * h.yi[lm*L+x]; }
                ch.Are[ch.ia(q,a)] += sr;
                ch.Aim[ch.ia(q,a)] += si;
              }
            }
          }
        }
      };

      struct ForcesHostFunctor
      {
        PolymlpView t;
        BondsView bv;
        ChunkView ch;
        OutputsView out;
        int c0 = 0;
        size_t jmax = 0;

        inline void operator () ( size_t a ) const
        {
          constexpr int L = kPolymlpLanes;
          const int i = ch.order[c0+a], ti = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[ti];
          HostLanes<true> h( t, ch.scratch + block_slot(a) * ch.nscratch, jmax );
          double fi[3] = { 0.0, 0.0, 0.0 }, vi[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          for(int nu=0;nu<t.nel;nu++)
          {
            const int qa = t.srcb[ ti * ( t.nel + 1 ) + nu ], qb = t.srcb[ ti * ( t.nel + 1 ) + nu + 1 ];
            if( qa == qb ) continue;
            const int cnt = HostLanes<true>::list( t, bv, b0, jnum, nu, h.jl );
            for(int k0=0;k0<cnt;k0+=L)
            {
              const int nb = cnt - k0 < L ? cnt - k0 : L;
              h.fill( t, bv, b0, k0, nb );
              double rad[L], ax[L], ay[L], az[L];
              for(int x=0;x<L;x++) { rad[x] = 0.0; ax[x] = 0.0; ay[x] = 0.0; az[x] = 0.0; }
              for(int q=qa;q<qb;q++)
              {
                const int sq = s0 + q, n = t.src_n[sq], lm = t.src_lm[sq];
                const double wr = ch.Wre[ch.ia(q,a)], wi = ch.Wim[ch.ia(q,a)];
                const double* __restrict__ gxr = h.gr + 3*lm*L; const double* __restrict__ gxi = h.gi + 3*lm*L;
#               pragma omp simd
                for(int x=0;x<L;x++)
                {
                  rad[x] += h.df[n*L+x] * ( wr * h.yr[lm*L+x] + wi * h.yi[lm*L+x] );
                  const double fv = h.f[n*L+x];
                  ax[x] += fv * ( wr * gxr[x] + wi * gxi[x] );
                  ay[x] += fv * ( wr * gxr[L+x] + wi * gxi[L+x] );
                  az[x] += fv * ( wr * gxr[2*L+x] + wi * gxi[2*L+x] );
                }
              }
              for(int x=0;x<nb;x++)
              {
                const int j = bv.j[ b0 + h.jl[k0+x] ];
                const double fv[3] = { rad[x] * h.ux[x] + ax[x], rad[x] * h.uy[x] + ay[x], rad[x] * h.uz[x] + az[x] };
                for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(j) + k ] , -fv[k] );
                for(int k=0;k<3;k++) fi[k] += fv[k];
                const double ddx = h.dx[x], ddy = h.dy[x], ddz = h.dz[x];
                vi[0] -= ddx * fv[0]; vi[1] -= ddy * fv[1]; vi[2] -= ddz * fv[2];
                vi[3] -= ddx * fv[1]; vi[4] -= ddx * fv[2]; vi[5] -= ddy * fv[2];
              }
            }
          }
          for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(i) + k ] , fi[k] );
          if( out.vatom != nullptr ) for(int k=0;k<6;k++) ONIKA_CU_ATOMIC_ADD( out.vatom[ 6 * size_t(i) + k ] , vi[k] );
        }
      };

    }
  }
}

namespace onika
{
  namespace parallel
  {
    template<> struct BlockParallelForFunctorTraits< exaStamp::mlip_polymlp::polymlp_onika::BasisFunctor > { static inline constexpr bool CudaCompatible = true; };
    template<> struct BlockParallelForFunctorTraits< exaStamp::mlip_polymlp::polymlp_onika::ForcesFunctor > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::mlip_polymlp::polymlp_onika::ReadoutFunctor > { static inline constexpr bool CudaCompatible = true; };
  }
}
