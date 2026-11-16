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

// YE3T (.yace, linear ACE) onika kernels, CPU and GPU, fp64. Per chunk of central atoms ordered by element:
//   1. basis    (a block per atom): its threads over a batch of neighbours (radial functions, Y_lm into the block
//                scratch), then over the sources (sums over the batch, in a fixed order) -> A [source][atom];
//   2. readout  (a thread per atom): rho and E_i along the monomial product tree, dE/dA into W [source][atom] by a
//                reverse walk of the tree;
//   3. forces   (a block per atom): its threads over the neighbours, radial and Y_lm derivatives contracted with W.
// Input: the CSR bonds (mlip_flat_bonds.h, r_ij = r_j - r_i). Output: forces on all atoms (ghosts included), per-atom
// energies, per-atom virials (the whole bond virial on the central atom). Block scratch: per chunk atom on a GPU, per
// OpenMP thread on the host (block size 1).

#include <exaStamp/mlip_ye3t/ye3t_math.h>
#include <exaStamp/mlip_utils/mlip_flat_bonds.h>

#include <onika/cuda/cuda.h>
#include <onika/parallel/parallel_for.h>
#include <onika/parallel/block_parallel_for.h>

#include <cstddef>
#include <cstdint>

namespace exaStamp
{
  namespace mlip_ye3t
  {
    namespace ye3t_onika
    {
      using mlip_utils::BondsView;
      using mlip_utils::OutputsView;

      // per-chunk arrays [source][atom of the chunk], blocked by groups of V atoms: element (q, a) at
      // (a >> vs) * sg + q * V + (a & vm). GPU: one group of the whole chunk (atom fastest, V = chunk, sg = 0);
      // host: groups of 8 atoms (the SIMD lanes of the host readout), sg = entries * 8
      struct ChunkView
      {
        unsigned vs = 62;             // log2 V (GPU: no group)
        size_t vm = ~size_t(0), V = 1;
        size_t sg = 0, cg = 0;        // group strides of the source / channel arrays
        ONIKA_HOST_DEVICE_FUNC inline size_t ia(int q, size_t a) const { return ( a >> vs ) * sg + q * V + ( a & vm ); }
        ONIKA_HOST_DEVICE_FUNC inline size_t ic(int c, size_t a) const { return ( a >> vs ) * cg + c * V + ( a & vm ); }
        const int* order = nullptr;   // chunk atom -> flat atom
        double *Are = nullptr, *Aim = nullptr, *Wre = nullptr, *Wim = nullptr;
        double *Cre = nullptr, *Cim = nullptr, *Gre = nullptr, *Gim = nullptr;   // channel values, d(fs rho)/d(channel)
        double* scratch = nullptr;    // [slot][nscratch]
        size_t nscratch = 0;
        int* flags = nullptr;         // [1] density guard violations
      };

      ONIKA_HOST_DEVICE_FUNC inline size_t block_slot(size_t a)
      {
        return ONIKA_CU_VALUE_IF_CUDA( a , size_t( ONIKA_CU_BLOCK_IDX ) );
      }

      // doubles of one neighbour record: basis (values), forces (values, derivatives, Y gradients)
      ONIKA_HOST_DEVICE_FUNC inline size_t basis_record(const Ye3tView& t) { return size_t(t.nbase_max) + t.ncontr_max + 2 * size_t(t.nlh) + 1; }
      ONIKA_HOST_DEVICE_FUNC inline size_t forces_record(const Ye3tView& t) { return 2 * size_t(t.nbase_max) + 2 * size_t(t.ncontr_max) + 8 * size_t(t.nlh); }

      // radial values (and derivatives) of a bond at r: base g_k into g (dg) from their spline (or directly below the
      // table / in analytic mode), contracted R_nl = sum_k c g_k into R (dR). The cubic Hermite spline is linear in its
      // node data: this is the spline of R_nl with the contracted nodes of PACE / ye3t-lammps, to round-off
      template<bool DER>
      ONIKA_HOST_DEVICE_FUNC inline void radials( const Ye3tView& t, const Ye3tBond& b, double r, double* g, double* dg, double* R, double* dR )
      {
        SplinePos sp;
        double dgl[64];
        double* dd = DER ? dg : dgl;
        if( !t.analytic && sp.locate( b, r ) )
        {
          const double* bs = t.pool + b.bspl;
          double v, d;
          for(int k=0;k<b.nbase;k++) { sp.eval( bs, b.nbase, k, v, d ); g[k] = v; dd[k] = d; }
        }
        else
        {
          if( b.nbase > 64 ) return;    // refused at init
          cheb_exp_cos( b, r, g, dd );
        }
        const double* c = t.pool + b.crad;
        for(int q=0;q<b.ncontr;q++)
        {
          double s = 0.0, ds = 0.0;
          for(int k=0;k<b.nbase;k++) { s += c[size_t(q)*b.nbase+k] * g[k]; ds += c[size_t(q)*b.nbase+k] * dd[k]; }
          R[q] = s;
          if( DER ) dR[q] = ds;
        }
      }

      struct BasisFunctor
      {
        Ye3tView t;
        BondsView bv;
        ChunkView ch;
        int c0 = 0;

        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t a ) const
        {
          const int B = ONIKA_CU_BLOCK_SIZE, tid = ONIKA_CU_THREAD_IDX;
          const int i = ch.order[c0+a], mu = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[mu], ns = t.src0[mu+1] - s0;
          const size_t PT = basis_record( t );
          double* scr = ch.scratch + block_slot(a) * ch.nscratch;

          for(int q=tid;q<ns;q+=B) { ch.Are[ch.ia(q,a)] = 0.0; ch.Aim[ch.ia(q,a)] = 0.0; }

          for(int k0=0;k0<jnum;k0+=B)
          {
            ONIKA_CU_BLOCK_SYNC();
            {
              double* rec = scr + size_t(tid) * PT;
              double* g = rec;
              double* R = g + t.nbase_max;
              double* yr = R + t.ncontr_max;
              double* yi = yr + t.nlh;
              double& flag = yi[t.nlh];
              flag = -1.0;
              const int jj = k0 + tid;
              if( jj < jnum )
              {
                const int j = bv.j[b0+jj], nu = bv.type[j];
                const Ye3tBond& b = t.bonds[ mu * t.nel + nu ];
                const double* d = bv.rij + 3 * ( b0 + jj );
                const double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                if( nu >= 0 && r2 < b.rc * b.rc && r2 > 0.0 )
                {
                  const double r = sqrt( r2 ), inv_r = 1.0 / r;
                  radials<false>( t, b, r, g, nullptr, R, nullptr );
                  ylm<false>( b.lb, t.ylm_a, t.ylm_b, d[0]*inv_r, d[1]*inv_r, d[2]*inv_r, inv_r, yr, yi, nullptr, nullptr );
                  flag = nu;
                }
              }
            }
            ONIKA_CU_BLOCK_SYNC();
            const int nb = jnum - k0 < B ? jnum - k0 : B;
            if( B == 1 )
            {
              // host: the sources of the neighbour's element only
              const double* g = scr;
              const double* R = g + t.nbase_max;
              const double* yr = R + t.ncontr_max;
              const double* yi = yr + t.nlh;
              const int nu = int( yi[t.nlh] );
              if( nu < 0 ) continue;
              const int qa = t.srcb[ mu * ( t.nel + 1 ) + nu ], qb = t.srcb[ mu * ( t.nel + 1 ) + nu + 1 ];
              for(int q=qa;q<qb;q++)
              {
                const int sq = s0 + q, rr = t.src_r[sq];
                if( t.src_kind[sq] == 0 ) ch.Are[ch.ia(q,a)] += g[rr];
                else
                {
                  const int lm = t.src_lm[sq];
                  ch.Are[ch.ia(q,a)] += R[rr] * yr[lm];
                  ch.Aim[ch.ia(q,a)] += R[rr] * yi[lm];
                }
              }
              continue;
            }
            for(int q=tid;q<ns;q+=B)
            {
              const int sq = s0 + q, nu = t.src_nu[sq], rr = t.src_r[sq], lm = t.src_lm[sq];
              const bool base = t.src_kind[sq] == 0;
              double sr = 0.0, si = 0.0;
              for(int tt=0;tt<nb;tt++)
              {
                const double* g = scr + size_t(tt) * PT;
                const double* R = g + t.nbase_max;
                const double* yr = R + t.ncontr_max;
                const double* yi = yr + t.nlh;
                if( int( yi[t.nlh] ) != nu ) continue;
                if( base ) sr += g[rr];
                else { sr += R[rr] * yr[lm]; si += R[rr] * yi[lm]; }
              }
              ch.Are[ch.ia(q,a)] += sr;
              ch.Aim[ch.ia(q,a)] += si;
            }
          }
        }
      };

      struct ReadoutFunctor
      {
        Ye3tView t;
        BondsView bv;
        ChunkView ch;
        OutputsView out;
        int c0 = 0;

        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t a ) const
        {
          const int i = ch.order[c0+a], mu = bv.type[i];
          const int s0 = t.src0[mu], ns = t.src0[mu+1] - s0, h0 = t.ch0[mu], nch = t.ch0[mu+1] - h0;
          const double fs = t.fs0[mu];
          // channel values
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c], tr = t.ch_tr[h0+c];
            const double cr = ch.Are[ch.ia(q,a)], ci = ch.Aim[ch.ia(q,a)];
            ch.Cre[ch.ic(c,a)] = tr == 0 ? cr : tr * cr;
            ch.Cim[ch.ic(c,a)] = tr == 0 ? ci : -tr * ci;
            ch.Gre[ch.ic(c,a)] = 0.0;
            ch.Gim[ch.ic(c,a)] = 0.0;
          }
          // Re rho: depth-first walk of the product tree, value of the path at each depth
          double pr[kYe3tMaxRank+1], pi[kYe3tMaxRank+1];
          pr[0] = 1.0; pi[0] = 0.0;
          double rho = 0.0;
          const int n0 = t.no0[mu], n1 = t.no0[mu+1];
          for(int n=n0;n<n1;n++)
          {
            const int d = t.node_depth[n], c = t.node_ch[n];
            const double cr = ch.Cre[ch.ic(c,a)], ci = ch.Cim[ch.ic(c,a)];
            pr[d] = pr[d-1] * cr - pi[d-1] * ci;
            pi[d] = pr[d-1] * ci + pi[d-1] * cr;
            rho += t.node_coef[n] * pr[d];
          }
          // reverse mode along the tree: adjoint of each node of the open path, closed when the walk leaves it:
          // G(channel) += adj * parent value, adj(parent) += adj * channel
          double ar[kYe3tMaxRank+1], ai[kYe3tMaxRank+1];
          int pc[kYe3tMaxRank+1];
          int top = 0;
          ar[0] = 0.0; ai[0] = 0.0;
          auto close = [&](int k)
          {
            const int c = pc[k];
            const double cr = ch.Cre[ch.ic(c,a)], ci = ch.Cim[ch.ic(c,a)];
            ch.Gre[ch.ic(c,a)] += ar[k] * pr[k-1] - ai[k] * pi[k-1];
            ch.Gim[ch.ic(c,a)] += ar[k] * pi[k-1] + ai[k] * pr[k-1];
            ar[k-1] += ar[k] * cr - ai[k] * ci;
            ai[k-1] += ar[k] * ci + ai[k] * cr;
          };
          for(int n=n0;n<n1;n++)
          {
            const int d = t.node_depth[n], c = t.node_ch[n];
            while( top >= d ) { close( top ); -- top; }
            const double cr = ch.Cre[ch.ic(c,a)], ci = ch.Cim[ch.ic(c,a)];
            pr[d] = pr[d-1] * cr - pi[d-1] * ci;
            pi[d] = pr[d-1] * ci + pi[d-1] * cr;
            ar[d] = fs * t.node_coef[n]; ai[d] = 0.0;
            pc[d] = c;
            top = d;
          }
          while( top >= 1 ) { close( top ); -- top; }
          // channel gradients folded into dE/d(source)
          for(int q=0;q<ns;q++) { ch.Wre[ch.ia(q,a)] = 0.0; ch.Wim[ch.ia(q,a)] = 0.0; }
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c], tr = t.ch_tr[h0+c];
            const double Gr = ch.Gre[ch.ic(c,a)], Gi = ch.Gim[ch.ic(c,a)];
            if( tr == 0 ) { ch.Wre[ch.ia(q,a)] += Gr; ch.Wim[ch.ia(q,a)] -= Gi; }
            else { ch.Wre[ch.ia(q,a)] += tr * Gr; ch.Wim[ch.ia(q,a)] += tr * Gi; }
          }
          out.eatom[i] += t.e0[mu] + fs * rho;
          if( rho >= t.rho_lim[mu] ) ONIKA_CU_ATOMIC_ADD( ch.flags[0] , 1 );
        }
      };

      // host readout of a group of kYe3tLanes atoms of one element, vectorised over the atoms (the scalar readout above for
      // a group mixing elements)
      static constexpr int kYe3tLanes = 8;
      struct ReadoutBatchFunctor
      {
        ReadoutFunctor one;
        int nc = 0;

        inline void operator () ( size_t g ) const
        {
          constexpr int L = kYe3tLanes;
          const Ye3tView& t = one.t;
          const ChunkView& ch = one.ch;
          const size_t a0 = g * L;
          const int nb = int( nc - a0 < size_t(L) ? nc - a0 : L );
          const int c0 = one.c0;
          const int mu = one.bv.type[ ch.order[c0+a0] ];
          bool same = nb == L;
          for(int b=1;b<nb && same;b++) same = one.bv.type[ ch.order[c0+a0+b] ] == mu;
          if( !same ) { for(int b=0;b<nb;b++) one( a0 + b ); return; }

          const int s0 = t.src0[mu], ns = t.src0[mu+1] - s0, h0 = t.ch0[mu], nch = t.ch0[mu+1] - h0;
          const double fs = t.fs0[mu];
          double* __restrict__ Cre = ch.Cre + ch.ic(0,a0);
          double* __restrict__ Cim = ch.Cim + ch.ic(0,a0);
          double* __restrict__ Gre = ch.Gre + ch.ic(0,a0);
          double* __restrict__ Gim = ch.Gim + ch.ic(0,a0);
          const double* __restrict__ Are = ch.Are + ch.ia(0,a0);
          const double* __restrict__ Aim = ch.Aim + ch.ia(0,a0);
          double* __restrict__ Wre = ch.Wre + ch.ia(0,a0);
          double* __restrict__ Wim = ch.Wim + ch.ia(0,a0);
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c], tr = t.ch_tr[h0+c];
            const double sr = tr == 0 ? 1.0 : tr, si = tr == 0 ? 1.0 : -tr;
#           pragma omp simd
            for(int b=0;b<L;b++) { Cre[c*L+b] = sr * Are[q*L+b]; Cim[c*L+b] = si * Aim[q*L+b]; Gre[c*L+b] = 0.0; Gim[c*L+b] = 0.0; }
          }
          double pr[kYe3tMaxRank+1][L], pi[kYe3tMaxRank+1][L], ar[kYe3tMaxRank+1][L], ai[kYe3tMaxRank+1][L], rho[L];
          for(int b=0;b<L;b++) { pr[0][b] = 1.0; pi[0][b] = 0.0; ar[0][b] = 0.0; ai[0][b] = 0.0; rho[b] = 0.0; }
          const int n0 = t.no0[mu], n1 = t.no0[mu+1];
          for(int n=n0;n<n1;n++)
          {
            const int d = t.node_depth[n], c = t.node_ch[n];
            const double cf = t.node_coef[n];
#           pragma omp simd
            for(int b=0;b<L;b++)
            {
              const double cr = Cre[c*L+b], ci = Cim[c*L+b];
              pr[d][b] = pr[d-1][b] * cr - pi[d-1][b] * ci;
              pi[d][b] = pr[d-1][b] * ci + pi[d-1][b] * cr;
              rho[b] += cf * pr[d][b];
            }
          }
          int pc[kYe3tMaxRank+1];
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
            const double cf = fs * t.node_coef[n];
#           pragma omp simd
            for(int b=0;b<L;b++)
            {
              const double cr = Cre[c*L+b], ci = Cim[c*L+b];
              pr[d][b] = pr[d-1][b] * cr - pi[d-1][b] * ci;
              pi[d][b] = pr[d-1][b] * ci + pi[d-1][b] * cr;
              ar[d][b] = cf; ai[d][b] = 0.0;
            }
            pc[d] = c;
            top = d;
          }
          while( top >= 1 ) { close( top ); -- top; }
          for(int q=0;q<ns;q++) for(int b=0;b<L;b++) { Wre[q*L+b] = 0.0; Wim[q*L+b] = 0.0; }
          for(int c=0;c<nch;c++)
          {
            const int q = t.ch_src[h0+c], tr = t.ch_tr[h0+c];
            const double sr = tr == 0 ? 1.0 : tr, si = tr == 0 ? -1.0 : tr;
#           pragma omp simd
            for(int b=0;b<L;b++) { Wre[q*L+b] += sr * Gre[c*L+b]; Wim[q*L+b] += si * Gim[c*L+b]; }
          }
          for(int b=0;b<L;b++)
          {
            const int i = ch.order[c0+a0+b];
            one.out.eatom[i] += t.e0[mu] + fs * rho[b];
            if( rho[b] >= t.rho_lim[mu] ) ONIKA_CU_ATOMIC_ADD( ch.flags[0] , 1 );
          }
        }
      };

      struct ForcesFunctor
      {
        Ye3tView t;
        BondsView bv;
        ChunkView ch;
        OutputsView out;
        int c0 = 0;

        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t a ) const
        {
          const int B = ONIKA_CU_BLOCK_SIZE, tid = ONIKA_CU_THREAD_IDX;
          const int i = ch.order[c0+a], mu = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[mu];
          double* rec = ch.scratch + block_slot(a) * ch.nscratch + size_t(tid) * forces_record( t );
          double* g = rec;
          double* dg = g + t.nbase_max;
          double* R = dg + t.nbase_max;
          double* dR = R + t.ncontr_max;
          double* yr = dR + t.ncontr_max;
          double* yi = yr + t.nlh;
          double* gr = yi + t.nlh;
          double* gi = gr + 3 * t.nlh;
          double fi[3] = { 0.0, 0.0, 0.0 }, vi[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          for(int jj=tid;jj<jnum;jj+=B)
          {
            const int j = bv.j[b0+jj], nu = bv.type[j];
            if( nu < 0 ) continue;
            const Ye3tBond& b = t.bonds[ mu * t.nel + nu ];
            const double* d = bv.rij + 3 * ( b0 + jj );
            const double dx = d[0], dy = d[1], dz = d[2], r2 = dx*dx + dy*dy + dz*dz;
            if( r2 >= b.rc * b.rc || r2 == 0.0 ) continue;
            const double r = sqrt( r2 ), inv_r = 1.0 / r, ux = dx * inv_r, uy = dy * inv_r, uz = dz * inv_r;
            radials<true>( t, b, r, g, dg, R, dR );
            ylm<true>( b.lb, t.ylm_a, t.ylm_b, ux, uy, uz, inv_r, yr, yi, gr, gi );
            double rad = 0.0, ang[3] = { 0.0, 0.0, 0.0 };
            const int qa = t.srcb[ mu * ( t.nel + 1 ) + nu ], qb = t.srcb[ mu * ( t.nel + 1 ) + nu + 1 ];
            for(int q=qa;q<qb;q++)
            {
              const int sq = s0 + q, rr = t.src_r[sq];
              const double wr = ch.Wre[ch.ia(q,a)];
              if( t.src_kind[sq] == 0 ) { rad += wr * dg[rr]; continue; }
              const double wi = ch.Wim[ch.ia(q,a)];
              const int lm = t.src_lm[sq];
              rad += dR[rr] * ( wr * yr[lm] + wi * yi[lm] );
              const double Rv = R[rr];
              for(int k=0;k<3;k++) ang[k] += Rv * ( wr * gr[3*lm+k] + wi * gi[3*lm+k] );
            }
            const double fv[3] = { rad * ux + ang[0], rad * uy + ang[1], rad * uz + ang[2] };    // dE_i / d r_ij
            for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(j) + k ] , -fv[k] );
            for(int k=0;k<3;k++) fi[k] += fv[k];
            // del (x) f with del = r_i - r_j, f = fv (xx yy zz xy xz yz)
            vi[0] -= dx * fv[0]; vi[1] -= dy * fv[1]; vi[2] -= dz * fv[2];
            vi[3] -= dx * fv[1]; vi[4] -= dx * fv[2]; vi[5] -= dy * fv[2];
          }
          for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(i) + k ] , fi[k] );
          if( out.vatom != nullptr ) for(int k=0;k<6;k++) ONIKA_CU_ATOMIC_ADD( out.vatom[ 6 * size_t(i) + k ] , vi[k] );
        }
      };

      // host basis and forces of one atom: its neighbours of each element in groups of kYe3tLanes, Y_lm and radial
      // functions of a group vectorised over the lanes. Thread scratch: neighbour list (jmax ints), then lane arrays.
      ONIKA_HOST_DEVICE_FUNC inline size_t host_lane_scratch(const Ye3tView& t, size_t jmax)
      {
        const size_t n = ( ( ( jmax + 1 ) / 2 + 2 ) & ~size_t(1) ) + size_t(kYe3tLanes) * ( 2 * size_t(t.nbase_max) + 2 * size_t(t.ncontr_max) + 8 * size_t(t.nlh) + 8 )
                       + 2 * size_t(t.nbase_max) + 2 * size_t(t.ncontr_max);
        return ( n + 1 ) & ~size_t(1);   // even: every thread slot 16-byte aligned
      }

      template<bool DER>
      struct HostLanes
      {
        static constexpr int L = kYe3tLanes;
        int* jl;
        double *g, *dg, *R, *dR, *yr, *yi, *gr, *gi, *ux, *uy, *uz, *ir, *dx, *dy, *dz, *tg, *tdg, *tR, *tdR;
        HostLanes(const Ye3tView& t, double* scr, size_t jmax)
        {
          jl = reinterpret_cast<int*>( scr );
          double* p = scr + ( ( ( jmax + 1 ) / 2 + 2 ) & ~size_t(1) );   // lane arrays 16-byte aligned (scr is)
          auto take = [&](size_t n) { double* x = p; p += n; return x; };
          g = take( L * t.nbase_max ); dg = take( L * t.nbase_max ); R = take( L * t.ncontr_max ); dR = take( L * t.ncontr_max );
          yr = take( L * t.nlh ); yi = take( L * t.nlh ); gr = take( 3 * L * t.nlh ); gi = take( 3 * L * t.nlh );
          ux = take( L ); uy = take( L ); uz = take( L ); ir = take( L ); dx = take( L ); dy = take( L ); dz = take( L ); take( L );
          tg = take( t.nbase_max ); tdg = take( t.nbase_max ); tR = take( t.ncontr_max ); tdR = take( t.ncontr_max );
        }
        // neighbours of element nu within the bond cutoff -> jl, count
        static int list( const Ye3tView& t, const BondsView& bv, int64_t b0, int jnum, int nu, double rc, int* jl )
        {
          int n = 0;
          for(int jj=0;jj<jnum;jj++)
          {
            if( bv.type[ bv.j[b0+jj] ] != nu ) continue;
            const double* d = bv.rij + 3 * ( b0 + jj );
            const double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            if( r2 < rc * rc && r2 > 0.0 ) jl[n++] = jj;
          }
          return n;
        }
        // geometry and radial functions of lanes 0 .. nb-1 (neighbours jl[k0 ..]), Y_lm of all lanes
        void fill( const Ye3tView& t, const Ye3tBond& b, const BondsView& bv, int64_t b0, int k0, int nb )
        {
          for(int x=0;x<L;x++)
          {
            if( x >= nb )
            {
              ux[x] = 0.0; uy[x] = 0.0; uz[x] = 1.0; ir[x] = 1.0; dx[x] = dy[x] = dz[x] = 0.0;
              for(int k=0;k<b.nbase;k++) { g[k*L+x] = 0.0; dg[k*L+x] = 0.0; }
              continue;
            }
            const double* d = bv.rij + 3 * ( b0 + jl[k0+x] );
            const double r = sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
            ir[x] = 1.0 / r; dx[x] = d[0]; dy[x] = d[1]; dz[x] = d[2];
            ux[x] = d[0] * ir[x]; uy[x] = d[1] * ir[x]; uz[x] = d[2] * ir[x];
            // base functions only (spline, or direct below the table / analytic mode); the contracted R_nl = sum_k c g_k
            // follow for all lanes at once (the cubic Hermite spline is linear in its node data: same values as the
            // spline of R_nl, whose nodes are the contractions of those of g_k)
            SplinePos sp;
            if( !t.analytic && sp.locate( b, r ) )
            {
              const double* bs = t.pool + b.bspl;
              double v, dv;
              for(int k=0;k<b.nbase;k++) { sp.eval( bs, b.nbase, k, v, dv ); g[k*L+x] = v; if( DER ) dg[k*L+x] = dv; }
            }
            else
            {
              cheb_exp_cos( b, r, tg, tdg );
              for(int k=0;k<b.nbase;k++) { g[k*L+x] = tg[k]; if( DER ) dg[k*L+x] = tdg[k]; }
            }
          }
          // R = C g, dR = C dg over the lanes, native 2-double vectors (lane arrays 16-byte aligned)
          {
            typedef double v2 __attribute__(( vector_size( 16 ) ));
            constexpr int W = L / 2;
            static_assert( L % 2 == 0 , "even lane count" );
            const double* c = t.pool + b.crad;
            const v2* G = reinterpret_cast<const v2*>( __builtin_assume_aligned( g , 16 ) );
            const v2* DG = reinterpret_cast<const v2*>( __builtin_assume_aligned( dg , 16 ) );
            v2* Rv = reinterpret_cast<v2*>( __builtin_assume_aligned( R , 16 ) );
            v2* dRv = reinterpret_cast<v2*>( __builtin_assume_aligned( dR , 16 ) );
            for(int q=0;q<b.ncontr;q++)
            {
              const double* cq = c + size_t(q) * b.nbase;
              v2 acc[W], dacc[W];
              for(int w=0;w<W;w++) { acc[w] = cq[0] * G[w]; if( DER ) dacc[w] = cq[0] * DG[w]; }
              for(int k=1;k<b.nbase;k++)
              {
                const double ck = cq[k];
                for(int w=0;w<W;w++) { acc[w] += ck * G[k*W+w]; if( DER ) dacc[w] += ck * DG[k*W+w]; }
              }
              for(int w=0;w<W;w++) { Rv[q*W+w] = acc[w]; if( DER ) dRv[q*W+w] = dacc[w]; }
            }
          }
          mlip_utils::ylm_lanes<L,DER>( b.lb, t.ylm_a, t.ylm_b, ux, uy, uz, ir, yr, yi, gr, gi );
        }
      };

      struct BasisHostFunctor
      {
        Ye3tView t;
        BondsView bv;
        ChunkView ch;
        int c0 = 0;
        size_t jmax = 0;

        inline void operator () ( size_t a ) const
        {
          constexpr int L = kYe3tLanes;
          const int i = ch.order[c0+a], mu = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[mu], ns = t.src0[mu+1] - s0;
          HostLanes<false> h( t, ch.scratch + block_slot(a) * ch.nscratch, jmax );
          for(int q=0;q<ns;q++) { ch.Are[ch.ia(q,a)] = 0.0; ch.Aim[ch.ia(q,a)] = 0.0; }
          for(int nu=0;nu<t.nel;nu++)
          {
            const Ye3tBond& b = t.bonds[ mu * t.nel + nu ];
            const int cnt = HostLanes<false>::list( t, bv, b0, jnum, nu, b.rc, h.jl );
            const int qa = t.srcb[ mu * ( t.nel + 1 ) + nu ], qb = t.srcb[ mu * ( t.nel + 1 ) + nu + 1 ];
            for(int k0=0;k0<cnt;k0+=L)
            {
              h.fill( t, b, bv, b0, k0, cnt - k0 < L ? cnt - k0 : L );
              for(int q=qa;q<qb;q++)
              {
                const int sq = s0 + q, rr = t.src_r[sq];
                double sr = 0.0, si = 0.0;
                if( t.src_kind[sq] == 0 )
                {
#                 pragma omp simd reduction(+:sr)
                  for(int x=0;x<L;x++) sr += h.g[rr*L+x];
                }
                else
                {
                  const int lm = t.src_lm[sq];
#                 pragma omp simd reduction(+:sr,si)
                  for(int x=0;x<L;x++) { sr += h.R[rr*L+x] * h.yr[lm*L+x]; si += h.R[rr*L+x] * h.yi[lm*L+x]; }
                }
                ch.Are[ch.ia(q,a)] += sr;
                ch.Aim[ch.ia(q,a)] += si;
              }
            }
          }
        }
      };

      struct ForcesHostFunctor
      {
        Ye3tView t;
        BondsView bv;
        ChunkView ch;
        OutputsView out;
        int c0 = 0;
        size_t jmax = 0;

        inline void operator () ( size_t a ) const
        {
          constexpr int L = kYe3tLanes;
          const int i = ch.order[c0+a], mu = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int s0 = t.src0[mu];
          HostLanes<true> h( t, ch.scratch + block_slot(a) * ch.nscratch, jmax );
          double fi[3] = { 0.0, 0.0, 0.0 }, vi[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          for(int nu=0;nu<t.nel;nu++)
          {
            const Ye3tBond& b = t.bonds[ mu * t.nel + nu ];
            const int cnt = HostLanes<true>::list( t, bv, b0, jnum, nu, b.rc, h.jl );
            const int qa = t.srcb[ mu * ( t.nel + 1 ) + nu ], qb = t.srcb[ mu * ( t.nel + 1 ) + nu + 1 ];
            for(int k0=0;k0<cnt;k0+=L)
            {
              const int nb = cnt - k0 < L ? cnt - k0 : L;
              h.fill( t, b, bv, b0, k0, nb );
              double rad[L], ax[L], ay[L], az[L];
              for(int x=0;x<L;x++) { rad[x] = 0.0; ax[x] = 0.0; ay[x] = 0.0; az[x] = 0.0; }
              for(int q=qa;q<qb;q++)
              {
                const int sq = s0 + q, rr = t.src_r[sq];
                const double wr = ch.Wre[ch.ia(q,a)];
                if( t.src_kind[sq] == 0 )
                {
#                 pragma omp simd
                  for(int x=0;x<L;x++) rad[x] += wr * h.dg[rr*L+x];
                  continue;
                }
                const double wi = ch.Wim[ch.ia(q,a)];
                const int lm = t.src_lm[sq];
                const double* __restrict__ gxr = h.gr + 3*lm*L; const double* __restrict__ gxi = h.gi + 3*lm*L;
#               pragma omp simd
                for(int x=0;x<L;x++)
                {
                  rad[x] += h.dR[rr*L+x] * ( wr * h.yr[lm*L+x] + wi * h.yi[lm*L+x] );
                  const double Rv = h.R[rr*L+x];
                  ax[x] += Rv * ( wr * gxr[x] + wi * gxi[x] );
                  ay[x] += Rv * ( wr * gxr[L+x] + wi * gxi[L+x] );
                  az[x] += Rv * ( wr * gxr[2*L+x] + wi * gxi[2*L+x] );
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

      struct ZeroIntFunctor
      {
        int* p = nullptr;
        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t k ) const { p[k] = 0; }
      };

    }
  }
}

namespace onika
{
  namespace parallel
  {
    template<> struct BlockParallelForFunctorTraits< exaStamp::mlip_ye3t::ye3t_onika::BasisFunctor > { static inline constexpr bool CudaCompatible = true; };
    template<> struct BlockParallelForFunctorTraits< exaStamp::mlip_ye3t::ye3t_onika::ForcesFunctor > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::mlip_ye3t::ye3t_onika::ReadoutFunctor > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::mlip_ye3t::ye3t_onika::ZeroIntFunctor > { static inline constexpr bool CudaCompatible = true; };
  }
}
