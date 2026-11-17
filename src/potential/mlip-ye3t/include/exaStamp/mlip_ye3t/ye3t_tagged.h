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

// Tagged-Cauchy correction of a YE3T composite model (ye3t-lammps model_family tagged_cauchy, slice v2 with a .yace
// backbone), flat tables and onika kernels, CPU and GPU, fp64. For a central atom i of element mu, over its edges
// j (r_ij = r_j - r_i, neighbour element nu):
//   components phi_c(r_ij) = R_n(r) Y_c(r_ij / r), n = radial channel of c (PACE ChebExpCos base g_n, analytic), Y_c a
//     real tesseral harmonic (a fixed linear form of the orthonormal complex Y_lm, m >= 0 stored), zero unless nu is
//     the neighbour element of c;
//   densities A_k = sum_j phi_c(k), moments M_m = sum_j prod_{c in key m} phi_c (same edge; empty key: M = z);
//   z = number of edges with r < max(tagged cutoff, backbone cutoff of (mu, nu));
//   E_i = offset[mu] + sum_terms beta[mu][f] coef prod A prod M F(z, p), F(z, p) = 0 if z < tag_count, else
//     prod_{q < tag_count - p} (z - p - q).
// The terms are a polynomial of the variables (A, M), kept as a prefix tree of their sorted factor lists per element,
// leaves weighted by F(z, p). One thread per central atom (GPU: one per atom of a chunk; host: OpenMP over atoms),
// per-slot scratch: pass 1 sums A, M and z, readout E and dE/dA, dE/dM (reverse mode along the tree), pass 2 the
// forces dE/dr_ij = sum_c dE/dphi_c dphi_c/dr_ij, dE/dphi_c from the densities and the moment tree in reverse.

#include <exaStamp/mlip_ye3t/ye3t_math.h>
#include <exaStamp/mlip_utils/mlip_flat_bonds.h>

#include <onika/cuda/cuda.h>
#include <onika/memory/allocator.h>
#include <onika/parallel/parallel_for.h>

#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace exaStamp
{
  namespace mlip_ye3t
  {
    static constexpr int kTaggedMaxDepth = 16;    // largest number of factors of a moment key or of a readout term

    struct Ye3tTaggedTables
    {
      std::vector<std::string> elements;
      double cutoff = 0.0;                            // tagged cutoff (edges counted within max(cutoff, backbone cutoff))
      int nel = 0;
      int ncomp = 0, nden = 0, nmom = 0, nvar = 0;    // nvar = nden + nmom
      int lmax = 0, nlh = 1;
      int tag_count = 0, pmax = 0;
      double rc_rad = 0.0;                            // radial cutoff (components vanish beyond)
      Ye3tBond rad;                                   // ChebExpCos parameters: rc, dcut, lambda, nbase = radial channels

      onika::memory::CudaMMVector<double> zcut2;      // [nel][nel] squared cutoff of the edge count z
      onika::memory::CudaMMVector<double> ylm_a, ylm_b;
      onika::memory::CudaMMVector<double> offset;     // [nel]

      // components: radial channel, neighbour element, l, offset of the (l + 1) pairs (re, im) of coefficients of
      // Re / Im Y_lm (m = 0 .. l, Y_00 = 1 normalisation) in cf
      onika::memory::CudaMMVector<int> comp_rad, comp_nu, comp_l, comp_cf;
      onika::memory::CudaMMVector<double> cf;
      onika::memory::CudaMMVector<int> den_comp;      // [nden] component of each density
      // host kernel: distinct angular forms (l, coefficient row of cf): ang_l, ang_cf; component -> form comp_ang
      int nang = 0;
      onika::memory::CudaMMVector<int> ang_l, ang_cf, comp_ang;
      // host kernel, moment tree: weight node -> first moment ending there handled per atom; first_child[n]: n is the
      // last child of its parent in pre-order (the first one met by the reverse sweep), has_child[n]
      onika::memory::CudaMMVector<int> mnode_first, mnode_haschild;
      // host kernel: slot of each internal node (forward values and reverse adjoints kept only for them), -1 for a
      // leaf; slot of the parent, -1 for a root
      int nslot = 0;
      onika::memory::CudaMMVector<int> mnode_slot, mnode_pslot;

      // moments with a non-empty key: prefix tree of the sorted component lists, pre-order; moments ending at node n:
      // mterm[mt0[n] .. mt0[n+1]). Empty keys (M = z): mzero
      onika::memory::CudaMMVector<int> mnode_comp, mnode_depth, mt0, mterm, mzero;
      onika::memory::CudaMMVector<int> mnode_parent;  // parent node, -1 for a root (host kernel: reverse sweep)

      // readout, per element: prefix tree of the sorted variable lists of the terms (variable v: density v < nden, else
      // moment v - nden), pre-order nodes ro0[mu] .. ro0[mu+1]; weights of the terms ending at node n:
      // rl[rl0[n] .. rl0[n+1]) = (p, coefficient beta * coef)
      onika::memory::CudaMMVector<int> ro0, rnode_var, rnode_depth, rl0, rl_p;
      onika::memory::CudaMMVector<double> rl_coef;
      onika::memory::CudaMMVector<int> r0terms;       // [nel+1] constant terms (no factor): rl index range r0t[mu] .. r0t[mu+1]
      onika::memory::CudaMMVector<int> r0t_p;
      onika::memory::CudaMMVector<double> r0t_coef;

      size_t bytes() const
      {
        return ( cf.size() + rl_coef.size() + r0t_coef.size() + zcut2.size() ) * 8
             + ( 4 * comp_rad.size() + den_comp.size() + 3 * mnode_comp.size() + mterm.size() + 3 * rnode_var.size() + rl_p.size() ) * 4;
      }
    };

    struct Ye3tTaggedView
    {
      int nel = 0, ncomp = 0, nden = 0, nmom = 0, nvar = 0, lmax = 0, nlh = 1, tag_count = 0, pmax = 0, nmzero = 0;
      double rc_rad = 0.0;
      Ye3tBond rad;
      const double *zcut2 = nullptr, *ylm_a = nullptr, *ylm_b = nullptr, *offset = nullptr, *cf = nullptr;
      const int *comp_rad = nullptr, *comp_nu = nullptr, *comp_l = nullptr, *comp_cf = nullptr, *den_comp = nullptr;
      int nmnode = 0, nang = 0;
      const int *ang_l = nullptr, *ang_cf = nullptr, *comp_ang = nullptr, *mnode_first = nullptr, *mnode_haschild = nullptr;
      int nslot = 0;
      const int *mnode_slot = nullptr, *mnode_pslot = nullptr;
      const int *mnode_comp = nullptr, *mnode_depth = nullptr, *mt0 = nullptr, *mterm = nullptr, *mzero = nullptr, *mnode_parent = nullptr;
      const int *ro0 = nullptr, *rnode_var = nullptr, *rnode_depth = nullptr, *rl0 = nullptr, *rl_p = nullptr;
      const double* rl_coef = nullptr;
      const int *r0terms = nullptr, *r0t_p = nullptr;
      const double* r0t_coef = nullptr;
    };

    inline Ye3tTaggedView make_tagged_view( const Ye3tTaggedTables& T )
    {
      Ye3tTaggedView v;
      v.nel = T.nel; v.ncomp = T.ncomp; v.nden = T.nden; v.nmom = T.nmom; v.nvar = T.nvar; v.lmax = T.lmax; v.nlh = T.nlh;
      v.tag_count = T.tag_count; v.pmax = T.pmax; v.nmzero = T.mzero.size(); v.rc_rad = T.rc_rad; v.rad = T.rad;
      v.zcut2 = T.zcut2.data(); v.ylm_a = T.ylm_a.data(); v.ylm_b = T.ylm_b.data(); v.offset = T.offset.data(); v.cf = T.cf.data();
      v.comp_rad = T.comp_rad.data(); v.comp_nu = T.comp_nu.data(); v.comp_l = T.comp_l.data(); v.comp_cf = T.comp_cf.data();
      v.den_comp = T.den_comp.data();
      v.nmnode = T.mnode_comp.size();
      v.nang = T.nang; v.ang_l = T.ang_l.data(); v.ang_cf = T.ang_cf.data(); v.comp_ang = T.comp_ang.data();
      v.mnode_first = T.mnode_first.data(); v.mnode_haschild = T.mnode_haschild.data();
      v.nslot = T.nslot; v.mnode_slot = T.mnode_slot.data(); v.mnode_pslot = T.mnode_pslot.data();
      v.mnode_comp = T.mnode_comp.data(); v.mnode_depth = T.mnode_depth.data(); v.mt0 = T.mt0.data(); v.mterm = T.mterm.data(); v.mzero = T.mzero.data(); v.mnode_parent = T.mnode_parent.data();
      v.ro0 = T.ro0.data(); v.rnode_var = T.rnode_var.data(); v.rnode_depth = T.rnode_depth.data(); v.rl0 = T.rl0.data(); v.rl_p = T.rl_p.data();
      v.rl_coef = T.rl_coef.data();
      v.r0terms = T.r0terms.data(); v.r0t_p = T.r0t_p.data(); v.r0t_coef = T.r0t_coef.data();
      return v;
    }

    struct Ye3tTaggedBuffers
    {
      onika::memory::CudaMMVector<double> scratch;
    };

    namespace ye3t_onika
    {
      using mlip_utils::BondsView;
      using mlip_utils::OutputsView;

      // scratch of one slot: A, M, dE/dA, dE/dM, dE/dphi from the densities, phi, dphi, dE/dphi, g, dg, Y, dY
      ONIKA_HOST_DEVICE_FUNC inline size_t tagged_record( const Ye3tTaggedView& t )
      {
        return 2 * size_t(t.nvar) + 6 * size_t(t.ncomp) + 2 * size_t(t.rad.nbase) + 8 * size_t(t.nlh);
      }

      // F(z, p): free count of the remaining tag positions
      ONIKA_HOST_DEVICE_FUNC inline double tagged_free_count( int z, int p, int tag_count )
      {
        if( z < tag_count ) return 0.0;
        double v = 1.0;
        for(int q=0;q<tag_count-p;q++) v *= double( z - p - q );
        return v;
      }

      struct TaggedFunctor
      {
        Ye3tTaggedView t;
        BondsView bv;
        OutputsView out;
        double* scratch = nullptr;
        size_t nscratch = 0;
        int c0 = 0;

        // components of one edge of neighbour element nu (and their gradients)
        template<bool DER>
        ONIKA_HOST_DEVICE_FUNC inline void components( int nu, const double* d, double r, double* phi, double* dphi,
                                                       double* g, double* dg, double* yr, double* yi, double* gr, double* gi ) const
        {
          const double inv_r = 1.0 / r, ux = d[0] * inv_r, uy = d[1] * inv_r, uz = d[2] * inv_r;
          cheb_exp_cos( t.rad, r, g, dg );
          mlip_utils::ylm<DER>( t.lmax, t.ylm_a, t.ylm_b, ux, uy, uz, inv_r, yr, yi, gr, gi );
          for(int c=0;c<t.ncomp;c++)
          {
            if( t.comp_nu[c] != nu ) { phi[c] = 0.0; if( DER ) { dphi[3*c] = 0.0; dphi[3*c+1] = 0.0; dphi[3*c+2] = 0.0; } continue; }
            const int l = t.comp_l[c], lm0 = mlip_utils::lmh( l, 0 );
            const double* w = t.cf + t.comp_cf[c];
            double y = 0.0, gy[3] = { 0.0, 0.0, 0.0 };
            for(int m=0;m<=l;m++)
            {
              const double wr = w[2*m], wi = w[2*m+1];
              y += wr * yr[lm0+m] + wi * yi[lm0+m];
              if( DER ) for(int k=0;k<3;k++) gy[k] += wr * gr[3*(lm0+m)+k] + wi * gi[3*(lm0+m)+k];
            }
            const int n = t.comp_rad[c];
            phi[c] = g[n] * y;
            if( DER )
            {
              const double u[3] = { ux, uy, uz };
              for(int k=0;k<3;k++) dphi[3*c+k] = dg[n] * u[k] * y + g[n] * gy[k];
            }
          }
        }

        ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t a ) const
        {
          const int i = bv.ilist[c0+a], ti = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          double* s = scratch + ONIKA_CU_VALUE_IF_CUDA( a , size_t( ONIKA_CU_BLOCK_IDX ) ) * nscratch;
          double* X = s;                       // [nvar] A then M
          double* GX = X + t.nvar;             // [nvar]
          double* GA = GX + t.nvar;            // [ncomp] dE/dphi from the densities
          double* phi = GA + t.ncomp;          // [ncomp]
          double* gphi = phi + t.ncomp;        // [ncomp]
          double* dphi = gphi + t.ncomp;       // [ncomp][3]
          double* g = dphi + 3 * t.ncomp;      // [nbase]
          double* dg = g + t.rad.nbase;
          double* yr = dg + t.rad.nbase;       // [nlh]
          double* yi = yr + t.nlh;
          double* gr = yi + t.nlh;             // [nlh][3]
          double* gi = gr + 3 * t.nlh;
          double* A = X;
          double* M = X + t.nden;
          const double rc2 = t.rc_rad * t.rc_rad;
          const double* zc2 = t.zcut2 + ti * t.nel;

          // pass 1: densities, moments, edge count
          for(int v=0;v<t.nvar;v++) { X[v] = 0.0; GX[v] = 0.0; }
          int z = 0;
          double pv[kTaggedMaxDepth+1];
          pv[0] = 1.0;
          for(int jj=0;jj<jnum;jj++)
          {
            const int j = bv.j[b0+jj], nu = bv.type[j];
            if( nu < 0 ) continue;
            const double* d = bv.rij + 3 * ( b0 + jj );
            const double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            if( r2 < zc2[nu] ) ++ z;
            if( r2 >= rc2 || r2 == 0.0 ) continue;
            components<false>( nu, d, sqrt( r2 ), phi, nullptr, g, dg, yr, yi, nullptr, nullptr );
            for(int k=0;k<t.nden;k++) A[k] += phi[ t.den_comp[k] ];
            for(int n=0;n<t.nmnode;n++)
            {
              const int dp = t.mnode_depth[n];
              pv[dp] = pv[dp-1] * phi[ t.mnode_comp[n] ];
              for(int q=t.mt0[n];q<t.mt0[n+1];q++) M[ t.mterm[q] ] += pv[dp];
            }
          }
          for(int q=0;q<t.nmzero;q++) M[ t.mzero[q] ] = double( z );

          // readout: energy, dE/dX
          double Fp[kTaggedMaxDepth+1];
          for(int p=0;p<=t.pmax;p++) Fp[p] = tagged_free_count( z, p, t.tag_count );
          double e = t.offset[ti];
          for(int q=t.r0terms[ti];q<t.r0terms[ti+1];q++) e += t.r0t_coef[q] * Fp[ t.r0t_p[q] ];
          {
            double ar[kTaggedMaxDepth+1];
            int pvar[kTaggedMaxDepth+1];
            int top = 0;
            ar[0] = 0.0;
            for(int n=t.ro0[ti];n<t.ro0[ti+1];n++)
            {
              const int dp = t.rnode_depth[n], v = t.rnode_var[n];
              while( top >= dp )
              {
                const int u = pvar[top];
                GX[u] += ar[top] * pv[top-1];
                ar[top-1] += ar[top] * X[u];
                -- top;
              }
              pv[dp] = pv[dp-1] * X[v];
              double w = 0.0;
              for(int q=t.rl0[n];q<t.rl0[n+1];q++) w += t.rl_coef[q] * Fp[ t.rl_p[q] ];
              e += w * pv[dp];
              ar[dp] = w;
              pvar[dp] = v;
              top = dp;
            }
            while( top >= 1 )
            {
              const int u = pvar[top];
              GX[u] += ar[top] * pv[top-1];
              ar[top-1] += ar[top] * X[u];
              -- top;
            }
          }
          out.eatom[i] += e;

          // pass 2: forces
          for(int c=0;c<t.ncomp;c++) GA[c] = 0.0;
          for(int k=0;k<t.nden;k++) GA[ t.den_comp[k] ] += GX[k];
          const double* GM = GX + t.nden;
          double fi[3] = { 0.0, 0.0, 0.0 }, vi[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          for(int jj=0;jj<jnum;jj++)
          {
            const int j = bv.j[b0+jj], nu = bv.type[j];
            if( nu < 0 ) continue;
            const double* d = bv.rij + 3 * ( b0 + jj );
            const double dx = d[0], dy = d[1], dz = d[2], r2 = dx*dx + dy*dy + dz*dz;
            if( r2 >= rc2 || r2 == 0.0 ) continue;
            components<true>( nu, d, sqrt( r2 ), phi, dphi, g, dg, yr, yi, gr, gi );
            for(int c=0;c<t.ncomp;c++) gphi[c] = GA[c];
            {
              double ar[kTaggedMaxDepth+1];
              int pc[kTaggedMaxDepth+1];
              int top = 0;
              ar[0] = 0.0;
              for(int n=0;n<t.nmnode;n++)
              {
                const int dp = t.mnode_depth[n], c = t.mnode_comp[n];
                while( top >= dp )
                {
                  gphi[ pc[top] ] += ar[top] * pv[top-1];
                  ar[top-1] += ar[top] * phi[ pc[top] ];
                  -- top;
                }
                pv[dp] = pv[dp-1] * phi[c];
                double w = 0.0;
                for(int q=t.mt0[n];q<t.mt0[n+1];q++) w += GM[ t.mterm[q] ];
                ar[dp] = w;
                pc[dp] = c;
                top = dp;
              }
              while( top >= 1 )
              {
                gphi[ pc[top] ] += ar[top] * pv[top-1];
                ar[top-1] += ar[top] * phi[ pc[top] ];
                -- top;
              }
            }
            double fv[3] = { 0.0, 0.0, 0.0 };    // dE_i / d r_ij
            for(int c=0;c<t.ncomp;c++)
            {
              const double w = gphi[c];
              if( w == 0.0 ) continue;
              fv[0] += w * dphi[3*c]; fv[1] += w * dphi[3*c+1]; fv[2] += w * dphi[3*c+2];
            }
            for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(j) + k ] , -fv[k] );
            for(int k=0;k<3;k++) fi[k] += fv[k];
            vi[0] -= dx * fv[0]; vi[1] -= dy * fv[1]; vi[2] -= dz * fv[2];
            vi[3] -= dx * fv[1]; vi[4] -= dx * fv[2]; vi[5] -= dy * fv[2];
          }
          for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(i) + k ] , fi[k] );
          if( out.vatom != nullptr ) for(int k=0;k<6;k++) ONIKA_CU_ATOMIC_ADD( out.vatom[ 6 * size_t(i) + k ] , vi[k] );
        }
      };

      static constexpr int kTaggedLanes = 8;
      // lane operations on kTaggedLanes doubles, separate arguments are never aliased
      namespace lanes
      {
        // lanes as native 2-double vectors (GCC / clang vector extension, 16-byte aligned: lane arrays start on 64 bytes)
        constexpr int L = kTaggedLanes, W = kTaggedLanes / 2;
        typedef double v2 __attribute__(( vector_size( 16 ) ));
        static_assert( kTaggedLanes % 2 == 0 , "even lane count" );
        inline v2* V( double* p ) { return reinterpret_cast<v2*>( __builtin_assume_aligned( p , 16 ) ); }
        inline const v2* V( const double* p ) { return reinterpret_cast<const v2*>( __builtin_assume_aligned( p , 16 ) ); }
        inline void zero( double* o ) { v2* a = V(o); for(int k=0;k<W;k++) a[k] = v2{ 0.0, 0.0 }; }
        inline void set( double* o, double c ) { v2* a = V(o); for(int k=0;k<W;k++) a[k] = v2{ c, c }; }
        inline void add( double* o, const double* x ) { v2* a = V(o); const v2* b = V(x); for(int k=0;k<W;k++) a[k] += b[k]; }
        inline void add_c( double* o, double c ) { v2* a = V(o); for(int k=0;k<W;k++) a[k] += c; }
        // o = x (root) and acc += x
        inline void copy_acc( double* o, double* acc, const double* x ) { v2* a = V(o); v2* s = V(acc); const v2* b = V(x); for(int k=0;k<W;k++) { a[k] = b[k]; s[k] += b[k]; } }
        // o = x y and acc += o
        inline void mul_acc( double* o, double* acc, const double* x, const double* y ) { v2* a = V(o); v2* s = V(acc); const v2* b = V(x); const v2* c = V(y); for(int k=0;k<W;k++) { const v2 v = b[k] * c[k]; a[k] = v; s[k] += v; } }
        // o += x y
        inline void fma( double* o, const double* x, const double* y ) { v2* a = V(o); const v2* b = V(x); const v2* c = V(y); for(int k=0;k<W;k++) a[k] += b[k] * c[k]; }
        // o = x y + c
        inline void mul_add_c( double* o, const double* x, const double* y, double c ) { v2* a = V(o); const v2* b = V(x); const v2* d = V(y); for(int k=0;k<W;k++) a[k] = b[k] * d[k] + c; }
        // o += c y
        inline void axpy( double* o, double c, const double* y ) { v2* a = V(o); const v2* b = V(y); for(int k=0;k<W;k++) a[k] += c * b[k]; }
        // o = c y + d
        inline void scale_add( double* o, double c, const double* y, double d ) { v2* a = V(o); const v2* b = V(y); for(int k=0;k<W;k++) a[k] = c * b[k] + d; }
        inline double sum( const double* x ) { v2 r = { 0.0, 0.0 }; const v2* b = V(x); for(int k=0;k<W;k++) r += b[k]; return r[0] + r[1]; }
      }
      inline double* tprof_stage() { static double t[4] = { 0, 0, 0, 0 }; return t; }   // YE3T_PROF (not thread safe: profile at 1 thread)

      // host scratch of one OpenMP thread, in doubles: lane arrays (component and node accumulators, reverse sweep,
      // per edge block the components, their gradients and the forward node values), then per-atom scalars
      inline size_t host_tagged_scratch( const Ye3tTaggedView& t, size_t jmax )
      {
        const size_t L = kTaggedLanes, nb = ( jmax + L - 1 ) / L + 1, C = t.ncomp, N = t.nmnode;
        const size_t NS = std::max( 1 , t.nslot );
        return ( C + N + nb * ( 6 * C + 2 * NS ) ) * L
             + 2 * size_t(t.nvar) + C + N + 4 * size_t(t.nang)
             + 2 * size_t(t.rad.nbase) + 8 * size_t(t.nlh) + ( jmax + 1 ) + L
             + ( 4 + 2 * size_t(t.rad.nbase) + 8 * size_t(t.nlh) + 4 * size_t(t.nang) + 1 ) * L;   // lane geometry; rounded to L by the driver
      }

      // host kernel: one atom per call, its edges by blocks of kTaggedLanes in SIMD lanes (padding lanes: zero components)
      struct TaggedHostFunctor
      {
        Ye3tTaggedView t;
        BondsView bv;
        OutputsView out;
        double* scratch = nullptr;
        size_t nscratch = 0;
        size_t jmax = 0;
        bool prof = false;

        inline void operator () ( size_t a ) const
        {
          constexpr int L = kTaggedLanes;
          const int i = bv.ilist[a], ti = bv.type[i];
          const int64_t b0 = bv.offset[i];
          const int jnum = int( bv.offset[i+1] - b0 );
          const int C = t.ncomp, N = t.nmnode;
          const size_t nbmax = ( jmax + L - 1 ) / L + 1;
          double* s = scratch + size_t( omp_get_thread_num() ) * nscratch;
          double* Cacc = s;                                        // [C][L]
          double* Nacc = Cacc + size_t(C) * L;                     // [N][L]
          const int NS = t.nslot > 0 ? t.nslot : 1;
          double* S = Nacc + size_t(N) * L;                        // [nb][NS][L] reverse adjoints of the internal nodes
          double* GP = S + nbmax * NS * L;                         // [nb][C][L]
          double* PH = GP + nbmax * C * L;                         // [nb][C][L]
          double* DPH = PH + nbmax * C * L;                        // [nb][3 C][L]
          double* PV = DPH + nbmax * 3 * C * L;                    // [nb][NS][L] forward values of the internal nodes
          double* X = PV + nbmax * NS * L;
          double* GX = X + t.nvar;
          double* GA = GX + t.nvar;
          double* GMn = GA + C;                                    // [N] dE/d(node moments)
          double* ya = GMn + N;                               // [nang] real harmonics of one edge, [nang][3] gradients
          double* g = ya + 4 * t.nang;
          double* dg = g + t.rad.nbase;
          double* yr = dg + t.rad.nbase; double* yi = yr + t.nlh; double* gr = yi + t.nlh; double* gi = gr + 3 * t.nlh;
          double* lu = gi + 3 * t.nlh;                              // [3][L] unit vectors
          double* lir = lu + 3 * L;                                 // [L] 1 / r
          double* lg = lir + L;                                     // [nbase][L]
          double* ldg = lg + t.rad.nbase * L;
          double* lyr = ldg + t.rad.nbase * L;                      // [nlh][L]
          double* lyi = lyr + t.nlh * L;
          double* lgr = lyi + t.nlh * L;                            // [nlh][3][L]
          double* lgi = lgr + 3 * t.nlh * L;
          double* lya = lgi + 3 * t.nlh * L;                        // [nang][4][L] real harmonics, gradients
          int* lnu = reinterpret_cast<int*>( lya + 4 * t.nang * L ); // [L]
          int* el = lnu + L;
          const double rc2 = t.rc_rad * t.rc_rad;
          const double* zc2 = t.zcut2 + ti * t.nel;
          // YE3T_PROF: cumulative stage times (components, pass 1, readout, pass 2), host only
          std::chrono::steady_clock::time_point T0;
          if( prof ) T0 = std::chrono::steady_clock::now();
          auto lap = [&](int k) { if( !prof ) return; const auto T1 = std::chrono::steady_clock::now(); tprof_stage()[k] += std::chrono::duration<double>( T1 - T0 ).count(); T0 = T1; };

          // edges: count z, keep those with nonzero components
          int z = 0, ne = 0;
          for(int jj=0;jj<jnum;jj++)
          {
            const int j = bv.j[b0+jj], nu = bv.type[j];
            if( nu < 0 ) continue;
            const double* d = bv.rij + 3 * ( b0 + jj );
            const double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            if( r2 < zc2[nu] ) ++ z;
            if( r2 < rc2 && r2 > 0.0 ) el[ne++] = jj;
          }
          const int nb = ( ne + L - 1 ) / L;

          // components and gradients of every edge (lane layout), zero in the padding: per block, geometry and radial
          // base per lane, harmonics and their real forms vectorised over the lanes, components masked by element
          for(int b=0;b<nb;b++)
          {
            double* ph = PH + size_t(b) * C * L;
            double* dph = DPH + size_t(b) * 3 * C * L;
            for(int l=0;l<L;l++)
            {
              const int e = b * L + l;
              if( e >= ne )
              {
                lu[l] = 0.0; lu[L+l] = 0.0; lu[2*L+l] = 1.0; lir[l] = 1.0; lnu[l] = -1;
                for(int k=0;k<t.rad.nbase;k++) { lg[k*L+l] = 0.0; ldg[k*L+l] = 0.0; }
                continue;
              }
              const int jj = el[e];
              lnu[l] = bv.type[ bv.j[b0+jj] ];
              const double* d = bv.rij + 3 * ( b0 + jj );
              const double r = sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] ), inv_r = 1.0 / r;
              lir[l] = inv_r; lu[l] = d[0] * inv_r; lu[L+l] = d[1] * inv_r; lu[2*L+l] = d[2] * inv_r;
              cheb_exp_cos( t.rad, r, g, dg );
              for(int k=0;k<t.rad.nbase;k++) { lg[k*L+l] = g[k]; ldg[k*L+l] = dg[k]; }
            }
            mlip_utils::ylm_lanes<L,true>( t.lmax, t.ylm_a, t.ylm_b, lu, lu + L, lu + 2 * L, lir, lyr, lyi, lgr, lgi );
            for(int q=0;q<t.nang;q++)
            {
              const int lc = t.ang_l[q], lm0 = mlip_utils::lmh( lc, 0 );
              const double* w = t.cf + t.ang_cf[q];
              double* __restrict__ y = lya + size_t(4*q) * L;
              for(int k=0;k<4*L;k++) y[k] = 0.0;
              for(int m=0;m<=lc;m++)
              {
                const double wr = w[2*m], wi = w[2*m+1];
                const int h = lm0 + m;
                const double* __restrict__ yrh = lyr + h * L; const double* __restrict__ yih = lyi + h * L;
                const double* __restrict__ grh = lgr + 3 * h * L; const double* __restrict__ gih = lgi + 3 * h * L;
#               pragma omp simd
                for(int l=0;l<L;l++)
                {
                  y[l] += wr * yrh[l] + wi * yih[l];
                  y[L+l] += wr * grh[l] + wi * gih[l];
                  y[2*L+l] += wr * grh[L+l] + wi * gih[L+l];
                  y[3*L+l] += wr * grh[2*L+l] + wi * gih[2*L+l];
                }
              }
            }
            for(int c=0;c<C;c++)
            {
              const int nu = t.comp_nu[c];
              const double* __restrict__ y = lya + size_t(4 * t.comp_ang[c]) * L;
              const double* __restrict__ gn = lg + t.comp_rad[c] * L;
              const double* __restrict__ dgn = ldg + t.comp_rad[c] * L;
              double* __restrict__ o = ph + c * L;
              double* __restrict__ o0 = dph + (3*c) * L; double* __restrict__ o1 = o0 + L; double* __restrict__ o2 = o1 + L;
#             pragma omp simd
              for(int l=0;l<L;l++)
              {
                const double mk = lnu[l] == nu ? 1.0 : 0.0;
                const double gv = mk * gn[l], dgv = mk * dgn[l] * y[l];
                o[l] = gv * y[l];
                o0[l] = dgv * lu[l] + gv * y[L+l];
                o1[l] = dgv * lu[L+l] + gv * y[2*L+l];
                o2[l] = dgv * lu[2*L+l] + gv * y[3*L+l];
              }
            }
          }
          lap(0);

          // pass 1: densities and moments, lane accumulators reduced once
          for(int c=0;c<C;c++) lanes::zero( Cacc + c * L );
          for(int n=0;n<N;n++) lanes::zero( Nacc + n * L );
          const size_t phs = size_t(C) * L, pvs_ = size_t(NS) * L;   // block strides
          for(int b=0;b<nb;b++) for(int c=0;c<C;c++) lanes::add( Cacc + c * L , PH + b * phs + c * L );
          for(int n=0;n<N;n++)   // node outer, blocks inner: the node data once for all blocks, independent blocks
          {
            const int ps = t.mnode_pslot[n], c = t.mnode_comp[n], sl = t.mnode_slot[n];
            double* na = Nacc + n * L;
            const double* ph = PH + c * L;
            if( sl >= 0 )
            {
              double* pv = PV + sl * L;
              if( ps < 0 ) for(int b=0;b<nb;b++) lanes::copy_acc( pv + b * pvs_ , na , ph + b * phs );
              else for(int b=0;b<nb;b++) lanes::mul_acc( pv + b * pvs_ , na , PV + ps * L + b * pvs_ , ph + b * phs );
            }
            else if( ps < 0 ) for(int b=0;b<nb;b++) lanes::add( na , ph + b * phs );
            else for(int b=0;b<nb;b++) lanes::fma( na , PV + ps * L + b * pvs_ , ph + b * phs );
          }
          auto hsum = [](const double* v) { return lanes::sum( v ); };
          for(int v=0;v<t.nvar;v++) { X[v] = 0.0; GX[v] = 0.0; }
          double* A = X;
          double* M = X + t.nden;
          for(int k=0;k<t.nden;k++) A[k] = hsum( Cacc + t.den_comp[k] * L );
          for(int n=0;n<N;n++)
          {
            if( t.mt0[n] == t.mt0[n+1] ) continue;
            const double sum = hsum( Nacc + n * L );
            for(int q=t.mt0[n];q<t.mt0[n+1];q++) M[ t.mterm[q] ] = sum;
          }
          for(int q=0;q<t.nmzero;q++) M[ t.mzero[q] ] = double( z );
          lap(1);

          // readout (scalar, depth-first walk of the term tree, the open path on a stack)
          double Fp[kTaggedMaxDepth+1], pvs[kTaggedMaxDepth+1];
          for(int p=0;p<=t.pmax;p++) Fp[p] = tagged_free_count( z, p, t.tag_count );
          double e = t.offset[ti];
          for(int q=t.r0terms[ti];q<t.r0terms[ti+1];q++) e += t.r0t_coef[q] * Fp[ t.r0t_p[q] ];
          {
            double ar[kTaggedMaxDepth+1];
            int pvar[kTaggedMaxDepth+1];
            int top = 0;
            pvs[0] = 1.0; ar[0] = 0.0;
            for(int n=t.ro0[ti];n<t.ro0[ti+1];n++)
            {
              const int dp = t.rnode_depth[n], v = t.rnode_var[n];
              while( top >= dp ) { const int u = pvar[top]; GX[u] += ar[top] * pvs[top-1]; ar[top-1] += ar[top] * X[u]; -- top; }
              pvs[dp] = pvs[dp-1] * X[v];
              double w = 0.0;
              for(int q=t.rl0[n];q<t.rl0[n+1];q++) w += t.rl_coef[q] * Fp[ t.rl_p[q] ];
              e += w * pvs[dp];
              ar[dp] = w; pvar[dp] = v; top = dp;
            }
            while( top >= 1 ) { const int u = pvar[top]; GX[u] += ar[top] * pvs[top-1]; ar[top-1] += ar[top] * X[u]; -- top; }
          }
          out.eatom[i] += e;
          lap(2);
          if( nb == 0 ) return;

          // pass 2: dE/dphi per edge by a reverse sweep of the moment tree (children before parents), forces
          for(int c=0;c<C;c++) GA[c] = 0.0;
          for(int k=0;k<t.nden;k++) GA[ t.den_comp[k] ] += GX[k];
          const double* GM = GX + t.nden;
          for(int n=0;n<N;n++) { double w = 0.0; for(int q=t.mt0[n];q<t.mt0[n+1];q++) w += GM[ t.mterm[q] ]; GMn[n] = w; }
          double fi[3] = { 0.0, 0.0, 0.0 }, vi[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          const size_t gps = size_t(C) * L;
          for(int b=0;b<nb;b++) for(int c=0;c<C;c++) lanes::set( GP + b * gps + c * L , GA[c] );
          for(int n=N-1;n>=0;n--)   // reverse sweep, node outer, blocks inner
          {
            const int ps = t.mnode_pslot[n], c = t.mnode_comp[n], sl = t.mnode_slot[n];
            double* gp = GP + c * L;
            const double* ph = PH + c * L;
            if( sl < 0 )   // leaf: adjoint = its moment weight, the same in every lane
            {
              const double g = GMn[n];
              if( ps < 0 ) { for(int b=0;b<nb;b++) lanes::add_c( gp + b * gps , g ); continue; }
              const double* pv = PV + ps * L;
              double* sp = S + ps * L;
              for(int b=0;b<nb;b++) lanes::axpy( gp + b * gps , g , pv + b * pvs_ );
              if( t.mnode_first[n] ) { const double gm = GMn[ t.mnode_parent[n] ]; for(int b=0;b<nb;b++) lanes::scale_add( sp + b * pvs_ , g , ph + b * phs , gm ); }
              else for(int b=0;b<nb;b++) lanes::axpy( sp + b * pvs_ , g , ph + b * phs );
              continue;
            }
            const double* sn = S + sl * L;
            if( ps < 0 ) { for(int b=0;b<nb;b++) lanes::add( gp + b * gps , sn + b * pvs_ ); continue; }
            const double* pv = PV + ps * L;
            double* sp = S + ps * L;
            for(int b=0;b<nb;b++) lanes::fma( gp + b * gps , sn + b * pvs_ , pv + b * pvs_ );
            if( t.mnode_first[n] ) { const double gm = GMn[ t.mnode_parent[n] ]; for(int b=0;b<nb;b++) lanes::mul_add_c( sp + b * pvs_ , sn + b * pvs_ , ph + b * phs , gm ); }
            else for(int b=0;b<nb;b++) lanes::fma( sp + b * pvs_ , sn + b * pvs_ , ph + b * phs );
          }
          for(int b=0;b<nb;b++)
          {
            const double* dph = DPH + size_t(b) * 3 * C * L;
            const double* GPb = GP + b * gps;
            double f0[L], f1[L], f2[L];
            lanes::zero( f0 ); lanes::zero( f1 ); lanes::zero( f2 );
            for(int c=0;c<C;c++) { lanes::fma( f0 , GPb + c * L , dph + (3*c) * L ); lanes::fma( f1 , GPb + c * L , dph + (3*c+1) * L ); lanes::fma( f2 , GPb + c * L , dph + (3*c+2) * L ); }
            for(int l=0;l<L && b*L+l<ne;l++)
            {
              const int jj = el[b*L+l], j = bv.j[b0+jj];
              const double* d = bv.rij + 3 * ( b0 + jj );
              const double x0 = f0[l], x1 = f1[l], x2 = f2[l];
              ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(j) ] , -x0 );
              ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(j) + 1 ] , -x1 );
              ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(j) + 2 ] , -x2 );
              fi[0] += x0; fi[1] += x1; fi[2] += x2;
              vi[0] -= d[0] * x0; vi[1] -= d[1] * x1; vi[2] -= d[2] * x2;
              vi[3] -= d[0] * x1; vi[4] -= d[0] * x2; vi[5] -= d[1] * x2;
            }
          }
          for(int k=0;k<3;k++) ONIKA_CU_ATOMIC_ADD( out.f[ 3 * size_t(i) + k ] , fi[k] );
          if( out.vatom != nullptr ) for(int k=0;k<6;k++) ONIKA_CU_ATOMIC_ADD( out.vatom[ 6 * size_t(i) + k ] , vi[k] );
          lap(3);
        }
      };
    }

    // tagged correction added into out (after or before the backbone, same bonds). gpu: the kernel runs on the GPU;
    // mem_gb: budget of the per-atom scratch (GPU); pec: launcher (mlip_onika_launch.h or the emulator)
    template<class LaunchT>
    inline void ye3t_tagged_compute( const Ye3tTaggedTables& T, Ye3tTaggedBuffers& buf, const mlip_utils::BondsView& bv, mlip_utils::OutputsView& out,
                                     bool gpu, double mem_gb, const LaunchT& pec )
    {
      using namespace ye3t_onika;
      const int inum = bv.inum;
      if( inum <= 0 ) return;
      const Ye3tTaggedView tv = make_tagged_view( T );
      static const bool prof = std::getenv( "YE3T_PROF" ) != nullptr;
      static double tprof = 0.0;
      const auto t0 = std::chrono::steady_clock::now();
      if( gpu )
      {
        const size_t rec = tagged_record( tv );
        const double budget = ( mem_gb > 0 ? mem_gb : 1.0 ) * 1e9 / 8.0;
        const int chunk = std::max( 1 , std::min( inum , int( std::min( budget / double(rec) , 2.0e9 ) ) ) );
        buf.scratch.resize( rec * size_t(chunk) );
        for(int c0=0;c0<inum;c0+=chunk)
        {
          const int nc = std::min( chunk , inum - c0 );
          pec.flat( nc , TaggedFunctor{ tv, bv, out, buf.scratch.data(), rec, c0 } , "ye3t_tagged" );
        }
      }
      else
      {
        size_t jmax = 0;
        for(int q=0;q<inum;q++) { const int i = bv.ilist[q]; jmax = std::max( jmax , size_t( bv.offset[i+1] - bv.offset[i] ) ); }
        const size_t rec = ( host_tagged_scratch( tv, jmax ) + kTaggedLanes - 1 ) / kTaggedLanes * kTaggedLanes;
        // the allocator only guarantees the SIMD alignment of the build (32 bytes with AVX2): over-allocate one cache
        // line and align the base to 64 bytes, so that every per-thread slot (rec is a multiple of 8 doubles) starts on a line
        buf.scratch.resize( rec * size_t( omp_get_max_threads() ) + 8 );
        double* scr = buf.scratch.data();
        scr += ( ( 64 - reinterpret_cast<uintptr_t>( scr ) % 64 ) % 64 ) / sizeof(double);
#       pragma omp parallel for schedule(dynamic,16)
        for(int q=0;q<inum;q++) TaggedHostFunctor{ tv, bv, out, scr, rec, jmax, prof }( q );
      }
      if( prof )
      {
        tprof += std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
        std::fprintf( stderr, "YE3T_PROF tagged %.3f s (comp %.3f pass1 %.3f readout %.3f pass2 %.3f)\n", tprof, tprof_stage()[0], tprof_stage()[1], tprof_stage()[2], tprof_stage()[3] );
      }
      if( gpu ) ONIKA_CU_DEVICE_SYNCHRONIZE();
    }

  }
}

namespace onika
{
  namespace parallel
  {
    template<> struct ParallelForFunctorTraits< exaStamp::mlip_ye3t::ye3t_onika::TaggedFunctor > { static inline constexpr bool CudaCompatible = true; };
  }
}
