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

// Strain Functional Descriptors (SFD), per-particle math, host and device.
//
// Reference: E. M. Kober, J. P. Tavenner, C. M. Adams, N. Mathew, "Strain functionals:
// A complete and symmetry-adapted set of descriptors to characterize atomistic
// configurations", arXiv:2402.04191.
//
// For a central particle i and its neighbors j (i itself included, weight 1), with
// r_ij = r_i - r_j in units of sigma and w_j = exp(-|r_ij|^2/2):
//   W = sum_j w_j, c = sum_j w_j r_ij / W, u_j = r_ij - c,
//   M_n = (1/W) sum_j w_j u_j^(x)n (n = 2, 3, 4, fully symmetric Cartesian tensors).
// Each M_n is split into its SO(3) irreducible parts (solid harmonics, Condon-Shortley
// phase, normalized so that the l=0 part of a unit scalar is 1). The descriptors are
// rotational invariants built from these parts with Clebsch-Gordan couplings, plus
// orientation factors (projections on the lab Z/X/Y axes).
//
// Implementation notes:
// - Raw weighted monomial sums m_pqr = sum_j w_j x^p y^q z^r (p+q+r <= 4) are accumulated in
//   a single neighbor pass, then shifted to the weighted centroid with the binomial formula.
// - Spherical tensors of rank L store the m >= 0 components only; a_{-m} = (-1)^m conj(a_m)
//   (the tensors come from real Cartesian tensors).
// - Clebsch-Gordan coefficients are evaluated with Racah's closed formula; all ranks are <= 4.

#include <onika/cuda/cuda.h>
#include <cmath>

namespace exaStamp
{
  namespace sfd
  {
    static constexpr int MAX_ORDER = 4;
    static constexpr int N_MONOMIALS = 35;   // monomials x^p y^q z^r with p+q+r <= 4
    static constexpr int N_INVARIANTS = 36;  // W, 1 + 3 + 7 + 12 invariants, 4 x 3 orientation factors
    static constexpr int N_SPHERICAL = 35;   // real components of the spherical tensors (see sfd_spherical_names)
    static constexpr double TINY = 1.0e-15;

    // ----------------------------------------------------------------------------------
    // small complex type usable on device
    // ----------------------------------------------------------------------------------
    struct Cplx
    {
      double re = 0.0;
      double im = 0.0;
    };
    ONIKA_HOST_DEVICE_FUNC inline Cplx operator + (Cplx a, Cplx b) { return { a.re+b.re , a.im+b.im }; }
    ONIKA_HOST_DEVICE_FUNC inline Cplx operator - (Cplx a, Cplx b) { return { a.re-b.re , a.im-b.im }; }
    ONIKA_HOST_DEVICE_FUNC inline Cplx operator - (Cplx a) { return { -a.re , -a.im }; }
    ONIKA_HOST_DEVICE_FUNC inline Cplx operator * (Cplx a, Cplx b) { return { a.re*b.re - a.im*b.im , a.re*b.im + a.im*b.re }; }
    ONIKA_HOST_DEVICE_FUNC inline Cplx operator * (double s, Cplx a) { return { s*a.re , s*a.im }; }
    ONIKA_HOST_DEVICE_FUNC inline Cplx conj(Cplx a) { return { a.re , -a.im }; }

    // spherical tensor of rank L, components m = 0..L
    template<int L> struct SphT
    {
      Cplx c[L+1] = {};
    };

    template<int L>
    ONIKA_HOST_DEVICE_FUNC inline Cplx component( const SphT<L>& a, int m )
    {
      if( m >= 0 ) return a.c[m];
      const Cplx z = conj( a.c[-m] );
      return ( (-m) & 1 ) ? -z : z;
    }

    // ----------------------------------------------------------------------------------
    // Clebsch-Gordan coefficients <j1 m1 j2 m2 | j m> (Racah formula)
    // ----------------------------------------------------------------------------------
    ONIKA_HOST_DEVICE_FUNC inline double factorial( int n )
    {
      double f = 1.0;
      for(int i=2;i<=n;i++) f *= i;
      return f;
    }

    ONIKA_HOST_DEVICE_FUNC inline double clebsch_gordan( int j1, int m1, int j2, int m2, int j, int m )
    {
      if( m1+m2 != m ) return 0.0;
      if( m1>j1 || -m1>j1 || m2>j2 || -m2>j2 || m>j || -m>j ) return 0.0;
      if( j > j1+j2 || j < j1-j2 || j < j2-j1 ) return 0.0;
      const double tri = (2*j+1) * factorial(j1+j2-j) * factorial(j1-j2+j) * factorial(-j1+j2+j) / factorial(j1+j2+j+1);
      const double norm = factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j+m) * factorial(j-m);
      int kmin = 0;
      if( j2-j-m1 > kmin ) kmin = j2-j-m1;
      if( j1-j+m2 > kmin ) kmin = j1-j+m2;
      int kmax = j1+j2-j;
      if( j1-m1 < kmax ) kmax = j1-m1;
      if( j2+m2 < kmax ) kmax = j2+m2;
      double s = 0.0;
      for(int k=kmin;k<=kmax;k++)
      {
        const double d = factorial(k) * factorial(j1+j2-j-k) * factorial(j1-m1-k) * factorial(j2+m2-k) * factorial(j-j2+m1+k) * factorial(j-j1-m2+k);
        s += ( (k&1) ? -1.0 : 1.0 ) / d;
      }
      return sqrt( tri * norm ) * s;
    }

    // invariant scalar product of two rank L tensors: sum_m <L m L -m|0 0> a_m b_-m,
    // i.e. ( a_0 b_0 + 2 sum_{m>0} Re(a_m conj(b_m)) ) / sqrt(2L+1)
    template<int L>
    ONIKA_HOST_DEVICE_FUNC inline double dot( const SphT<L>& a, const SphT<L>& b )
    {
      double s = a.c[0].re * b.c[0].re - a.c[0].im * b.c[0].im;
      for(int m=1;m<=L;m++) s += 2.0 * ( a.c[m].re * b.c[m].re + a.c[m].im * b.c[m].im );
      return s / sqrt( double(2*L+1) );
    }

    // rank L coupling of a (rank L1) and b (rank L2): sum_{m1+m2=m} <L1 m1 L2 m2|L m> a_m1 b_m2
    template<int L, int L1, int L2>
    ONIKA_HOST_DEVICE_FUNC inline SphT<L> couple( const SphT<L1>& a, const SphT<L2>& b )
    {
      SphT<L> r;
      for(int m=0;m<=L;m++)
      {
        Cplx s = {};
        for(int m1=-L1;m1<=L1;m1++)
        {
          const int m2 = m - m1;
          if( m2 < -L2 || m2 > L2 ) continue;
          s = s + clebsch_gordan(L1,m1,L2,m2,L,m) * ( component(a,m1) * component(b,m2) );
        }
        r.c[m] = s;
      }
      return r;
    }

    template<int L> ONIKA_HOST_DEVICE_FUNC inline double norm( const SphT<L>& a ) { return sqrt( dot(a,a) ); }
    ONIKA_HOST_DEVICE_FUNC inline double max_tiny( double x ) { return ( x > TINY ) ? x : TINY; }

    // orientation factors of a rank L tensor: projections on the lab unit tensors
    // Z (m=0 component), X and Y (real and imaginary m=L components), scaled by (2L+1)^(1/4)
    // so that a tensor aligned with one of them gives 1. For L=1 the m=1 unit tensors carry
    // the Condon-Shortley sign, so that X/Y match the x/y axes.
    template<int L>
    ONIKA_HOST_DEVICE_FUNC inline void orientation( const SphT<L>& a, double inv_norm, double* out )
    {
      const double sgn = ( L == 1 ) ? -1.0 : 1.0;
      const double s = sqrt( sqrt( double(2*L+1) ) ) * inv_norm;
      SphT<L> e;
      e.c[0] = Cplx{ 1.0 , 0.0 };
      out[0] = s * dot( a , e );
      e.c[0] = Cplx{};
      e.c[L] = Cplx{ sgn / sqrt(2.0) , 0.0 };
      out[1] = s * dot( a , e );
      e.c[L] = Cplx{ 0.0 , sgn / sqrt(2.0) };
      out[2] = s * dot( a , e );
    }

    // ----------------------------------------------------------------------------------
    // monomial indexing: degree n block starts at n(n+1)(n+2)/6, inside a block the
    // monomial x^p y^q z^r is at s(s+1)/2 + r with s = q + r
    // ----------------------------------------------------------------------------------
    ONIKA_HOST_DEVICE_FUNC inline constexpr int monomial_index( int p, int q, int r )
    {
      const int n = p+q+r;
      const int s = q+r;
      return n*(n+1)*(n+2)/6 + s*(s+1)/2 + r;
    }

    // adds w x^p y^q z^r to all monomials
    ONIKA_HOST_DEVICE_FUNC inline void accumulate_monomials( double* m, double w, double x, double y, double z )
    {
      double px[MAX_ORDER+1], py[MAX_ORDER+1], pz[MAX_ORDER+1];
      px[0] = py[0] = pz[0] = 1.0;
      for(int k=1;k<=MAX_ORDER;k++) { px[k] = px[k-1]*x; py[k] = py[k-1]*y; pz[k] = pz[k-1]*z; }
      for(int n=0;n<=MAX_ORDER;n++)
      for(int s=0;s<=n;s++)
      for(int r=0;r<=s;r++)
      {
        const int q = s-r, p = n-s;
        m[ monomial_index(p,q,r) ] += w * px[p] * py[q] * pz[r];
      }
    }

    ONIKA_HOST_DEVICE_FUNC inline double binomial( int n, int k )
    {
      return factorial(n) / ( factorial(k) * factorial(n-k) );
    }

    // ----------------------------------------------------------------------------------
    // spherical parts of symmetric Cartesian tensors
    // ----------------------------------------------------------------------------------
    // rank 1 (vector)
    ONIKA_HOST_DEVICE_FUNC inline SphT<1> spherical_rank1( double vx, double vy, double vz )
    {
      SphT<1> t;
      t.c[0] = Cplx{ sqrt(2.0) * vz , 0.0 };
      t.c[1] = Cplx{ -vx , -vy };
      return t;
    }

    // traceless part (rank 2) of a symmetric matrix
    ONIKA_HOST_DEVICE_FUNC inline SphT<2> spherical_rank2( double xx, double yy, double zz, double xy, double xz, double yz )
    {
      SphT<2> t;
      t.c[0] = Cplx{ ( 2.0*zz - xx - yy ) / sqrt(3.0) , 0.0 };
      t.c[1] = Cplx{ -sqrt(2.0) * xz , -sqrt(2.0) * yz };
      t.c[2] = Cplx{ ( xx - yy ) / sqrt(2.0) , sqrt(2.0) * xy };
      return t;
    }

    // all spherical parts of the centered moments
    struct SFDSpherical
    {
      double w = 0.0;   // total weight
      SphT<1> t1;       // centroid offset
      double s2 = 0.0;  // order 2 scalar part (isotropic Gaussian reference removed)
      SphT<2> t2;       // order 2, rank 2
      SphT<1> t31;      // order 3, rank 1 (trace vector)
      SphT<3> t3;       // order 3, rank 3
      double s4 = 0.0;  // order 4 scalar part (isotropic Gaussian reference removed)
      SphT<2> t42;      // order 4, rank 2 (order 2 contribution removed)
      SphT<4> t4;       // order 4, rank 4
    };

    // raw: monomial sums m_pqr of r_ij/sigma, self term included in m_000
    ONIKA_HOST_DEVICE_FUNC inline SFDSpherical spherical_decomposition( const double* raw )
    {
      SFDSpherical out;
      const double W = raw[0];
      const double iW = 1.0 / W;
      const double c[3] = { raw[monomial_index(1,0,0)]*iW , raw[monomial_index(0,1,0)]*iW , raw[monomial_index(0,0,1)]*iW };

      // centered normalized moments mu_pqr = (1/W) sum_j w_j (x-cx)^p (y-cy)^q (z-cz)^r, degrees 2..4
      double mu[N_MONOMIALS];
      for(int i=0;i<N_MONOMIALS;i++) mu[i] = 0.0;
      double nc[3][MAX_ORDER+1];
      for(int d=0;d<3;d++) { nc[d][0] = 1.0; for(int k=1;k<=MAX_ORDER;k++) nc[d][k] = nc[d][k-1] * (-c[d]); }
      for(int n=2;n<=MAX_ORDER;n++)
      for(int s=0;s<=n;s++)
      for(int r=0;r<=s;r++)
      {
        const int q = s-r, p = n-s;
        double v = 0.0;
        for(int a=0;a<=p;a++)
        for(int b=0;b<=q;b++)
        for(int e=0;e<=r;e++)
        {
          v += binomial(p,a) * binomial(q,b) * binomial(r,e) * nc[0][p-a] * nc[1][q-b] * nc[2][r-e] * raw[ monomial_index(a,b,e) ];
        }
        mu[ monomial_index(p,q,r) ] = v * iW;
      }
#     define SFD_MU(p,q,r) mu[ monomial_index(p,q,r) ]

      out.w = W;
      out.t1 = spherical_rank1( c[0] , c[1] , c[2] );

      // order 2
      const double Sxx = SFD_MU(2,0,0), Syy = SFD_MU(0,2,0), Szz = SFD_MU(0,0,2);
      const double Sxy = SFD_MU(1,1,0), Sxz = SFD_MU(1,0,1), Syz = SFD_MU(0,1,1);
      const double trS = Sxx + Syy + Szz;
      out.s2 = sqrt(2.0/3.0) * trS - sqrt(1.5);
      out.t2 = spherical_rank2( Sxx, Syy, Szz, Sxy, Sxz, Syz );

      // order 3: trace vector V_a = sum_b T_abb, and traceless rank 3 part
      const double Vx = SFD_MU(3,0,0) + SFD_MU(1,2,0) + SFD_MU(1,0,2);
      const double Vy = SFD_MU(2,1,0) + SFD_MU(0,3,0) + SFD_MU(0,1,2);
      const double Vz = SFD_MU(2,0,1) + SFD_MU(0,2,1) + SFD_MU(0,0,3);
      out.t31 = spherical_rank1( Vx, Vy, Vz );
      for(int m=0;m<=1;m++) out.t31.c[m] = sqrt(0.4) * out.t31.c[m];
      out.t3.c[0] = Cplx{ sqrt(2.0/15.0) * ( 2.0*SFD_MU(0,0,3) - 3.0*SFD_MU(2,0,1) - 3.0*SFD_MU(0,2,1) ) , 0.0 };
      out.t3.c[1] = Cplx{ ( SFD_MU(3,0,0) + SFD_MU(1,2,0) - 4.0*SFD_MU(1,0,2) ) / sqrt(10.0)
                        , ( SFD_MU(2,1,0) + SFD_MU(0,3,0) - 4.0*SFD_MU(0,1,2) ) / sqrt(10.0) };
      out.t3.c[2] = Cplx{ SFD_MU(2,0,1) - SFD_MU(0,2,1) , 2.0*SFD_MU(1,1,1) };
      out.t3.c[3] = Cplx{ ( 3.0*SFD_MU(1,2,0) - SFD_MU(3,0,0) ) / sqrt(6.0)
                        , ( SFD_MU(0,3,0) - 3.0*SFD_MU(2,1,0) ) / sqrt(6.0) };

      // order 4: trace matrix D_ab = sum_c Q_abcc, scalar part, rank 2 and rank 4 parts
      const double Dxx = SFD_MU(4,0,0) + SFD_MU(2,2,0) + SFD_MU(2,0,2);
      const double Dyy = SFD_MU(2,2,0) + SFD_MU(0,4,0) + SFD_MU(0,2,2);
      const double Dzz = SFD_MU(2,0,2) + SFD_MU(0,2,2) + SFD_MU(0,0,4);
      const double Dxy = SFD_MU(3,1,0) + SFD_MU(1,3,0) + SFD_MU(1,1,2);
      const double Dxz = SFD_MU(3,0,1) + SFD_MU(1,2,1) + SFD_MU(1,0,3);
      const double Dyz = SFD_MU(2,1,1) + SFD_MU(0,3,1) + SFD_MU(0,1,3);
      // for an isotropic unit Gaussian, trS = 3 and trD = 15: both scalar parts vanish
      out.s4 = sqrt(2.0/15.0) * ( Dxx + Dyy + Dzz ) - sqrt(10.0/3.0) * trS + sqrt(15.0/8.0);
      const SphT<2> d2 = spherical_rank2( Dxx, Dyy, Dzz, Dxy, Dxz, Dyz );
      for(int m=0;m<=2;m++) out.t42.c[m] = sqrt(2.0/7.0) * d2.c[m] - sqrt(3.5) * out.t2.c[m];

      out.t4.c[0] = Cplx{ ( 3.0*SFD_MU(4,0,0) + 3.0*SFD_MU(0,4,0) + 8.0*SFD_MU(0,0,4)
                          + 6.0*SFD_MU(2,2,0) - 24.0*SFD_MU(2,0,2) - 24.0*SFD_MU(0,2,2) ) / sqrt(420.0) , 0.0 };
      out.t4.c[1] = Cplx{ ( 12.0*SFD_MU(3,0,1) + 12.0*SFD_MU(1,2,1) - 16.0*SFD_MU(1,0,3) ) / sqrt(336.0)
                        , ( 12.0*SFD_MU(2,1,1) + 12.0*SFD_MU(0,3,1) - 16.0*SFD_MU(0,1,3) ) / sqrt(336.0) };
      out.t4.c[2] = Cplx{ ( 2.0*( SFD_MU(0,4,0) - SFD_MU(4,0,0) ) + 12.0*( SFD_MU(2,0,2) - SFD_MU(0,2,2) ) ) / sqrt(42.0)
                        , ( 24.0*SFD_MU(1,1,2) - 4.0*SFD_MU(3,1,0) - 4.0*SFD_MU(1,3,0) ) / sqrt(42.0) };
      out.t4.c[3] = Cplx{ ( 12.0*SFD_MU(1,2,1) - 4.0*SFD_MU(3,0,1) ) / sqrt(48.0)
                        , ( 4.0*SFD_MU(0,3,1) - 12.0*SFD_MU(2,1,1) ) / sqrt(48.0) };
      out.t4.c[4] = Cplx{ ( SFD_MU(4,0,0) - 6.0*SFD_MU(2,2,0) + SFD_MU(0,4,0) ) / sqrt(24.0)
                        , ( 4.0*SFD_MU(3,1,0) - 4.0*SFD_MU(1,3,0) ) / sqrt(24.0) };
#     undef SFD_MU
      return out;
    }

    // 35 real values: w, t1, s2, t2, t31, t3, s4, t42, t4 (m = 0 real part, then re/im for m > 0)
    template<int L>
    ONIKA_HOST_DEVICE_FUNC inline int write_spherical( const SphT<L>& t, double* out )
    {
      out[0] = t.c[0].re;
      for(int m=1;m<=L;m++) { out[2*m-1] = t.c[m].re; out[2*m] = t.c[m].im; }
      return 2*L+1;
    }

    ONIKA_HOST_DEVICE_FUNC inline void spherical_values( const SFDSpherical& s, double* out )
    {
      int k = 0;
      out[k++] = s.w;
      k += write_spherical( s.t1 , out+k );
      out[k++] = s.s2;
      k += write_spherical( s.t2 , out+k );
      k += write_spherical( s.t31 , out+k );
      k += write_spherical( s.t3 , out+k );
      out[k++] = s.s4;
      k += write_spherical( s.t42 , out+k );
      k += write_spherical( s.t4 , out+k );
    }

    // ----------------------------------------------------------------------------------
    // rotational invariants and orientation factors, 36 values (see sfd_invariant_names)
    // ----------------------------------------------------------------------------------
    ONIKA_HOST_DEVICE_FUNC inline void invariants( const SFDSpherical& s, double* out )
    {
      double* o = out;
      *o++ = s.w;

      // order 1
      const double n1 = norm( s.t1 );
      *o++ = n1;

      // order 2: norm, skewness, scalar part
      const double n2 = norm( s.t2 );
      const double i2 = 1.0 / max_tiny( n2 );
      *o++ = n2;
      *o++ = dot( s.t2 , couple<2>( s.t2 , s.t2 ) ) * i2 * i2;
      *o++ = s.s2;

      // order 3
      {
        const SphT<3>& a = s.t3;
        const SphT<1>& b = s.t31;
        const double n3 = norm( a );
        const double i3 = 1.0 / max_tiny( n3 );
        const SphT<2> a2 = couple<2>( a , a );
        const double n32 = norm( a2 );
        const SphT<2> a22 = couple<2>( a2 , a2 );
        const SphT<4> a4 = couple<4>( a , a );
        const double n34 = norm( a4 );
        const SphT<4> a44 = couple<4>( a4 , a4 );
        const double n31 = norm( b );
        const SphT<2> b2 = couple<2>( b , b );
        const double n12 = norm( b2 );
        const double n322 = norm( a22 );
        const double sc = sqrt( n3 * n31 );
        *o++ = n3;
        *o++ = n32 * i3;
        *o++ = dot( a2 , a22 ) * i3 / max_tiny( n32*n32 );
        *o++ = dot( a4 , a44 ) * i3 / max_tiny( n34*n34 );
        *o++ = n31;
        *o++ = dot( a2 , b2 ) * sc / max_tiny( n32*n12 );
        *o++ = dot( a22 , b2 ) * sc / max_tiny( n322*n12 );
      }

      // order 4
      {
        const SphT<4>& a = s.t4;
        const SphT<2>& b = s.t42;
        const double n4 = norm( a );
        const double i4 = 1.0 / max_tiny( n4 );
        const SphT<2> a2 = couple<2>( a , a );
        const double n42 = norm( a2 );
        const double i42 = 1.0 / max_tiny( n42 );
        const SphT<2> a22 = couple<2>( a2 , a2 );
        const SphT<4> a4 = couple<4>( a , a );
        const SphT<4> a44 = couple<4>( a , a4 );
        const SphT<4> a4d4 = couple<4>( a4 , a4 );
        const double n44 = max_tiny( dot( a4 , a4 ) );
        const SphT<4> a224 = couple<4>( a2 , a2 );
        const double nb = norm( b );
        const double ib = 1.0 / max_tiny( nb );
        const SphT<2> b22 = couple<2>( b , b );
        const double sc = sqrt( n4 * nb );
        *o++ = n4;
        *o++ = n42 * i4;
        *o++ = dot( a2 , a22 ) * i4 * i42 * i42;
        *o++ = dot( a4 , a4d4 ) * i4 / n44;
        *o++ = dot( a , a224 ) * i4 * i42 * sqrt( i42 );
        *o++ = dot( a44 , a4d4 ) * i4 * i4 / n44;
        *o++ = nb;
        *o++ = dot( b , b22 ) * ib * ib;
        *o++ = s.s4;
        *o++ = dot( a2 , b ) * sc / max_tiny( n42 * nb );
        *o++ = dot( a22 , b ) * sc / max_tiny( n42 * n42 * nb );
        *o++ = dot( a2 , b22 ) * sc / max_tiny( n42 * nb * nb );
      }

      // orientation factors
      orientation( s.t1 , 1.0 / max_tiny( n1 ) , o ); o += 3;
      orientation( s.t2 , i2 , o ); o += 3;
      orientation( s.t3 , 1.0 / max_tiny( norm(s.t3) ) , o ); o += 3;
      orientation( s.t4 , 1.0 / max_tiny( norm(s.t4) ) , o ); o += 3;
    }

    // ----------------------------------------------------------------------------------
    // column names (host only)
    // ----------------------------------------------------------------------------------
    inline const char* const* sfd_invariant_names()
    {
      static const char* const names[N_INVARIANTS] = {
        "w", "P1_I0",
        "P2_I0", "P2_I1", "P2_I2",
        "P3_I0", "P3_I1", "P3_I2", "P3_I3", "P3_I4", "P3_I5", "P3_I6",
        "P4_I0", "P4_I1", "P4_I2", "P4_I3", "P4_I4", "P4_I5", "P4_I6", "P4_I7", "P4_I8", "P4_I9", "P4_I10", "P4_I11",
        "O1_I0", "O1_I1", "O1_I2", "O2_I0", "O2_I1", "O2_I2", "O3_I0", "O3_I1", "O3_I2", "O4_I0", "O4_I1", "O4_I2" };
      return names;
    }

    inline const char* const* sfd_spherical_names()
    {
      static const char* const names[N_SPHERICAL] = {
        "P00_0", "P11_0", "P11_1r", "P11_1i",
        "P20_0", "P22_0", "P22_1r", "P22_1i", "P22_2r", "P22_2i",
        "P31a_0", "P31a_1r", "P31a_1i",
        "P33a_0", "P33a_1r", "P33a_1i", "P33a_2r", "P33a_2i", "P33a_3r", "P33a_3i",
        "P40_0", "P42a_0", "P42a_1r", "P42a_1i", "P42a_2r", "P42a_2i",
        "P44a_0", "P44a_1r", "P44a_1i", "P44a_2r", "P44a_2i", "P44a_3r", "P44a_3i", "P44a_4r", "P44a_4i" };
      return names;
    }

  }
}
