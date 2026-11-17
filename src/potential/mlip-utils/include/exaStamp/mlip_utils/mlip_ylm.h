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

// Complex spherical harmonics Y_lm (m >= 0), Condon-Shortley phase, normalised to Y_00 = 1 (orthonormal ones:
// multiply by 1 / sqrt(4 pi)), and their gradient with respect to the bond vector, for device / host kernels.
// Storage index lmh(l, m) = l (l + 1) / 2 + m. Y_{l,-m} = (-1)^m conj(Y_lm).

#include <onika/cuda/cuda.h>
#include <cmath>

namespace exaStamp
{
  namespace mlip_utils
  {
    ONIKA_HOST_DEVICE_FUNC inline int lmh(int l, int m) { return l * ( l + 1 ) / 2 + m; }

    // Y_lm (m >= 0, l <= lmax), Condon-Shortley phase, Y_00 = 1, at the unit vector u, and dY_lm / d(bond vector)
    // (bond length r). Y[lmh] (re, im) in yr / yi, gradients [lmh][3] in gr / gi.
    // Recurrences on the unit sphere: Y_mm = a_mm (ux + i uy) Y_{m-1,m-1}, Y_lm = a_lm uz Y_{l-1,m} - b_lm Y_{l-2,m}; the
    // partial derivatives along u of these polynomial extensions are projected on the tangent plane, divided by r.
    template<bool GRAD>
    ONIKA_HOST_DEVICE_FUNC inline void ylm( int lmax, const double* A, const double* B, double ux, double uy, double uz, double inv_r,
                                            double* yr, double* yi, double* gr, double* gi )
    {
      yr[0] = 1.0; yi[0] = 0.0;
      if( GRAD ) for(int k=0;k<3;k++) { gr[k] = 0.0; gi[k] = 0.0; }
      for(int m=0;m<=lmax;m++)
      {
        const int d = lmh(m,m);
        if( m > 0 )
        {
          const int p = lmh(m-1,m-1);
          const double a = A[d];
          const double pr = yr[p], pi = yi[p];
          yr[d] = a * ( ux * pr - uy * pi );
          yi[d] = a * ( ux * pi + uy * pr );
          if( GRAD )
          {
            for(int k=0;k<3;k++)
            {
              const double qr = gr[3*p+k], qi = gi[3*p+k];
              gr[3*d+k] = a * ( ux * qr - uy * qi );
              gi[3*d+k] = a * ( ux * qi + uy * qr );
            }
            gr[3*d+0] += a * pr;  gi[3*d+0] += a * pi;    // d/dux (ux + i uy) = 1
            gr[3*d+1] -= a * pi;  gi[3*d+1] += a * pr;    // d/duy (ux + i uy) = i
          }
        }
        for(int l=m+1;l<=lmax;l++)
        {
          const int q = lmh(l,m), q1 = lmh(l-1,m);
          const double a = A[q], bb = B[q];
          const bool two = l - 2 >= m;
          const int q2 = two ? lmh(l-2,m) : 0;
          yr[q] = a * uz * yr[q1] - ( two ? bb * yr[q2] : 0.0 );
          yi[q] = a * uz * yi[q1] - ( two ? bb * yi[q2] : 0.0 );
          if( GRAD )
          {
            for(int k=0;k<3;k++)
            {
              gr[3*q+k] = a * uz * gr[3*q1+k] - ( two ? bb * gr[3*q2+k] : 0.0 );
              gi[3*q+k] = a * uz * gi[3*q1+k] - ( two ? bb * gi[3*q2+k] : 0.0 );
            }
            gr[3*q+2] += a * yr[q1];
            gi[3*q+2] += a * yi[q1];
          }
        }
      }
      if( GRAD )
      {
        const int n = lmh(lmax,lmax) + 1;
        for(int q=0;q<n;q++)
        {
          const double sr = ux * gr[3*q] + uy * gr[3*q+1] + uz * gr[3*q+2];
          const double si = ux * gi[3*q] + uy * gi[3*q+1] + uz * gi[3*q+2];
          gr[3*q] = ( gr[3*q] - sr * ux ) * inv_r;  gr[3*q+1] = ( gr[3*q+1] - sr * uy ) * inv_r;  gr[3*q+2] = ( gr[3*q+2] - sr * uz ) * inv_r;
          gi[3*q] = ( gi[3*q] - si * ux ) * inv_r;  gi[3*q+1] = ( gi[3*q+1] - si * uy ) * inv_r;  gi[3*q+2] = ( gi[3*q+2] - si * uz ) * inv_r;
        }
      }
    }

    // host: Y_lm (and gradients) of L bond vectors at once, vectorised over the lanes (2-double vectors, L even, no
    // alignment requirement); arrays [lmh][L] (yr, yi) and [lmh][3][L] (gr, gi). Unused lanes must hold a valid unit
    // vector (their results are ignored).
    template<int L, bool GRAD>
    inline void ylm_lanes( int lmax, const double* A, const double* B, const double* ux_, const double* uy_, const double* uz_, const double* inv_r_,
                           double* __restrict__ yr_, double* __restrict__ yi_, double* __restrict__ gr_, double* __restrict__ gi_ )
    {
      static_assert( L % 2 == 0 , "ylm_lanes: even lane count" );
      typedef double v2 __attribute__(( vector_size( 16 ), aligned( 8 ) ));
      constexpr int W = L / 2;
      auto V = [](const double* p) { return reinterpret_cast<const v2*>( p ); };
      auto M = [](double* p) { return reinterpret_cast<v2*>( p ); };
      const v2* ux = V( ux_ ); const v2* uy = V( uy_ ); const v2* uz = V( uz_ ); const v2* ir = V( inv_r_ );
      v2* yr = M( yr_ ); v2* yi = M( yi_ ); v2* gr = M( gr_ ); v2* gi = M( gi_ );
      const v2 one = { 1.0, 1.0 }, zero = { 0.0, 0.0 };
      for(int w=0;w<W;w++) { yr[w] = one; yi[w] = zero; }
      if( GRAD ) for(int k=0;k<3*W;k++) { gr[k] = zero; gi[k] = zero; }
      for(int m=0;m<=lmax;m++)
      {
        const int d = lmh(m,m);
        if( m > 0 )
        {
          const int p = lmh(m-1,m-1);
          const double a = A[d];
          for(int w=0;w<W;w++)
          {
            const v2 pr = yr[p*W+w], pi = yi[p*W+w];
            yr[d*W+w] = a * ( ux[w] * pr - uy[w] * pi );
            yi[d*W+w] = a * ( ux[w] * pi + uy[w] * pr );
            if( GRAD )
            {
              for(int k=0;k<3;k++)
              {
                const v2 qr = gr[(3*p+k)*W+w], qi = gi[(3*p+k)*W+w];
                gr[(3*d+k)*W+w] = a * ( ux[w] * qr - uy[w] * qi );
                gi[(3*d+k)*W+w] = a * ( ux[w] * qi + uy[w] * qr );
              }
              gr[(3*d+0)*W+w] += a * pr;  gi[(3*d+0)*W+w] += a * pi;
              gr[(3*d+1)*W+w] -= a * pi;  gi[(3*d+1)*W+w] += a * pr;
            }
          }
        }
        for(int l=m+1;l<=lmax;l++)
        {
          const int q = lmh(l,m), q1 = lmh(l-1,m);
          const double a = A[q], bb = B[q];
          const bool two = l - 2 >= m;
          const int q2 = two ? lmh(l-2,m) : q1;
          const double b2 = two ? bb : 0.0;
          for(int w=0;w<W;w++)
          {
            yr[q*W+w] = a * uz[w] * yr[q1*W+w] - b2 * yr[q2*W+w];
            yi[q*W+w] = a * uz[w] * yi[q1*W+w] - b2 * yi[q2*W+w];
            if( GRAD )
            {
              for(int k=0;k<3;k++)
              {
                gr[(3*q+k)*W+w] = a * uz[w] * gr[(3*q1+k)*W+w] - b2 * gr[(3*q2+k)*W+w];
                gi[(3*q+k)*W+w] = a * uz[w] * gi[(3*q1+k)*W+w] - b2 * gi[(3*q2+k)*W+w];
              }
              gr[(3*q+2)*W+w] += a * yr[q1*W+w];
              gi[(3*q+2)*W+w] += a * yi[q1*W+w];
            }
          }
        }
      }
      if( GRAD )
      {
        const int n = lmh(lmax,lmax) + 1;
        for(int q=0;q<n;q++)
        {
          v2* x = gr + 3*q*W;
          v2* y = gi + 3*q*W;
          for(int w=0;w<W;w++)
          {
            const v2 sr = ux[w] * x[w] + uy[w] * x[W+w] + uz[w] * x[2*W+w];
            const v2 si = ux[w] * y[w] + uy[w] * y[W+w] + uz[w] * y[2*W+w];
            x[w] = ( x[w] - sr * ux[w] ) * ir[w]; x[W+w] = ( x[W+w] - sr * uy[w] ) * ir[w]; x[2*W+w] = ( x[2*W+w] - sr * uz[w] ) * ir[w];
            y[w] = ( y[w] - si * ux[w] ) * ir[w]; y[W+w] = ( y[W+w] - si * uy[w] ) * ir[w]; y[2*W+w] = ( y[2*W+w] - si * uz[w] ) * ir[w];
          }
        }
      }
    }

    // Y_lm recurrence coefficients (host): A[lmh(m,m)] = -sqrt((2m+1)/(2m)), A[lmh(l,m)], B[lmh(l,m)] of the
    // normalised associated Legendre recurrence in cos(theta)
    inline void ylm_coefficients(int lmax, double* A, double* B)
    {
      A[0] = 1.0; B[0] = 0.0;
      for(int l=1;l<=lmax;l++)
      {
        for(int m=0;m<=l;m++)
        {
          const int q = lmh(l,m);
          const double L = l, M = m;
          if( m == l ) { A[q] = -std::sqrt( ( 2.0 * L + 1.0 ) / ( 2.0 * L ) ); B[q] = 0.0; continue; }
          A[q] = std::sqrt( ( 2.0 * L - 1.0 ) * ( 2.0 * L + 1.0 ) / ( ( L - M ) * ( L + M ) ) );
          B[q] = ( l - m >= 2 ) ? std::sqrt( ( 2.0 * L + 1.0 ) * ( L + M - 1.0 ) * ( L - M - 1.0 ) / ( ( 2.0 * L - 3.0 ) * ( L - M ) * ( L + M ) ) ) : 0.0;
        }
      }
    }

  }
}
