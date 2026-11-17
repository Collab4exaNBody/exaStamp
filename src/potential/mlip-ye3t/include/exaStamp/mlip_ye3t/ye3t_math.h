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

// Device / host math of the YE3T kernels: ChebExpCos radial base, uniform cubic splines (spherical harmonics:
// mlip_ylm.h).

#include <exaStamp/mlip_ye3t/ye3t_tables.h>
#include <exaStamp/mlip_utils/mlip_ylm.h>

#include <onika/cuda/cuda.h>
#include <cmath>

namespace exaStamp
{
  namespace mlip_ye3t
  {
    using mlip_utils::lmh;
    using mlip_utils::ylm;
    using mlip_utils::ylm_coefficients;

    // radial base g_k(r), dg_k/dr, k < b.nbase (r < b.rc)
    ONIKA_HOST_DEVICE_FUNC inline void cheb_exp_cos( const Ye3tBond& b, double r, double* g, double* dg )
    {
      const double pi = M_PI;
      const double rc = b.rc;
      const double s = r / rc;
      double fc = 0.5 * ( 1.0 + cos( pi * s ) );
      double dfc = -0.5 * pi * sin( pi * s ) / rc;
      if( b.dcut > 0.0 && r > rc - b.dcut )
      {
        const double ph = pi * ( r - ( rc - b.dcut ) ) / b.dcut;
        const double sw = 0.5 * ( 1.0 + cos( ph ) );
        const double dsw = -0.5 * pi * sin( ph ) / b.dcut;
        dfc = dfc * sw + fc * dsw;
        fc = fc * sw;
      }
      g[0] = fc;
      dg[0] = dfc;
      if( b.nbase == 1 ) return;
      const double lam = b.lambda;
      const double e = exp( -lam * ( s - 1.0 ) );
      const double den = expm1( lam );
      const double x = 1.0 - 2.0 * ( e - 1.0 ) / den;
      const double dx = 2.0 * lam * e / ( den * rc );
      double t0 = 1.0, t1 = x, d0 = 0.0, d1 = dx;
      for(int k=1;k<b.nbase;k++)
      {
        if( k > 1 )
        {
          const double t2 = 2.0 * x * t1 - t0;
          const double d2 = 2.0 * ( dx * t1 + x * d1 ) - d0;
          t0 = t1; t1 = t2; d0 = d1; d1 = d2;
        }
        const double v = 0.5 * ( 1.0 - t1 ), dv = -0.5 * d1;
        g[k] = v * fc;
        dg[k] = dv * fc + v * dfc;
      }
    }

    // location of r in the spline table of a bond: interval k (>= 1), local coordinate; false below the table
    struct SplinePos
    {
      int k = 0;
      double x = 0.0, x2 = 0.0, x3 = 0.0;
      double scale = 0.0;
      ONIKA_HOST_DEVICE_FUNC inline bool locate(const Ye3tBond& b, double r)
      {
        scale = b.scale;
        const double t = r * scale;
        k = int( t );
        if( k <= 0 ) return false;
        if( k > b.nint ) k = b.nint;
        x = t - k;
        x2 = x * x;
        x3 = x2 * x;
        return true;
      }
      // value, derivative of function f of a table with nf functions
      ONIKA_HOST_DEVICE_FUNC inline void eval(const double* spl, int nf, int f, double& v, double& d) const
      {
        const double* c = spl + ( size_t(k) * nf + f ) * 4;
        v = c[0] + c[1] * x + c[2] * x2 + c[3] * x3;
        d = ( c[1] + 2.0 * c[2] * x + 3.0 * c[3] * x2 ) * scale;
      }
    };

  }
}
