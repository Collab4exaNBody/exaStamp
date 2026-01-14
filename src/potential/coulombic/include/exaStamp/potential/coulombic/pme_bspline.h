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

#include <onika/cuda/cuda.h>
#include <cmath>

namespace exaStamp
{
  //============================================================================
  // Cardinal B-Spline Functions
  //============================================================================
  //
  // The cardinal B-spline M_n(u) of order n is defined recursively:
  //   M_1(u) = 1 if 0 <= u < 1, else 0
  //   M_n(u) = u/(n-1) * M_{n-1}(u) + (n-u)/(n-1) * M_{n-1}(u-1)
  //
  // For PME, we need M_n(u) evaluated at fractional grid positions.
  // The spline is non-zero only for u in [0, n).
  //
  //============================================================================

  //============================================================================
  // Compute B-spline coefficients for a single fractional coordinate
  // 
  // Input:  u = fractional coordinate (0 <= u < 1 typically)
  //         order = spline order (4, 5, or 6)
  // Output: theta[0..order-1] = M_order(u-j) for j = 0..order-1
  //         dtheta[0..order-1] = dM_order/du (derivatives)
  //============================================================================
  
  template<int ORDER>
  ONIKA_HOST_DEVICE_FUNC
  inline void compute_bspline_coeffs(
    double u,           // Fractional position within cell (0 to 1)
    double* theta,      // Output: spline values [ORDER]
    double* dtheta)     // Output: spline derivatives [ORDER]
  {
    // Initialize spline values using the recursion
    // Start with order 1 (box function)
    theta[ORDER - 1] = 0.0;
    theta[0] = 1.0 - u;
    theta[1] = u;
    
    // Build up spline using recursion
    for (int k = 2; k < ORDER - 1; ++k)
    {
      double div = 1.0 / static_cast<double>(k);
      theta[k] = u * theta[k-1] * div;
      
      for (int j = k - 1; j > 0; --j)
      {
        theta[j] = ((u + j) * theta[j-1] + (k + 1 - j - u) * theta[j]) * div;
      }
      
      theta[0] = (1.0 - u) * theta[0] * div;
    }
    
    // Compute derivatives before final recursion step
    // dM_n(u)/du = M_{n-1}(u) - M_{n-1}(u-1)
    dtheta[0] = -theta[0];
    for (int j = 1; j < ORDER; ++j)
    {
      dtheta[j] = theta[j-1] - theta[j];
    }
    
    // Final recursion step for spline values
    {
      double div = 1.0 / static_cast<double>(ORDER - 1);
      theta[ORDER - 1] = u * theta[ORDER - 2] * div;
      
      for (int j = ORDER - 2; j > 0; --j)
      {
        theta[j] = ((u + j) * theta[j-1] + (ORDER - j - u) * theta[j]) * div;
      }
      
      theta[0] = (1.0 - u) * theta[0] * div;
    }
  }

  //============================================================================
  // Specialized versions for common spline orders (compile-time optimization)
  //============================================================================
  
  // Order 4 (Cubic B-spline) - Most common choice
  ONIKA_HOST_DEVICE_FUNC
  inline void compute_bspline_order4(double u, double* theta, double* dtheta)
  {
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double one_minus_u = 1.0 - u;
    const double one_minus_u2 = one_minus_u * one_minus_u;
    const double one_minus_u3 = one_minus_u2 * one_minus_u;
    
    // M_4(u), M_4(u-1), M_4(u-2), M_4(u-3)
    // Evaluated at shifted positions for the 4 grid points
    const double c1_6 = 1.0 / 6.0;
    const double c2_3 = 2.0 / 3.0;
    
    theta[0] = c1_6 * one_minus_u3;
    theta[1] = c2_3 - u2 + 0.5 * u3;
    theta[2] = c2_3 - one_minus_u2 + 0.5 * one_minus_u3;
    theta[3] = c1_6 * u3;
    
    // Derivatives
    dtheta[0] = -0.5 * one_minus_u2;
    dtheta[1] = -2.0 * u + 1.5 * u2;
    dtheta[2] = 2.0 * one_minus_u - 1.5 * one_minus_u2;
    dtheta[3] = 0.5 * u2;
  }

  // Order 5 (Quartic B-spline)
  ONIKA_HOST_DEVICE_FUNC
  inline void compute_bspline_order5(double u, double* theta, double* dtheta)
  {
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double u4 = u3 * u;
    
    const double t = 1.0 - u;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double t4 = t3 * t;
    
    const double c1_24 = 1.0 / 24.0;
    const double c1_6 = 1.0 / 6.0;
    const double c11_24 = 11.0 / 24.0;
    
    theta[0] = c1_24 * t4;
    theta[1] = c1_24 * (1.0 + 4.0*t + 6.0*t2 + 4.0*t3 - 4.0*t4);
    theta[2] = c1_24 * (11.0 + 4.0*(u - t) + 6.0*(u2 - t2) - 4.0*(u3 + t3) + (u4 + t4));
    theta[3] = c1_24 * (1.0 + 4.0*u + 6.0*u2 + 4.0*u3 - 4.0*u4);
    theta[4] = c1_24 * u4;
    
    // Derivatives
    dtheta[0] = -c1_6 * t3;
    dtheta[1] = c1_6 * (-1.0 - 3.0*t - 3.0*t2 + 4.0*t3);
    dtheta[2] = c1_6 * (1.0 + 3.0*(u - t) - 3.0*(u2 + t2) + (u3 - t3));
    dtheta[3] = c1_6 * (1.0 + 3.0*u + 3.0*u2 - 4.0*u3);
    dtheta[4] = c1_6 * u3;
  }

  // Order 6 (Quintic B-spline) - Higher accuracy
  ONIKA_HOST_DEVICE_FUNC
  inline void compute_bspline_order6(double u, double* theta, double* dtheta)
  {
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double u4 = u3 * u;
    const double u5 = u4 * u;
    
    const double t = 1.0 - u;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double t4 = t3 * t;
    const double t5 = t4 * t;
    
    const double c1_120 = 1.0 / 120.0;
    const double c1_24 = 1.0 / 24.0;
    const double c26_120 = 26.0 / 120.0;
    const double c66_120 = 66.0 / 120.0;
    
    theta[0] = c1_120 * t5;
    theta[1] = c1_120 * (1.0 + 5.0*t + 10.0*t2 + 10.0*t3 + 5.0*t4 - 5.0*t5);
    theta[2] = c1_120 * (26.0 + 10.0*(u - t) + 10.0*(t2 - u2) + 5.0*(u3 - t3) - 5.0*(u4 - t4) + (t5 - u5));
    theta[3] = c1_120 * (26.0 + 10.0*(t - u) + 10.0*(u2 - t2) + 5.0*(t3 - u3) - 5.0*(t4 - u4) + (u5 - t5));
    theta[4] = c1_120 * (1.0 + 5.0*u + 10.0*u2 + 10.0*u3 + 5.0*u4 - 5.0*u5);
    theta[5] = c1_120 * u5;
    
    // Derivatives
    dtheta[0] = -c1_24 * t4;
    dtheta[1] = c1_24 * (-1.0 - 4.0*t - 6.0*t2 - 4.0*t3 + 5.0*t4);
    dtheta[2] = c1_24 * (2.0 - 4.0*(u + t) + 3.0*(u2 - t2) + 4.0*(u3 + t3) - (u4 + t4));
    dtheta[3] = c1_24 * (-2.0 + 4.0*(u + t) + 3.0*(t2 - u2) - 4.0*(u3 + t3) + (u4 + t4));
    dtheta[4] = c1_24 * (1.0 + 4.0*u + 6.0*u2 + 4.0*u3 - 5.0*u4);
    dtheta[5] = c1_24 * u4;
  }

  //============================================================================
  // Generic dispatcher based on runtime order
  //============================================================================
  
  ONIKA_HOST_DEVICE_FUNC
  inline void compute_bspline(int order, double u, double* theta, double* dtheta)
  {
    switch (order)
    {
      case 4:
        compute_bspline_order4(u, theta, dtheta);
        break;
      case 5:
        compute_bspline_order5(u, theta, dtheta);
        break;
      case 6:
        compute_bspline_order6(u, theta, dtheta);
        break;
      default:
        // Fallback to generic (slower)
        compute_bspline_coeffs<8>(u, theta, dtheta); // Max supported order
        break;
    }
  }

  //============================================================================
  // Compute B-spline modulus for influence function correction
  //
  // The B-spline modulus corrects for the discrete interpolation:
  // |b(m)|² where b(m) = Σ_{k=0}^{order-1} M_order(k+1) * exp(2πi*m*k/M)
  //============================================================================
  
  inline void compute_bspline_moduli(
    int order,
    int mesh_size,
    double* moduli)  // Output array of size mesh_size
  {
    // Compute spline values at integer points
    double theta[8];  // Max order = 8
    double dtheta[8]; // Not used, but needed for function call
    
    // M_n(k) for k = 1, 2, ..., n-1 (non-zero values)
    double spline_vals[8];
    for (int k = 0; k < order - 1; ++k)
    {
      // M_order is non-zero at integer points 1, 2, ..., order-1
      double u = static_cast<double>(k + 1) / static_cast<double>(order);
      compute_bspline(order, u, theta, dtheta);
      spline_vals[k] = theta[k];
    }
    
    // Compute modulus for each frequency
    for (int m = 0; m < mesh_size; ++m)
    {
      double sum_cos = 0.0;
      double sum_sin = 0.0;
      
      const double arg = 2.0 * M_PI * m / mesh_size;
      
      for (int k = 0; k < order; ++k)
      {
        // Use the recursive definition to get proper values
        double frac = static_cast<double>(k) / static_cast<double>(order);
        compute_bspline(order, frac, theta, dtheta);
        
        double phase = arg * k;
        sum_cos += theta[0] * cos(phase);
        sum_sin += theta[0] * sin(phase);
      }
      
      // |b(m)|² 
      double mod_sq = sum_cos * sum_cos + sum_sin * sum_sin;
      
      // Avoid division by zero at m=0
      moduli[m] = (mod_sq > 1e-30) ? mod_sq : 1.0;
    }
    
    // Normalize so that moduli[0] = 1
    double norm = moduli[0];
    for (int m = 0; m < mesh_size; ++m)
    {
      moduli[m] /= norm;
    }
  }

  //============================================================================
  // Alternative: Analytical B-spline moduli for common orders
  //============================================================================
  
  inline void compute_bspline_moduli_analytical(
    int order,
    int mesh_size,
    double* moduli)
  {
    for (int m = 0; m < mesh_size; ++m)
    {
      if (m == 0)
      {
        moduli[m] = 1.0;
        continue;
      }
      
      // General formula: |b(m)|² = |sin(πm/M) / (πm/M)|^(2*order)
      // This is the Fourier transform of the B-spline
      double arg = M_PI * m / mesh_size;
      double sinc = sin(arg) / arg;
      moduli[m] = pow(sinc, 2 * order);
    }
  }

}
