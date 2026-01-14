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

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_yaml.h>
#include <cmath>
#include <array>
#include <cufft.h>

#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>

namespace exaStamp
{
  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  //============================================================================
  // Physical Constants (same as ewald.h for consistency)
  //============================================================================
  namespace pme_constants
  {
    static constexpr double epsilonZero = 0x1.337f14f782bb5p-21;
    static constexpr double fpe0 = 0x1.e303a50b0d1ap-18;
    static constexpr double ONE_4PI_EPS0 = 0.25 / M_PI / epsilonZero;
  }

  //============================================================================
  // PME Grid Dimensions
  //============================================================================
  struct PMEGridSize
  {
    int nx = 64;   // Mesh points in x
    int ny = 64;   // Mesh points in y
    int nz = 64;   // Mesh points in z
    
    ONIKA_HOST_DEVICE_FUNC
    inline int total() const { return nx * ny * nz; }
    
    ONIKA_HOST_DEVICE_FUNC
    inline int total_complex() const { return nx * ny * (nz/2 + 1); }
  };

  //============================================================================
  // PME Parameters
  //============================================================================
  struct alignas(DEFAULT_ALIGNMENT) PMEParameters
  {
    // Grid dimensions
    PMEGridSize grid_size;
    
    // B-spline order (4 = cubic, 5 = quartic, 6 = quintic)
    int spline_order = 4;
    
    // Ewald splitting parameter (alpha)
    double alpha = 0.0;
    
    // Real-space cutoff
    double rcut = 0.0;
    
    // Accuracy target
    double accuracy_relative = 1.0e-5;
    
    // Domain information
    double volume = 0.0;
    Vec3d box_size = {0.0, 0.0, 0.0};
    Vec3d inv_box_size = {0.0, 0.0, 0.0};
    
    // Grid spacing
    Vec3d grid_spacing = {0.0, 0.0, 0.0};
    Vec3d inv_grid_spacing = {0.0, 0.0, 0.0};
    
    // Precomputed influence function (B-spline modulated Green's function)
    // Size: nx * ny * (nz/2 + 1) for R2C transform
    onika::memory::CudaMMVector<double> influence_function;
    
    // B-spline modulation factors for correction
    onika::memory::CudaMMVector<double> bspline_moduli_x;
    onika::memory::CudaMMVector<double> bspline_moduli_y;
    onika::memory::CudaMMVector<double> bspline_moduli_z;
    
    // erfc approximation constants (for short-range)
    double EWALD_P = 0.3275911;
    double A1 = 0.254829592;
    double A2 = -0.284496736;
    double A3 = 1.421413741;
    double A4 = -1.453152027;
    double A5 = 1.061405429;
  };

  //============================================================================
  // Read-only PME parameters for GPU kernels
  //============================================================================
  struct ReadOnlyPMEParameters
  {
    // Grid dimensions
    int nx, ny, nz;
    int spline_order;
    
    // Physical parameters
    double alpha;
    double volume;
    
    // Grid information
    double inv_box_x, inv_box_y, inv_box_z;
    double grid_spacing_x, grid_spacing_y, grid_spacing_z;
    
    // Precomputed arrays (device pointers)
    const double* __restrict__ influence_function;
    const double* __restrict__ bspline_moduli_x;
    const double* __restrict__ bspline_moduli_y;
    const double* __restrict__ bspline_moduli_z;
    
    // erfc constants
    double EWALD_P, A1, A2, A3, A4, A5;
    
    ReadOnlyPMEParameters() = default;
    ReadOnlyPMEParameters(const ReadOnlyPMEParameters&) = default;
    
    inline ReadOnlyPMEParameters(const PMEParameters& p)
      : nx(p.grid_size.nx)
      , ny(p.grid_size.ny)
      , nz(p.grid_size.nz)
      , spline_order(p.spline_order)
      , alpha(p.alpha)
      , volume(p.volume)
      , inv_box_x(p.inv_box_size.x)
      , inv_box_y(p.inv_box_size.y)
      , inv_box_z(p.inv_box_size.z)
      , grid_spacing_x(p.grid_spacing.x)
      , grid_spacing_y(p.grid_spacing.y)
      , grid_spacing_z(p.grid_spacing.z)
      , influence_function(p.influence_function.data())
      , bspline_moduli_x(p.bspline_moduli_x.data())
      , bspline_moduli_y(p.bspline_moduli_y.data())
      , bspline_moduli_z(p.bspline_moduli_z.data())
      , EWALD_P(p.EWALD_P)
      , A1(p.A1)
      , A2(p.A2)
      , A3(p.A3)
      , A4(p.A4)
      , A5(p.A5)
    {}
  };

  //============================================================================
  // PME Mesh Data (working storage)
  //============================================================================
  struct PMEMesh
  {
    // Real-space charge mesh
    onika::memory::CudaMMVector<double> charge_mesh;
    
    // Complex k-space mesh (after R2C FFT)
    // Using cufftDoubleComplex for compatibility with cuFFT
    onika::memory::CudaMMVector<cufftDoubleComplex> kspace_mesh;
    
    // Force meshes (x, y, z components)
    onika::memory::CudaMMVector<double> force_mesh_x;
    onika::memory::CudaMMVector<double> force_mesh_y;
    onika::memory::CudaMMVector<double> force_mesh_z;
    
    // cuFFT plans
    cufftHandle fft_plan_forward;   // R2C
    cufftHandle fft_plan_inverse;   // C2R
    bool plans_initialized = false;
    
    // Grid dimensions (cached)
    int nx = 0, ny = 0, nz = 0;
    
    void resize(int nx_, int ny_, int nz_)
    {
      nx = nx_; ny = ny_; nz = nz_;
      
      const size_t real_size = nx * ny * nz;
      const size_t complex_size = nx * ny * (nz/2 + 1);
      
      charge_mesh.resize(real_size);
      kspace_mesh.resize(complex_size);
      force_mesh_x.resize(real_size);
      force_mesh_y.resize(real_size);
      force_mesh_z.resize(real_size);
    }
    
    void init_fft_plans()
    {
      if (plans_initialized) return;
      
      // Create R2C (Real to Complex) forward plan
      cufftPlan3d(&fft_plan_forward, nx, ny, nz, CUFFT_D2Z);
      
      // Create C2R (Complex to Real) inverse plan
      cufftPlan3d(&fft_plan_inverse, nx, ny, nz, CUFFT_Z2D);
      
      plans_initialized = true;
    }
    
    void destroy_fft_plans()
    {
      if (!plans_initialized) return;
      
      cufftDestroy(fft_plan_forward);
      cufftDestroy(fft_plan_inverse);
      plans_initialized = false;
    }
    
    ~PMEMesh()
    {
      destroy_fft_plans();
    }
  };

  //============================================================================
  // Per-particle B-spline data (precomputed for efficiency)
  //============================================================================
  struct ParticleBSplineData
  {
    // Fractional coordinates (0 to grid_size)
    onika::memory::CudaMMVector<double> frac_x;
    onika::memory::CudaMMVector<double> frac_y;
    onika::memory::CudaMMVector<double> frac_z;
    
    // Base grid indices
    onika::memory::CudaMMVector<int> grid_idx_x;
    onika::memory::CudaMMVector<int> grid_idx_y;
    onika::memory::CudaMMVector<int> grid_idx_z;
    
    // B-spline coefficients: [n_particles][spline_order]
    // Stored as flattened arrays
    onika::memory::CudaMMVector<double> theta_x;  // M_n(frac_x)
    onika::memory::CudaMMVector<double> theta_y;  // M_n(frac_y)
    onika::memory::CudaMMVector<double> theta_z;  // M_n(frac_z)
    
    // B-spline derivatives (for force calculation)
    onika::memory::CudaMMVector<double> dtheta_x; // dM_n/du
    onika::memory::CudaMMVector<double> dtheta_y;
    onika::memory::CudaMMVector<double> dtheta_z;
    
    size_t n_particles = 0;
    int spline_order = 0;
    
    void resize(size_t n, int order)
    {
      n_particles = n;
      spline_order = order;
      
      frac_x.resize(n);
      frac_y.resize(n);
      frac_z.resize(n);
      
      grid_idx_x.resize(n);
      grid_idx_y.resize(n);
      grid_idx_z.resize(n);
      
      const size_t coeff_size = n * order;
      theta_x.resize(coeff_size);
      theta_y.resize(coeff_size);
      theta_z.resize(coeff_size);
      
      dtheta_x.resize(coeff_size);
      dtheta_y.resize(coeff_size);
      dtheta_z.resize(coeff_size);
    }
  };

  //============================================================================
  // Short-range energy/force computation (same as Ewald)
  //============================================================================
  template<class PMEParamsT>
  ONIKA_HOST_DEVICE_FUNC static inline void pme_compute_short_range(const PMEParamsT& p, double c, double r, double& e, double& de)
  {
    const double cf = pme_constants::ONE_4PI_EPS0 * c;
    const double grij = p.alpha * r;
    const double expm2 = exp(-grij * grij);
    const double t = 1.0 / (1.0 + p.EWALD_P * grij);
    const double erfc = t * (p.A1 + t * (p.A2 + t * (p.A3 + t * (p.A4 + t * p.A5)))) * expm2;
    const double invr = 1. / r;
    e = cf * erfc * invr;
    de = - (cf * 2. * p.alpha / std::sqrt(M_PI) * expm2 + e) * invr;
    
    // // Energy: q1*q2 * erfc(alpha*r) / r
    // const double inv_r = 1.0 / r;
    // e = pme_constants::ONE_4PI_EPS0 * q1q2 * erfc_val * inv_r;
    
    // // Force: -dE/dr = q1*q2 * [erfc(alpha*r)/r² + 2*alpha/sqrt(pi) * exp(-alpha²r²)/r]
    // const double two_alpha_sqrtpi = 2.0 * alpha / std::sqrt(M_PI);
    // de = -pme_constants::ONE_4PI_EPS0 * q1q2 * (erfc_val * inv_r + two_alpha_sqrtpi * exp_term) * inv_r;
  }

  //============================================================================
  // PME Parameter Initialization (full implementation)
  //============================================================================
  inline void pme_init_parameters(
    const Vec3d& box_size,
    const PMEGridSize& grid_size,
    int spline_order,
    double alpha,
    double rcut,
    double accuracy_relative,
    uint64_t natoms,
    double qsq,
    PMEParameters& p,
    ::exanb::LogStreamWrapper& ldbg)
  {
    p.grid_size = grid_size;
    p.spline_order = spline_order;
    p.rcut = rcut;
    p.accuracy_relative = accuracy_relative;
    
    p.box_size = box_size;
    p.volume = box_size.x * box_size.y * box_size.z;
    
    if (p.volume == 0.0) return;
    
    p.inv_box_size = Vec3d{1.0/box_size.x, 1.0/box_size.y, 1.0/box_size.z};
    
    p.grid_spacing = Vec3d{
      box_size.x / grid_size.nx,
      box_size.y / grid_size.ny,
      box_size.z / grid_size.nz
    };
    
    p.inv_grid_spacing = Vec3d{
      static_cast<double>(grid_size.nx) / box_size.x,
      static_cast<double>(grid_size.ny) / box_size.y,
      static_cast<double>(grid_size.nz) / box_size.z
    };
    
    // Auto-compute alpha if not specified
    if (alpha <= 0.0)
    {
      // Optimal alpha balances real-space and k-space errors
      // Heuristic: alpha = 5 / rcut for typical accuracy
      // More sophisticated: alpha = sqrt(-log(accuracy)) / rcut
      double beta = -std::log(accuracy_relative);
      alpha = std::sqrt(beta) / rcut;
      
      // Clamp to reasonable range
      alpha = std::max(0.1 / rcut, std::min(10.0 / rcut, alpha));
    }
    p.alpha = alpha;
    
    ldbg << "PME Parameters:" << std::endl;
    ldbg << "  Grid size: " << grid_size.nx << " x " << grid_size.ny << " x " << grid_size.nz << std::endl;
    ldbg << "  Spline order: " << spline_order << std::endl;
    ldbg << "  Alpha (Ewald parameter): " << alpha << std::endl;
    ldbg << "  Real-space cutoff: " << rcut << std::endl;
    ldbg << "  Grid spacing: " << p.grid_spacing.x << " x " << p.grid_spacing.y << " x " << p.grid_spacing.z << std::endl;
    
    // Compute B-spline moduli for influence function correction
    const int nx = grid_size.nx;
    const int ny = grid_size.ny;
    const int nz = grid_size.nz;
    
    p.bspline_moduli_x.resize(nx);
    p.bspline_moduli_y.resize(ny);
    p.bspline_moduli_z.resize(nz);
    
    // Use analytical formula for B-spline moduli
    // |b(m)|² = |sin(πm/M) / (πm/M)|^(2*order)
    auto compute_moduli = [](int order, int mesh_size, double* moduli)
    {
      for (int m = 0; m < mesh_size; ++m)
      {
        if (m == 0)
        {
          moduli[m] = 1.0;
        }
        else
        {
          double arg = M_PI * m / mesh_size;
          double sinc = std::sin(arg) / arg;
          moduli[m] = std::pow(sinc, 2 * order);
        }
      }
    };
    
    compute_moduli(spline_order, nx, p.bspline_moduli_x.data());
    compute_moduli(spline_order, ny, p.bspline_moduli_y.data());
    compute_moduli(spline_order, nz, p.bspline_moduli_z.data());
    
    // Compute influence function (precomputed for efficiency)
    // For R2C transform, only store positive kz frequencies
    const int nz_complex = nz / 2 + 1;
    const size_t influence_size = static_cast<size_t>(nx) * ny * nz_complex;
    p.influence_function.resize(influence_size);
    
    const double inv_vol = 1.0 / p.volume;
    const double alpha_sq = alpha * alpha;
    const double four_alpha_sq = 4.0 * alpha_sq;
    
    // Reciprocal lattice vectors
    const double bx = 2.0 * M_PI / box_size.x;
    const double by = 2.0 * M_PI / box_size.y;
    const double bz = 2.0 * M_PI / box_size.z;
    
    // Coulomb constant
    constexpr double ONE_4PI_EPS0 = 0.25 / M_PI / pme_constants::epsilonZero;
    
    size_t idx = 0;
    for (int ix = 0; ix < nx; ++ix)
    {
      // k_x index: 0, 1, ..., nx/2, -nx/2+1, ..., -1
      int mx = (ix <= nx/2) ? ix : ix - nx;
      double kx = mx * bx;
      double kx_sq = kx * kx;
      double bmod_x = p.bspline_moduli_x[ix];
      
      for (int iy = 0; iy < ny; ++iy)
      {
        int my = (iy <= ny/2) ? iy : iy - ny;
        double ky = my * by;
        double ky_sq = ky * ky;
        double bmod_xy = bmod_x * p.bspline_moduli_y[iy];
        
        for (int iz = 0; iz < nz_complex; ++iz)
        {
          int mz = iz;  // Only positive frequencies for R2C
          double kz = mz * bz;
          double kz_sq = kz * kz;
          
          double k_sq = kx_sq + ky_sq + kz_sq;
          
          if (k_sq < 1e-12)
          {
            // k = 0 term: set to zero (net charge assumed neutral)
            p.influence_function[idx] = 0.0;
          }
          else
          {
            // G(k) = (1/πV) * exp(-k²/4α²) / k² / |b(k)|²
            double bmod = bmod_xy * p.bspline_moduli_z[iz];
            
            // Avoid division by very small moduli
            if (bmod < 1e-10) bmod = 1e-10;
            
            double exp_term = std::exp(-k_sq / four_alpha_sq);
            double G = inv_vol * exp_term / (k_sq * bmod) * ONE_4PI_EPS0;
            
            p.influence_function[idx] = G;
          }
          
          ++idx;
        }
      }
    }
    
    ldbg << "  Influence function computed (" << influence_size << " k-points)" << std::endl;
    
    // Compute some statistics for logging
    double G_max = 0.0, G_sum = 0.0;
    for (size_t i = 0; i < influence_size; ++i)
    {
      double G = std::abs(p.influence_function[i]);
      G_max = std::max(G_max, G);
      G_sum += G;
    }
    
    ldbg << "  Influence function: max = " << G_max << ", sum = " << G_sum << std::endl;
  }

  // Overload without debug stream
  inline void pme_init_parameters(
    const Vec3d& box_size,
    const PMEGridSize& grid_size,
    int spline_order,
    double alpha,
    double rcut,
    double accuracy_relative,
    uint64_t natoms,
    double qsq,
    PMEParameters& p)
  {
    pme_init_parameters(box_size, grid_size, spline_order, alpha, rcut, 
                        accuracy_relative, natoms, qsq, p, ::exanb::ldbg << "");
  }

}

//============================================================================
// YAML conversion for PMEParameters
//============================================================================
namespace YAML
{
  using exaStamp::PMEParameters;
  using exaStamp::PMEGridSize;
  using onika::physics::Quantity;
  using exanb::Vec3d;

  template<> struct convert<PMEGridSize>
  {
    static bool decode(const Node& node, PMEGridSize& v)
    {
      if (node.IsSequence() && node.size() == 3)
      {
        v.nx = node[0].as<int>();
        v.ny = node[1].as<int>();
        v.nz = node[2].as<int>();
        return true;
      }
      if (node.IsScalar())
      {
        int n = node.as<int>();
        v.nx = v.ny = v.nz = n;
        return true;
      }
      return false;
    }
  };

  template<> struct convert<PMEParameters>
  {
    static bool decode(const Node& node, PMEParameters& v)
    {
      if (!node.IsMap()) return false;
      
      PMEGridSize grid_size;
      if (node["mesh_size"])
        grid_size = node["mesh_size"].as<PMEGridSize>();
      
      int spline_order = 4;
      if (node["spline_order"])
        spline_order = node["spline_order"].as<int>();
      
      double alpha = 0.0;
      if (node["alpha"])
        alpha = node["alpha"].as<Quantity>().convert();
      
      double rcut = 10.0;
      if (node["rcut"])
        rcut = node["rcut"].as<Quantity>().convert();
      
      double accuracy = 1.0e-5;
      if (node["accuracy_relative"])
        accuracy = node["accuracy_relative"].as<Quantity>().convert();
      
      Vec3d box_size = node["size"].as<Vec3d>();
      
      exaStamp::pme_init_parameters(box_size, grid_size, spline_order, 
                                     alpha, rcut, accuracy, 1, 1.0, v);
      return true;
    }
  };
}
