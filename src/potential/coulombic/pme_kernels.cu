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

#include <exaStamp/potential/coulombic/pme.h>
#include <exaStamp/potential/coulombic/pme_bspline.h>
#include <exaStamp/potential/coulombic/pme_kernels.h>
#include <onika/cuda/cuda.h>
#include <cufft.h>

namespace exaStamp
{
  //============================================================================
  // Constants
  //============================================================================
  
  static constexpr int PME_MAX_ORDER = 8;
  static constexpr int PME_BLOCK_SIZE = 256;

  //============================================================================
  // KERNEL 1: Compute B-spline coefficients for all particles
  //============================================================================
  
  __global__ void pme_compute_bsplines_kernel(
    const double* __restrict__ rx,    // Particle x positions
    const double* __restrict__ ry,    // Particle y positions  
    const double* __restrict__ rz,    // Particle z positions
    int n_particles,
    int spline_order,
    int nx, int ny, int nz,           // Grid dimensions
    double inv_box_x, double inv_box_y, double inv_box_z,
    // Outputs:
    double* __restrict__ frac_x,
    double* __restrict__ frac_y,
    double* __restrict__ frac_z,
    int* __restrict__ grid_idx_x,
    int* __restrict__ grid_idx_y,
    int* __restrict__ grid_idx_z,
    double* __restrict__ theta_x,     // [n_particles * spline_order]
    double* __restrict__ theta_y,
    double* __restrict__ theta_z,
    double* __restrict__ dtheta_x,
    double* __restrict__ dtheta_y,
    double* __restrict__ dtheta_z)
  {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;
    
    // Get particle position
    double px = rx[i];
    double py = ry[i];
    double pz = rz[i];
    
    // Convert to fractional coordinates [0, grid_size)
    // First get fractional box coordinates [0, 1)
    double fx = px * inv_box_x;
    double fy = py * inv_box_y;
    double fz = pz * inv_box_z;
    
    // Handle periodic boundaries
    fx -= floor(fx);
    fy -= floor(fy);
    fz -= floor(fz);
    
    // Scale to grid coordinates
    fx *= nx;
    fy *= ny;
    fz *= nz;
    
    // Grid index (lower-left corner of interpolation stencil)
    int ix = static_cast<int>(fx) - spline_order/2 + 1;
    int iy = static_cast<int>(fy) - spline_order/2 + 1;
    int iz = static_cast<int>(fz) - spline_order/2 + 1;
    
    // Fractional part for B-spline evaluation
    double ux = fx - floor(fx);
    double uy = fy - floor(fy);
    double uz = fz - floor(fz);
    
    // Store grid indices (with periodic wrapping handled later)
    frac_x[i] = fx;
    frac_y[i] = fy;
    frac_z[i] = fz;
    grid_idx_x[i] = ix;
    grid_idx_y[i] = iy;
    grid_idx_z[i] = iz;
    
    // Compute B-spline coefficients
    double theta_local[PME_MAX_ORDER];
    double dtheta_local[PME_MAX_ORDER];
    
    // X direction
    compute_bspline(spline_order, ux, theta_local, dtheta_local);
    for (int k = 0; k < spline_order; ++k)
    {
      theta_x[i * spline_order + k] = theta_local[k];
      dtheta_x[i * spline_order + k] = dtheta_local[k];
    }
    
    // Y direction
    compute_bspline(spline_order, uy, theta_local, dtheta_local);
    for (int k = 0; k < spline_order; ++k)
    {
      theta_y[i * spline_order + k] = theta_local[k];
      dtheta_y[i * spline_order + k] = dtheta_local[k];
    }
    
    // Z direction
    compute_bspline(spline_order, uz, theta_local, dtheta_local);
    for (int k = 0; k < spline_order; ++k)
    {
      theta_z[i * spline_order + k] = theta_local[k];
      dtheta_z[i * spline_order + k] = dtheta_local[k];
    }
  }

  // //============================================================================
  // // KERNEL 2: Spread charges onto mesh
  // //============================================================================
  
  // __global__ void pme_spread_charges_kernel(
  //   const double* __restrict__ charges,
  //   const int* __restrict__ grid_idx_x,
  //   const int* __restrict__ grid_idx_y,
  //   const int* __restrict__ grid_idx_z,
  //   const double* __restrict__ theta_x,
  //   const double* __restrict__ theta_y,
  //   const double* __restrict__ theta_z,
  //   int n_particles,
  //   int spline_order,
  //   int nx, int ny, int nz,
  //   double* __restrict__ charge_mesh)
  // {
  //   const int i = blockIdx.x * blockDim.x + threadIdx.x;
  //   if (i >= n_particles) return;
    
  //   const double q = charges[i];
  //   const int ix0 = grid_idx_x[i];
  //   const int iy0 = grid_idx_y[i];
  //   const int iz0 = grid_idx_z[i];
    
  //   // Loop over spline support
  //   for (int kx = 0; kx < spline_order; ++kx)
  //   {
  //     // Periodic wrap for x index
  //     int ix = ix0 + kx;
  //     if (ix < 0) ix += nx;
  //     if (ix >= nx) ix -= nx;
      
  //     const double wx = theta_x[i * spline_order + kx];
      
  //     for (int ky = 0; ky < spline_order; ++ky)
  //     {
  //       // Periodic wrap for y index
  //       int iy = iy0 + ky;
  //       if (iy < 0) iy += ny;
  //       if (iy >= ny) iy -= ny;
        
  //       const double wy = theta_y[i * spline_order + ky];
  //       const double wxy = wx * wy;
        
  //       for (int kz = 0; kz < spline_order; ++kz)
  //       {
  //         // Periodic wrap for z index
  //         int iz = iz0 + kz;
  //         if (iz < 0) iz += nz;
  //         if (iz >= nz) iz -= nz;
          
  //         const double wz = theta_z[i * spline_order + kz];
  //         const double weight = q * wxy * wz;
          
  //         // Linear index into mesh
  //         const int mesh_idx = ix * ny * nz + iy * nz + iz;
          
  //         // Atomic add to handle race conditions
  //         atomicAdd(&charge_mesh[mesh_idx], weight);
  //       }
  //     }
  //   }
  // }

  // //============================================================================
  // // KERNEL 3: Apply influence function in k-space (after FFT)
  // //============================================================================
  
  // __global__ void pme_kspace_solve_kernel(
  //   cufftDoubleComplex* __restrict__ kspace_mesh,
  //   const double* __restrict__ influence_function,
  //   int nx, int ny, int nz_complex)  // nz_complex = nz/2 + 1
  // {
  //   const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  //   const int total = nx * ny * nz_complex;
  //   if (idx >= total) return;
    
  //   // Apply influence function (multiply in k-space)
  //   const double G = influence_function[idx];
  //   kspace_mesh[idx].x *= G;
  //   kspace_mesh[idx].y *= G;
  // }

  // //============================================================================
  // // KERNEL 4: Compute k-space forces (gradient of potential)
  // //============================================================================
  
  // __global__ void pme_kspace_forces_kernel(
  //   const cufftDoubleComplex* __restrict__ kspace_potential,
  //   cufftDoubleComplex* __restrict__ kspace_fx,
  //   cufftDoubleComplex* __restrict__ kspace_fy,
  //   cufftDoubleComplex* __restrict__ kspace_fz,
  //   int nx, int ny, int nz,
  //   double box_x, double box_y, double box_z)
  // {
  //   const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  //   const int nz_complex = nz / 2 + 1;
  //   const int total = nx * ny * nz_complex;
  //   if (idx >= total) return;
    
  //   // Decode 3D indices
  //   const int iz = idx % nz_complex;
  //   const int iy = (idx / nz_complex) % ny;
  //   const int ix = idx / (nz_complex * ny);
    
  //   // Compute k-vector components
  //   // For R2C transform, frequencies are:
  //   // kx: 0, 1, ..., nx/2, -nx/2+1, ..., -1  (use ix < nx/2 ? ix : ix - nx)
  //   // ky: similar
  //   // kz: 0, 1, ..., nz/2 (only positive frequencies stored)
    
  //   int mx = (ix <= nx/2) ? ix : ix - nx;
  //   int my = (iy <= ny/2) ? iy : iy - ny;
  //   int mz = iz;  // Always positive for R2C
    
  //   // k-vector in reciprocal space units
  //   double kx = 2.0 * M_PI * mx / box_x;
  //   double ky = 2.0 * M_PI * my / box_y;
  //   double kz = 2.0 * M_PI * mz / box_z;
    
  //   // Force = -i * k * potential
  //   // F(k) = -i * k * φ(k)
  //   // For complex multiplication: -i * (a + bi) = b - ai
    
  //   const double pot_re = kspace_potential[idx].x;
  //   const double pot_im = kspace_potential[idx].y;
    
  //   // Fx = -i * kx * pot = kx * pot_im - i * kx * pot_re
  //   kspace_fx[idx].x = kx * pot_im;
  //   kspace_fx[idx].y = -kx * pot_re;
    
  //   // Fy = -i * ky * pot
  //   kspace_fy[idx].x = ky * pot_im;
  //   kspace_fy[idx].y = -ky * pot_re;
    
  //   // Fz = -i * kz * pot
  //   kspace_fz[idx].x = kz * pot_im;
  //   kspace_fz[idx].y = -kz * pot_re;
  // }

  // //============================================================================
  // // KERNEL 5: Interpolate forces back to particles
  // //============================================================================
  
  // __global__ void pme_interpolate_forces_kernel(
  //   const double* __restrict__ force_mesh_x,
  //   const double* __restrict__ force_mesh_y,
  //   const double* __restrict__ force_mesh_z,
  //   const double* __restrict__ charges,
  //   const int* __restrict__ grid_idx_x,
  //   const int* __restrict__ grid_idx_y,
  //   const int* __restrict__ grid_idx_z,
  //   const double* __restrict__ theta_x,
  //   const double* __restrict__ theta_y,
  //   const double* __restrict__ theta_z,
  //   const double* __restrict__ dtheta_x,
  //   const double* __restrict__ dtheta_y,
  //   const double* __restrict__ dtheta_z,
  //   int n_particles,
  //   int spline_order,
  //   int nx, int ny, int nz,
  //   double grid_spacing_x,
  //   double grid_spacing_y,
  //   double grid_spacing_z,
  //   double* __restrict__ fx,
  //   double* __restrict__ fy,
  //   double* __restrict__ fz)
  // {
  //   const int i = blockIdx.x * blockDim.x + threadIdx.x;
  //   if (i >= n_particles) return;
    
  //   const double q = charges[i];
  //   const int ix0 = grid_idx_x[i];
  //   const int iy0 = grid_idx_y[i];
  //   const int iz0 = grid_idx_z[i];
    
  //   double force_x = 0.0;
  //   double force_y = 0.0;
  //   double force_z = 0.0;
    
  //   // Interpolate potential gradient to particle position
  //   for (int kx = 0; kx < spline_order; ++kx)
  //   {
  //     int ix = ix0 + kx;
  //     if (ix < 0) ix += nx;
  //     if (ix >= nx) ix -= nx;
      
  //     const double wx = theta_x[i * spline_order + kx];
  //     const double dwx = dtheta_x[i * spline_order + kx];
      
  //     for (int ky = 0; ky < spline_order; ++ky)
  //     {
  //       int iy = iy0 + ky;
  //       if (iy < 0) iy += ny;
  //       if (iy >= ny) iy -= ny;
        
  //       const double wy = theta_y[i * spline_order + ky];
  //       const double dwy = dtheta_y[i * spline_order + ky];
        
  //       for (int kz = 0; kz < spline_order; ++kz)
  //       {
  //         int iz = iz0 + kz;
  //         if (iz < 0) iz += nz;
  //         if (iz >= nz) iz -= nz;
          
  //         const double wz = theta_z[i * spline_order + kz];
  //         const double dwz = dtheta_z[i * spline_order + kz];
          
  //         const int mesh_idx = ix * ny * nz + iy * nz + iz;
          
  //         // Get potential at this grid point
  //         const double phi_x = force_mesh_x[mesh_idx];
  //         const double phi_y = force_mesh_y[mesh_idx];
  //         const double phi_z = force_mesh_z[mesh_idx];
          
  //         // Accumulate force using B-spline weights
  //         // F = -q * ∇φ, where gradient uses B-spline derivatives
  //         force_x += phi_x * wx * wy * wz;
  //         force_y += phi_y * wx * wy * wz;
  //         force_z += phi_z * wx * wy * wz;
  //       }
  //     }
  //   }
    
  //   // Apply charge and grid spacing scaling
  //   // The force from PME is already in correct units after this
  //   fx[i] += q * force_x;
  //   fy[i] += q * force_y;
  //   fz[i] += q * force_z;
  // }

  // //============================================================================
  // // KERNEL 6: Compute reciprocal space energy
  // //============================================================================
  
  // __global__ void pme_compute_energy_kernel(
  //   const cufftDoubleComplex* __restrict__ kspace_mesh,
  //   const double* __restrict__ influence_function,
  //   int nx, int ny, int nz_complex,
  //   double* __restrict__ energy_partial)  // Block reduction output
  // {
  //   __shared__ double sdata[PME_BLOCK_SIZE];
    
  //   const int tid = threadIdx.x;
  //   const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  //   const int total = nx * ny * nz_complex;
    
  //   // Each thread computes its contribution
  //   double local_energy = 0.0;
  //   if (idx < total)
  //   {
  //     // Energy = 0.5 * Σ_k |ρ(k)|² * G(k)
  //     // But we already applied G(k), so mesh contains φ(k) = ρ(k) * G(k)
  //     // We need the original |ρ(k)|² * G(k)
  //     // This kernel assumes it's called before kspace_solve
      
  //     const double re = kspace_mesh[idx].x;
  //     const double im = kspace_mesh[idx].y;
  //     const double rho_sq = re * re + im * im;
  //     const double G = influence_function[idx];
      
  //     // Account for complex conjugate symmetry in R2C
  //     // k=0 and k=nz/2 are not duplicated; others count twice
  //     const int iz = idx % nz_complex;
  //     double weight = (iz == 0 || iz == nz_complex - 1) ? 1.0 : 2.0;
      
  //     local_energy = 0.5 * weight * rho_sq * G;
  //   }
    
  //   sdata[tid] = local_energy;
  //   __syncthreads();
    
  //   // Block reduction
  //   for (int s = blockDim.x / 2; s > 0; s >>= 1)
  //   {
  //     if (tid < s)
  //     {
  //       sdata[tid] += sdata[tid + s];
  //     }
  //     __syncthreads();
  //   }
    
  //   // Write block result
  //   if (tid == 0)
  //   {
  //     energy_partial[blockIdx.x] = sdata[0];
  //   }
  // }

  //============================================================================
  // HOST FUNCTIONS: Launch kernels
  //============================================================================
  
  void pme_compute_bsplines(
    const double* rx, const double* ry, const double* rz,
    int n_particles,
    const PMEParameters& params,
    ParticleBSplineData& bspline_data,
    cudaStream_t stream)
  {
    if (n_particles == 0) return;
    
    bspline_data.resize(n_particles, params.spline_order);
    
    const int block_size = PME_BLOCK_SIZE;
    const int num_blocks = (n_particles + block_size - 1) / block_size;
    
    pme_compute_bsplines_kernel<<<num_blocks, block_size, 0, stream>>>(
      rx, ry, rz,
      n_particles,
      params.spline_order,
      params.grid_size.nx, params.grid_size.ny, params.grid_size.nz,
      params.inv_box_size.x, params.inv_box_size.y, params.inv_box_size.z,
      bspline_data.frac_x.data(),
      bspline_data.frac_y.data(),
      bspline_data.frac_z.data(),
      bspline_data.grid_idx_x.data(),
      bspline_data.grid_idx_y.data(),
      bspline_data.grid_idx_z.data(),
      bspline_data.theta_x.data(),
      bspline_data.theta_y.data(),
      bspline_data.theta_z.data(),
      bspline_data.dtheta_x.data(),
      bspline_data.dtheta_y.data(),
      bspline_data.dtheta_z.data()
    );
  }

  // void pme_spread_charges(
  //   const double* charges,
  //   int n_particles,
  //   const PMEParameters& params,
  //   const ParticleBSplineData& bspline_data,
  //   PMEMesh& mesh,
  //   cudaStream_t stream)
  // {
  //   if (n_particles == 0) return;
    
  //   // Clear charge mesh
  //   cudaMemsetAsync(mesh.charge_mesh.data(), 0, 
  //                   mesh.charge_mesh.size() * sizeof(double), stream);
    
  //   const int block_size = PME_BLOCK_SIZE;
  //   const int num_blocks = (n_particles + block_size - 1) / block_size;
    
  //   pme_spread_charges_kernel<<<num_blocks, block_size, 0, stream>>>(
  //     charges,
  //     bspline_data.grid_idx_x.data(),
  //     bspline_data.grid_idx_y.data(),
  //     bspline_data.grid_idx_z.data(),
  //     bspline_data.theta_x.data(),
  //     bspline_data.theta_y.data(),
  //     bspline_data.theta_z.data(),
  //     n_particles,
  //     params.spline_order,
  //     params.grid_size.nx, params.grid_size.ny, params.grid_size.nz,
  //     mesh.charge_mesh.data()
  //   );
  // }

  // void pme_solve_kspace(
  //   const PMEParameters& params,
  //   PMEMesh& mesh,
  //   cudaStream_t stream)
  // {
  //   const int nx = params.grid_size.nx;
  //   const int ny = params.grid_size.ny;
  //   const int nz = params.grid_size.nz;
  //   const int nz_complex = nz / 2 + 1;
  //   const int total_complex = nx * ny * nz_complex;
    
  //   // Set cuFFT stream
  //   cufftSetStream(mesh.fft_plan_forward, stream);
  //   cufftSetStream(mesh.fft_plan_inverse, stream);
    
  //   // Forward FFT: real charge mesh -> complex k-space
  //   cufftExecD2Z(mesh.fft_plan_forward,
  //                mesh.charge_mesh.data(),
  //                mesh.kspace_mesh.data());
    
  //   // Apply influence function
  //   const int block_size = PME_BLOCK_SIZE;
  //   const int num_blocks = (total_complex + block_size - 1) / block_size;
    
  //   pme_kspace_solve_kernel<<<num_blocks, block_size, 0, stream>>>(
  //     mesh.kspace_mesh.data(),
  //     params.influence_function.data(),
  //     nx, ny, nz_complex
  //   );
  // }

  // void pme_compute_forces_kspace(
  //   const PMEParameters& params,
  //   PMEMesh& mesh,
  //   cudaStream_t stream)
  // {
  //   const int nx = params.grid_size.nx;
  //   const int ny = params.grid_size.ny;
  //   const int nz = params.grid_size.nz;
  //   const int nz_complex = nz / 2 + 1;
  //   const int total_complex = nx * ny * nz_complex;
    
  //   // Allocate temporary k-space force arrays
  //   static onika::memory::CudaMMVector<cufftDoubleComplex> kspace_fx, kspace_fy, kspace_fz;
  //   if (kspace_fx.size() < total_complex)
  //   {
  //     kspace_fx.resize(total_complex);
  //     kspace_fy.resize(total_complex);
  //     kspace_fz.resize(total_complex);
  //   }
    
  //   // Compute k-space forces (gradient)
  //   const int block_size = PME_BLOCK_SIZE;
  //   const int num_blocks = (total_complex + block_size - 1) / block_size;
    
  //   pme_kspace_forces_kernel<<<num_blocks, block_size, 0, stream>>>(
  //     mesh.kspace_mesh.data(),
  //     kspace_fx.data(),
  //     kspace_fy.data(),
  //     kspace_fz.data(),
  //     nx, ny, nz,
  //     params.box_size.x, params.box_size.y, params.box_size.z
  //   );
    
  //   // Inverse FFT for each force component
  //   cufftSetStream(mesh.fft_plan_inverse, stream);
    
  //   cufftExecZ2D(mesh.fft_plan_inverse, kspace_fx.data(), mesh.force_mesh_x.data());
  //   cufftExecZ2D(mesh.fft_plan_inverse, kspace_fy.data(), mesh.force_mesh_y.data());
  //   cufftExecZ2D(mesh.fft_plan_inverse, kspace_fz.data(), mesh.force_mesh_z.data());
    
  //   // Note: cuFFT does not normalize inverse transform
  //   // Normalization by 1/(nx*ny*nz) should be applied in interpolation
  // }

  // void pme_interpolate_forces(
  //   const double* charges,
  //   int n_particles,
  //   const PMEParameters& params,
  //   const ParticleBSplineData& bspline_data,
  //   const PMEMesh& mesh,
  //   double* fx, double* fy, double* fz,
  //   cudaStream_t stream)
  // {
  //   if (n_particles == 0) return;
    
  //   const int block_size = PME_BLOCK_SIZE;
  //   const int num_blocks = (n_particles + block_size - 1) / block_size;
    
  //   pme_interpolate_forces_kernel<<<num_blocks, block_size, 0, stream>>>(
  //     mesh.force_mesh_x.data(),
  //     mesh.force_mesh_y.data(),
  //     mesh.force_mesh_z.data(),
  //     charges,
  //     bspline_data.grid_idx_x.data(),
  //     bspline_data.grid_idx_y.data(),
  //     bspline_data.grid_idx_z.data(),
  //     bspline_data.theta_x.data(),
  //     bspline_data.theta_y.data(),
  //     bspline_data.theta_z.data(),
  //     bspline_data.dtheta_x.data(),
  //     bspline_data.dtheta_y.data(),
  //     bspline_data.dtheta_z.data(),
  //     n_particles,
  //     params.spline_order,
  //     params.grid_size.nx, params.grid_size.ny, params.grid_size.nz,
  //     params.grid_spacing.x, params.grid_spacing.y, params.grid_spacing.z,
  //     fx, fy, fz
  //   );
  // }

}
