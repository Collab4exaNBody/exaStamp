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

#include "pme.h"
#include <cuda_runtime.h>

namespace exaStamp
{
  //============================================================================
  // PME Kernel Host-Side Launch Functions
  //============================================================================
  
  /**
   * @brief Compute B-spline coefficients for all particles
   * 
   * @param rx, ry, rz    Particle positions (device pointers)
   * @param n_particles   Number of particles
   * @param params        PME parameters
   * @param bspline_data  Output: B-spline coefficients and grid indices
   * @param stream        CUDA stream
   */
  void pme_compute_bsplines(
    const double* rx, const double* ry, const double* rz,
    int n_particles,
    const PMEParameters& params,
    ParticleBSplineData& bspline_data,
    cudaStream_t stream = 0);

  /**
   * @brief Spread particle charges onto the mesh using B-spline interpolation
   * 
   * @param charges       Particle charges (device pointer)
   * @param n_particles   Number of particles
   * @param params        PME parameters
   * @param bspline_data  Precomputed B-spline coefficients
   * @param mesh          Output: Charge mesh (will be cleared first)
   * @param stream        CUDA stream
   */
  void pme_spread_charges(
    const double* charges,
    int n_particles,
    const PMEParameters& params,
    const ParticleBSplineData& bspline_data,
    PMEMesh& mesh,
    cudaStream_t stream = 0);

  /**
   * @brief Solve in k-space: FFT + apply influence function
   * 
   * @param params        PME parameters (contains influence function)
   * @param mesh          Input/Output: charge mesh -> potential mesh
   * @param stream        CUDA stream
   */
  void pme_solve_kspace(
    const PMEParameters& params,
    PMEMesh& mesh,
    cudaStream_t stream = 0);

  /**
   * @brief Compute forces in k-space and inverse FFT
   * 
   * @param params        PME parameters
   * @param mesh          Input: k-space potential, Output: real-space forces
   * @param stream        CUDA stream
   */
  void pme_compute_forces_kspace(
    const PMEParameters& params,
    PMEMesh& mesh,
    cudaStream_t stream = 0);

  /**
   * @brief Interpolate forces from mesh back to particles
   * 
   * @param charges       Particle charges (device pointer)
   * @param n_particles   Number of particles
   * @param params        PME parameters
   * @param bspline_data  Precomputed B-spline coefficients
   * @param mesh          Force mesh (from pme_compute_forces_kspace)
   * @param fx, fy, fz    Output: Particle forces (accumulated)
   * @param stream        CUDA stream
   */
  void pme_interpolate_forces(
    const double* charges,
    int n_particles,
    const PMEParameters& params,
    const ParticleBSplineData& bspline_data,
    const PMEMesh& mesh,
    double* fx, double* fy, double* fz,
    cudaStream_t stream = 0);

  /**
   * @brief Compute reciprocal space energy
   * 
   * @param params        PME parameters
   * @param mesh          k-space mesh (before applying influence function)
   * @param stream        CUDA stream
   * @return              Reciprocal space energy
   */
  double pme_compute_energy(
    const PMEParameters& params,
    const PMEMesh& mesh,
    cudaStream_t stream = 0);

}
