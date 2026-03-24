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

// Implementation of the 2-body kernel descriptor / potential from:
//   Dézaphie et al., Computational Materials Science 246 (2025) 113459
//   "Designing hybrid descriptors for improved machine learning models
//    in atomistic materials science simulations"
//
// The local energy contribution of a pair at distance r is given by Eq. (19):
//
//   e(r) = Delta * sum_{k=1}^{K2b} w_k * exp( -(r - s_k)^2 / (2*sigma^2) )
//
// where s_k = r_min + (k-1) * (r_cut - r_min) / (K2b - 1)  is a regular radial grid,
// sigma controls the Gaussian width, Delta is an energy scale, and w_k are the
// learned ML weights (applied externally by the framework).
//
// This function evaluates the *weighted* kernel sum (full Eq. 19):
//   e  = Delta * sum_k  w_k * exp( -(r - s_k)^2 / (2*sigma^2) )
//   de = d(e)/d(r)
//      = -Delta * sum_k  w_k * (r - s_k)/sigma^2 * exp( -(r - s_k)^2 / (2*sigma^2) )
//
// The weight vector w (size K2b) is stored in K2bPotentialParameters and loaded
// from the YAML config alongside the other hyperparameters.

#pragma once

#include <cassert>
#include <cmath>
#include <vector>
#include <yaml-cpp/yaml.h>

#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{

struct K2bPotentialParameters
{
  // Number of radial grid points (K_2b in the paper). Default: 40 (paper value).
  int    n_rbf    = 40;

  // Minimum interatomic distance for the radial grid start [Angstrom].
  // Points below r_min receive negligible kernel weight due to the Gaussian decay.
  double r_min    = 0.0;

  // Cutoff radius: upper bound of the regular radial grid [Angstrom].
  double r_cut    = 6.0;

  // Gaussian width sigma [Angstrom]. Controls smoothness of each basis function.
  // Paper value: sigma = 0.2 Å.
  double sigma    = 0.2;

  // Energy scale Delta [eV]. Calibrates the typical energy range of the kernel.
  // Paper value: Delta = 2 eV.
  double delta    = 2.0;

  // Learned ML weight vector w_2b (size = n_rbf).
  // w[k] is the coefficient for the k-th Gaussian basis function (Eq. 19).
  // Must contain exactly n_rbf entries before calling k2b_potential_compute_force.
  std::vector<double> w = {};
};

} // namespace exaStamp


// ---------------------------------------------------------------------------
// YAML decode: read K2bPotentialParameters from a config node.
// Expected YAML keys (all optional except `w`, fall back to defaults above):
//
//   k2b:
//     n_rbf: 40
//     r_min: 0.0
//     r_cut: 6.0
//     sigma: 0.2
//     delta: 2.0
//     w: [w0, w1, ..., w_{n_rbf-1}]   # must match n_rbf entries
// ---------------------------------------------------------------------------
namespace YAML
{
  template<> struct convert<exaStamp::K2bPotentialParameters>
  {
    static bool decode(const Node& node, exaStamp::K2bPotentialParameters& v)
    {
      if( node["n_rbf"] ) v.n_rbf  = node["n_rbf"].as<int>();
      if( node["r_min"] ) v.r_min  = node["r_min"].as<double>();
      if( node["r_cut"] ) v.r_cut  = node["r_cut"].as<double>();
      if( node["sigma"] ) v.sigma  = node["sigma"].as<double>();
      if( node["delta"] ) v.delta  = node["delta"].as<double>();
      if( node["w"] )
      {
        v.w = node["w"].as< std::vector<double> >();
        // Consistency check: w must have exactly n_rbf entries.
        if( static_cast<int>(v.w.size()) != v.n_rbf )
          return false;
      }
      return true;
    }
  };
} // namespace YAML


namespace exaStamp
{

// ---------------------------------------------------------------------------
// k2b_potential_compute_force
//
// Evaluates the full weighted 2-body kernel energy and its radial derivative
// (Eq. 19 of Dézaphie et al. 2025):
//
//   e  = Delta * sum_{k=0}^{K-1}  w[k] * exp( -(r-s_k)^2 / (2*sigma^2) )
//   de = d(e)/d(r)
//      = -Delta * sum_{k=0}^{K-1}  w[k] * (r-s_k)/sigma^2
//                                        * exp( -(r-s_k)^2 / (2*sigma^2) )
//
// Parameters:
//   params      : K2bPotentialParameters  (n_rbf, r_min, r_cut, sigma, delta, w)
//   pair_params : PairPotentialMinimalParameters  (not used; purely radial)
//   r           : interatomic distance [Angstrom], must be > 0
//   e  (out)    : pair energy contribution [eV]
//   de (out)    : d(e)/d(r) [eV/Angstrom]
//                 The framework derives the force as  F = -de/r * r_vec
// ---------------------------------------------------------------------------
ONIKA_HOST_DEVICE_FUNC inline void k2b_potential_compute_force(
    const K2bPotentialParameters&        params,
    const PairPotentialMinimalParameters& /*pair_params*/,
    double r,
    double& e,
    double& de)
{
  assert( r > 0. );
  assert( static_cast<int>(params.w.size()) == params.n_rbf );

  e  = 0.0;
  de = 0.0;

  const int     K       = params.n_rbf;
  const double  r_min   = params.r_min;
  const double  r_cut   = params.r_cut;
  const double  sigma   = params.sigma;
  const double  delta   = params.delta;
  const double* w       = params.w.data();

  // Pre-compute shared constants.
  const double inv_2sig2 = 0.5 / (sigma * sigma);  // 1 / (2 sigma^2)
  const double inv_sig2  = 1.0 / (sigma * sigma);  // 1 / sigma^2

  // Grid spacing: K points equally spaced on [r_min, r_cut].
  // s_k = r_min + k * h,  k = 0 ... K-1
  const double h = (K > 1) ? (r_cut - r_min) / static_cast<double>(K - 1) : 0.0;

  // Accumulate weighted Gaussian basis functions.
  for (int k = 0; k < K; ++k)
  {
    const double s_k = r_min + k * h;           // grid centre for basis k
    const double dr  = r - s_k;                 // displacement from centre
    const double arg = dr * dr * inv_2sig2;     // (r - s_k)^2 / (2 sigma^2)

    // Skip negligible contributions (exp(-20) ≈ 2e-9).
    if (arg > 20.0) continue;

    const double g = std::exp(-arg);             // Gaussian kernel value

    e  += w[k] * g;
    de -= w[k] * dr * inv_sig2 * g;             // d/dr [w_k * exp(...)]
  }

  // Apply the energy scale Delta (Eq. 19).
  e  *= delta;
  de *= delta;
}

} // namespace exaStamp


#define USTAMP_POTENTIAL_NAME     k2b
#define USTAMP_POTENTIAL_PARAMS   K2bPotentialParameters
#define USTAMP_POTENTIAL_COMPUTE  k2b_potential_compute_force

#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0) // PairPotentialParameters not used in computation

#define USTAMP_POTENTIAL_ENABLE_CUDA 1
