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

// Gaussian mixture model (full covariances) on strain functional descriptor columns.
// Same conventions as scikit-learn's GaussianMixture, so that models trained with
// scripts/python/sfd_gmm_train.py and with sfd_gmm_fit are interchangeable:
//   - features are optionally standardized: x = (row[col] - center) / scale
//   - log p_k(x) = log w_k + log det(P_k) - D/2 log(2 pi) - 1/2 |(x - mu_k)^T P_k|^2
//     with P_k the upper-triangular "precision Cholesky" factor, Sigma_k^-1 = P_k P_k^T

#include <onika/cuda/cuda.h>
#include <onika/memory/allocator.h>
#include <exaStamp/analysis_particle/strain_functional/sfd_math.h>
#include <string>
#include <vector>
#include <cmath>

namespace exaStamp
{
  namespace sfd
  {
    static constexpr int GMM_MAX_FEATURES = N_INVARIANTS;
    static constexpr int GMM_MAX_COMPONENTS = 16;

    // host model, all arrays flattened (component-major)
    struct SFDGMMModel
    {
      std::vector<std::string> features;            // D column names of sfd_invariants
      onika::memory::CudaMMVector<int> columns;     // D column indices
      onika::memory::CudaMMVector<double> center;   // D (0 if not standardized)
      onika::memory::CudaMMVector<double> scale;    // D (1 if not standardized)
      onika::memory::CudaMMVector<double> weights;  // K
      onika::memory::CudaMMVector<double> means;    // K*D
      onika::memory::CudaMMVector<double> covariances; // K*D*D
      // derived
      onika::memory::CudaMMVector<double> prec_chol;   // K*D*D, upper triangular, row-major
      onika::memory::CudaMMVector<double> log_norm;    // K: log w_k + log det P_k - D/2 log(2 pi)

      inline int n_features() const { return static_cast<int>( features.size() ); }
      inline int n_components() const { return static_cast<int>( weights.size() ); }
    };

    // device view
    struct GMMView
    {
      int K = 0;
      int D = 0;
      const int * __restrict__ columns = nullptr;
      const double * __restrict__ center = nullptr;
      const double * __restrict__ scale = nullptr;
      const double * __restrict__ means = nullptr;
      const double * __restrict__ prec_chol = nullptr;
      const double * __restrict__ log_norm = nullptr;
    };

    inline GMMView gmm_view( const SFDGMMModel& m )
    {
      return GMMView{ m.n_components(), m.n_features(), m.columns.data(), m.center.data(), m.scale.data()
                    , m.means.data(), m.prec_chol.data(), m.log_norm.data() };
    }

    // standardized features of one particle (row = its N_INVARIANTS values)
    ONIKA_HOST_DEVICE_FUNC inline void gmm_features( const GMMView& g, const double* row, double* x )
    {
      for(int d=0;d<g.D;d++) x[d] = ( row[ g.columns[d] ] - g.center[d] ) / g.scale[d];
    }

    // per-component weighted log densities lp[k] = log(w_k N(x | mu_k, Sigma_k)), returns log sum_k exp(lp[k])
    ONIKA_HOST_DEVICE_FUNC inline double gmm_log_prob( const GMMView& g, const double* x, double* lp )
    {
      double mx = -1.0e300;
      for(int k=0;k<g.K;k++)
      {
        const double* mu = g.means + k*g.D;
        const double* P = g.prec_chol + k*g.D*g.D;
        double q = 0.0;
        for(int j=0;j<g.D;j++)
        {
          double y = 0.0;
          for(int i=0;i<=j;i++) y += ( x[i] - mu[i] ) * P[ i*g.D + j ];
          q += y*y;
        }
        lp[k] = g.log_norm[k] - 0.5 * q;
        if( lp[k] > mx ) mx = lp[k];
      }
      double s = 0.0;
      for(int k=0;k<g.K;k++) s += exp( lp[k] - mx );
      return mx + log( s );
    }

    // most probable component and its posterior probability
    ONIKA_HOST_DEVICE_FUNC inline int gmm_classify( const GMMView& g, const double* row, double& posterior )
    {
      double x[GMM_MAX_FEATURES];
      double lp[GMM_MAX_COMPONENTS];
      gmm_features( g, row, x );
      const double lse = gmm_log_prob( g, x, lp );
      int best = 0;
      for(int k=1;k<g.K;k++) if( lp[k] > lp[best] ) best = k;
      posterior = exp( lp[best] - lse );
      return best;
    }
  }
}
