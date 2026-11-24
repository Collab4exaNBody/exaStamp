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

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/memory/allocator.h>

#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exaStamp/analysis_particle/strain_functional/sfd_gmm.h>

#include <yaml-cpp/yaml.h>
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

// Gaussian mixture classification of particles from their strain functional descriptors,
// following the unsupervised analysis of Kober et al. (arXiv:2402.04191, section 3.6):
// sfd_gmm_fit fits the mixture on the current configuration (EM, MPI parallel, host),
// sfd_gmm_classify assigns each particle to its most probable component (GPU compatible).
namespace exaStamp
{
  using namespace exanb;
  using sfd::SFDGMMModel;

  namespace sfd_gmm_detail
  {
    static constexpr double LOG_2PI = 1.8378770664093453;

    // feature names -> sfd_invariants column indices
    inline std::vector<int> resolve_columns( const std::vector<std::string>& features )
    {
      const char* const* names = sfd::sfd_invariant_names();
      std::vector<int> cols;
      for(const auto& f : features)
      {
        const auto it = std::find_if( names, names + sfd::N_INVARIANTS, [&f](const char* s){ return f == s; } );
        if( it == names + sfd::N_INVARIANTS ) fatal_error() << "sfd_gmm: unknown feature '" << f << "'" << std::endl;
        cols.push_back( static_cast<int>( it - names ) );
      }
      return cols;
    }

    inline std::vector<std::string> default_features()
    {
      // the 23 rotational invariants of orders 1 to 4 (no total weight, no orientation factors)
      const char* const* names = sfd::sfd_invariant_names();
      return std::vector<std::string>( names + 1 , names + 24 );
    }

    // Sigma = L L^T, P = (L^-1)^T (upper triangular), log_norm = log w + log det P - D/2 log(2 pi)
    inline void update_precisions( SFDGMMModel& m )
    {
      const int K = m.n_components(), D = m.n_features();
      m.prec_chol.assign( size_t(K)*D*D , 0.0 );
      m.log_norm.assign( K , 0.0 );
      std::vector<double> L( D*D ), Linv( D*D );
      for(int k=0;k<K;k++)
      {
        const double* S = m.covariances.data() + size_t(k)*D*D;
        std::fill( L.begin(), L.end(), 0.0 );
        for(int j=0;j<D;j++)
        {
          double s = S[j*D+j];
          for(int p=0;p<j;p++) s -= L[j*D+p]*L[j*D+p];
          if( !( s > 0.0 ) )
          {
            fatal_error() << "sfd_gmm: covariance of component " << k << " is not positive definite;"
                          << " increase reg_covar or use fewer components/features" << std::endl;
          }
          L[j*D+j] = std::sqrt( s );
          for(int i=j+1;i<D;i++)
          {
            double t = S[i*D+j];
            for(int p=0;p<j;p++) t -= L[i*D+p]*L[j*D+p];
            L[i*D+j] = t / L[j*D+j];
          }
        }
        // Linv: lower triangular inverse by forward substitution
        std::fill( Linv.begin(), Linv.end(), 0.0 );
        for(int c=0;c<D;c++)
        {
          for(int i=c;i<D;i++)
          {
            double t = ( i == c ) ? 1.0 : 0.0;
            for(int p=c;p<i;p++) t -= L[i*D+p] * Linv[p*D+c];
            Linv[i*D+c] = t / L[i*D+i];
          }
        }
        double log_det = 0.0;
        double* P = m.prec_chol.data() + size_t(k)*D*D;
        for(int i=0;i<D;i++)
        {
          log_det -= std::log( L[i*D+i] );
          for(int j=0;j<D;j++) P[i*D+j] = Linv[j*D+i];
        }
        m.log_norm[k] = std::log( m.weights[k] ) + log_det - 0.5 * D * LOG_2PI;
      }
    }

    // ---------------------------------------------------------------------------------
    // YAML model I/O (same layout as scripts/python/sfd_gmm_train.py)
    // ---------------------------------------------------------------------------------
    inline void write_model( const SFDGMMModel& m, const std::string& file_name )
    {
      const int K = m.n_components(), D = m.n_features();
      YAML::Emitter out;
      out.SetDoublePrecision( 17 );
      out << YAML::BeginMap;
      out << YAML::Key << "features" << YAML::Value << YAML::Flow << m.features;
      out << YAML::Key << "center" << YAML::Value << YAML::Flow << std::vector<double>( m.center.begin(), m.center.end() );
      out << YAML::Key << "scale" << YAML::Value << YAML::Flow << std::vector<double>( m.scale.begin(), m.scale.end() );
      out << YAML::Key << "weights" << YAML::Value << YAML::Flow << std::vector<double>( m.weights.begin(), m.weights.end() );
      out << YAML::Key << "means" << YAML::Value << YAML::BeginSeq;
      for(int k=0;k<K;k++) out << YAML::Flow << std::vector<double>( m.means.begin() + size_t(k)*D , m.means.begin() + size_t(k+1)*D );
      out << YAML::EndSeq;
      out << YAML::Key << "covariances" << YAML::Value << YAML::BeginSeq;
      for(int k=0;k<K;k++)
      {
        out << YAML::BeginSeq;
        for(int i=0;i<D;i++) out << YAML::Flow << std::vector<double>( m.covariances.begin() + (size_t(k)*D+i)*D , m.covariances.begin() + (size_t(k)*D+i+1)*D );
        out << YAML::EndSeq;
      }
      out << YAML::EndSeq;
      out << YAML::EndMap;
      std::ofstream f( file_name );
      f << "# Gaussian mixture model on strain functional descriptors (sfd_gmm_fit / sfd_gmm_train.py)\n" << out.c_str() << "\n";
    }

    inline SFDGMMModel read_model( const std::string& file_name )
    {
      const YAML::Node n = YAML::LoadFile( file_name );
      SFDGMMModel m;
      m.features = n["features"].as< std::vector<std::string> >();
      const int D = m.n_features();
      const auto w = n["weights"].as< std::vector<double> >();
      const int K = static_cast<int>( w.size() );
      if( D < 1 || D > sfd::GMM_MAX_FEATURES ) fatal_error() << "sfd_gmm: " << file_name << ": 1 to " << sfd::GMM_MAX_FEATURES << " features expected" << std::endl;
      if( K < 1 || K > sfd::GMM_MAX_COMPONENTS ) fatal_error() << "sfd_gmm: " << file_name << ": 1 to " << sfd::GMM_MAX_COMPONENTS << " components expected" << std::endl;
      const auto cols = resolve_columns( m.features );
      m.columns.assign( cols.begin(), cols.end() );
      m.center.assign( D, 0.0 );
      m.scale.assign( D, 1.0 );
      if( n["center"] ) { const auto v = n["center"].as< std::vector<double> >(); if( int(v.size()) != D ) fatal_error() << "sfd_gmm: bad center size" << std::endl; m.center.assign( v.begin(), v.end() ); }
      if( n["scale"] ) { const auto v = n["scale"].as< std::vector<double> >(); if( int(v.size()) != D ) fatal_error() << "sfd_gmm: bad scale size" << std::endl; m.scale.assign( v.begin(), v.end() ); }
      m.weights.assign( w.begin(), w.end() );
      const auto mu = n["means"].as< std::vector< std::vector<double> > >();
      const auto cov = n["covariances"].as< std::vector< std::vector< std::vector<double> > > >();
      if( int(mu.size()) != K || int(cov.size()) != K ) fatal_error() << "sfd_gmm: " << file_name << ": means/covariances do not match weights" << std::endl;
      for(int k=0;k<K;k++)
      {
        if( int(mu[k].size()) != D || int(cov[k].size()) != D ) fatal_error() << "sfd_gmm: " << file_name << ": bad means/covariances size" << std::endl;
        m.means.insert( m.means.end(), mu[k].begin(), mu[k].end() );
        for(int i=0;i<D;i++)
        {
          if( int(cov[k][i].size()) != D ) fatal_error() << "sfd_gmm: " << file_name << ": bad covariance size" << std::endl;
          m.covariances.insert( m.covariances.end(), cov[k][i].begin(), cov[k][i].end() );
        }
      }
      update_precisions( m );
      return m;
    }

    // ---------------------------------------------------------------------------------
    // fit helpers
    // ---------------------------------------------------------------------------------
    inline uint64_t splitmix64( uint64_t x )
    {
      x += 0x9E3779B97F4A7C15ULL;
      x = ( x ^ ( x >> 30 ) ) * 0xBF58476D1CE4E5B9ULL;
      x = ( x ^ ( x >> 27 ) ) * 0x94D049BB133111EBULL;
      return x ^ ( x >> 31 );
    }

    // uniform in [0,1) from the raw 64 bit engine output (same on every platform)
    inline double uniform01( std::mt19937_64& rng ) { return double( rng() >> 11 ) * 0x1.0p-53; }

    // data owned by this rank
    struct LocalData
    {
      int D = 0;
      size_t n = 0;
      std::vector<double> x;    // n*D standardized features
      std::vector<uint64_t> id; // n particle ids
    };

    // k-means++ seeding + Lloyd iterations on a sample (rank 0), returns K*D centers
    inline std::vector<double> kmeans( const std::vector<double>& s, int D, int K, uint64_t seed )
    {
      const size_t n = s.size() / D;
      std::mt19937_64 rng( seed );
      std::vector<double> c( size_t(K)*D );
      auto dist2 = [&](const double* a, const double* b){ double r=0; for(int d=0;d<D;d++){ const double t=a[d]-b[d]; r+=t*t; } return r; };
      const size_t first = static_cast<size_t>( rng() % n );
      std::copy( s.begin() + first*D, s.begin() + (first+1)*D, c.begin() );
      std::vector<double> dmin( n );
      for(size_t i=0;i<n;i++) dmin[i] = dist2( &s[i*D], &c[0] );
      for(int k=1;k<K;k++)
      {
        const double tot = std::accumulate( dmin.begin(), dmin.end(), 0.0 );
        double r = uniform01( rng ) * tot;
        size_t pick = n-1;
        for(size_t i=0;i<n;i++) { r -= dmin[i]; if( r < 0.0 ) { pick = i; break; } }
        std::copy( s.begin() + pick*D, s.begin() + (pick+1)*D, c.begin() + size_t(k)*D );
        for(size_t i=0;i<n;i++) dmin[i] = std::min( dmin[i], dist2( &s[i*D], &c[size_t(k)*D] ) );
      }
      std::vector<int> label( n, -1 );
      for(int it=0;it<300;it++)
      {
        bool changed = false;
        for(size_t i=0;i<n;i++)
        {
          int best = 0; double bd = dist2( &s[i*D], &c[0] );
          for(int k=1;k<K;k++) { const double dd = dist2( &s[i*D], &c[size_t(k)*D] ); if( dd < bd ) { bd = dd; best = k; } }
          if( best != label[i] ) { label[i] = best; changed = true; }
        }
        if( !changed ) break;
        std::vector<double> sum( size_t(K)*D, 0.0 ); std::vector<size_t> cnt( K, 0 );
        for(size_t i=0;i<n;i++) { ++cnt[label[i]]; for(int d=0;d<D;d++) sum[size_t(label[i])*D+d] += s[i*D+d]; }
        for(int k=0;k<K;k++) if( cnt[k] > 0 ) for(int d=0;d<D;d++) c[size_t(k)*D+d] = sum[size_t(k)*D+d] / double(cnt[k]);
      }
      return c;
    }

    // M-step from responsibilities resp (n*K) : weights, means, covariances (+ reg_covar on the diagonal)
    inline void m_step( const LocalData& data, const std::vector<double>& resp, int K, double reg_covar, MPI_Comm comm, SFDGMMModel& m )
    {
      const int D = data.D;
      const size_t n = data.n;
      const double eps10 = 10.0 * std::numeric_limits<double>::epsilon();
      // pass 1: nk, sum r x
      std::vector<double> s1( size_t(K)*(1+D), 0.0 );
#     pragma omp parallel
      {
        std::vector<double> loc( s1.size(), 0.0 );
#       pragma omp for schedule(static)
        for(size_t i=0;i<n;i++)
        {
          for(int k=0;k<K;k++)
          {
            const double r = resp[i*K+k];
            double* a = &loc[size_t(k)*(1+D)];
            a[0] += r;
            for(int d=0;d<D;d++) a[1+d] += r * data.x[i*D+d];
          }
        }
#       pragma omp critical
        for(size_t j=0;j<s1.size();j++) s1[j] += loc[j];
      }
      MPI_Allreduce( MPI_IN_PLACE, s1.data(), int(s1.size()), MPI_DOUBLE, MPI_SUM, comm );
      std::vector<double> nk( K );
      m.means.assign( size_t(K)*D, 0.0 );
      for(int k=0;k<K;k++)
      {
        nk[k] = s1[size_t(k)*(1+D)] + eps10;
        for(int d=0;d<D;d++) m.means[size_t(k)*D+d] = s1[size_t(k)*(1+D)+1+d] / nk[k];
      }
      // pass 2: covariances about the new means
      std::vector<double> s2( size_t(K)*D*D, 0.0 );
#     pragma omp parallel
      {
        std::vector<double> loc( s2.size(), 0.0 );
        std::vector<double> diff( D );
#       pragma omp for schedule(static)
        for(size_t i=0;i<n;i++)
        {
          for(int k=0;k<K;k++)
          {
            const double r = resp[i*K+k];
            if( r == 0.0 ) continue;
            for(int d=0;d<D;d++) diff[d] = data.x[i*D+d] - m.means[size_t(k)*D+d];
            double* a = &loc[size_t(k)*D*D];
            for(int p=0;p<D;p++) { const double rp = r*diff[p]; for(int q=p;q<D;q++) a[p*D+q] += rp*diff[q]; }
          }
        }
#       pragma omp critical
        for(size_t j=0;j<s2.size();j++) s2[j] += loc[j];
      }
      MPI_Allreduce( MPI_IN_PLACE, s2.data(), int(s2.size()), MPI_DOUBLE, MPI_SUM, comm );
      m.covariances.assign( size_t(K)*D*D, 0.0 );
      double ntot = 0.0;
      for(int k=0;k<K;k++) ntot += nk[k];
      m.weights.assign( K, 0.0 );
      for(int k=0;k<K;k++)
      {
        double* C = m.covariances.data() + size_t(k)*D*D;
        for(int p=0;p<D;p++) for(int q=p;q<D;q++) { C[p*D+q] = C[q*D+p] = s2[size_t(k)*D*D+p*D+q] / nk[k]; }
        for(int p=0;p<D;p++) C[p*D+p] += reg_covar;
        m.weights[k] = nk[k] / ntot;
      }
      update_precisions( m );
    }

    // E-step: responsibilities and mean log-likelihood (global)
    inline double e_step( const LocalData& data, const SFDGMMModel& m, std::vector<double>* resp, std::vector<int>* labels, double ntot, MPI_Comm comm )
    {
      const sfd::GMMView g = sfd::gmm_view( m );
      const int K = g.K, D = g.D;
      double sum_lse = 0.0;
#     pragma omp parallel for schedule(static) reduction(+:sum_lse)
      for(size_t i=0;i<data.n;i++)
      {
        double lp[sfd::GMM_MAX_COMPONENTS];
        const double lse = sfd::gmm_log_prob( g, &data.x[i*D], lp );
        sum_lse += lse;
        if( resp != nullptr ) for(int k=0;k<K;k++) (*resp)[i*K+k] = std::exp( lp[k] - lse );
        if( labels != nullptr ) { int b=0; for(int k=1;k<K;k++) if( lp[k] > lp[b] ) b=k; (*labels)[i] = b; }
      }
      MPI_Allreduce( MPI_IN_PLACE, &sum_lse, 1, MPI_DOUBLE, MPI_SUM, comm );
      return sum_lse / ntot;
    }

    struct FitScores
    {
      int K = 0;
      int n_iter = 0;
      bool converged = false;
      double log_likelihood = 0.0; // mean per particle
      double bic = 0.0;
      double chi = 0.0;            // Calinski-Harabasz
      double dbi = 0.0;            // Davies-Bouldin
    };

    // Calinski-Harabasz and Davies-Bouldin scores of the hard labels (sklearn definitions)
    inline void cluster_scores( const LocalData& data, const std::vector<int>& labels, int K, double ntot, MPI_Comm comm, FitScores& sc )
    {
      const int D = data.D;
      std::vector<double> s( size_t(K)*(1+D) + D, 0.0 ); // per cluster count + sum, then global sum
      for(size_t i=0;i<data.n;i++)
      {
        double* a = &s[size_t(labels[i])*(1+D)];
        a[0] += 1.0;
        for(int d=0;d<D;d++) { a[1+d] += data.x[i*D+d]; s[size_t(K)*(1+D)+d] += data.x[i*D+d]; }
      }
      MPI_Allreduce( MPI_IN_PLACE, s.data(), int(s.size()), MPI_DOUBLE, MPI_SUM, comm );
      std::vector<double> cent( size_t(K)*D, 0.0 ), mean( D );
      std::vector<double> cnt( K );
      int n_nonempty = 0;
      for(int d=0;d<D;d++) mean[d] = s[size_t(K)*(1+D)+d] / ntot;
      for(int k=0;k<K;k++)
      {
        cnt[k] = s[size_t(k)*(1+D)];
        if( cnt[k] > 0 ) { ++ n_nonempty; for(int d=0;d<D;d++) cent[size_t(k)*D+d] = s[size_t(k)*(1+D)+1+d] / cnt[k]; }
      }
      // within dispersion and mean distance to centroid per cluster
      std::vector<double> w( K+1, 0.0 );
      for(size_t i=0;i<data.n;i++)
      {
        const int k = labels[i];
        double d2 = 0.0;
        for(int d=0;d<D;d++) { const double t = data.x[i*D+d] - cent[size_t(k)*D+d]; d2 += t*t; }
        w[K] += d2;
        w[k] += std::sqrt( d2 );
      }
      MPI_Allreduce( MPI_IN_PLACE, w.data(), K+1, MPI_DOUBLE, MPI_SUM, comm );
      double extra = 0.0;
      for(int k=0;k<K;k++) if( cnt[k] > 0 ) { double d2=0; for(int d=0;d<D;d++){ const double t=cent[size_t(k)*D+d]-mean[d]; d2+=t*t; } extra += cnt[k]*d2; }
      const double intra = w[K];
      sc.chi = ( n_nonempty < 2 ) ? 0.0 : ( ( intra == 0.0 ) ? 1.0 : extra * ( ntot - n_nonempty ) / ( intra * ( n_nonempty - 1 ) ) );
      // Davies-Bouldin over the non-empty clusters
      double db = 0.0;
      for(int a=0;a<K;a++)
      {
        if( cnt[a] <= 0 ) continue;
        const double sa = w[a] / cnt[a];
        double worst = 0.0;
        for(int b=0;b<K;b++)
        {
          if( b == a || cnt[b] <= 0 ) continue;
          double d2=0; for(int d=0;d<D;d++){ const double t=cent[size_t(a)*D+d]-cent[size_t(b)*D+d]; d2+=t*t; }
          const double sep = std::sqrt( d2 );
          const double r = ( sep > 0.0 ) ? ( sa + w[b]/cnt[b] ) / sep : 0.0;
          worst = std::max( worst, r );
        }
        db += worst;
      }
      sc.dbi = ( n_nonempty < 2 ) ? 0.0 : db / n_nonempty;
    }

    // components sorted by decreasing weight (class 0 = most populated)
    inline void sort_components( SFDGMMModel& m )
    {
      const int K = m.n_components(), D = m.n_features();
      std::vector<int> order( K );
      std::iota( order.begin(), order.end(), 0 );
      std::stable_sort( order.begin(), order.end(), [&m](int a, int b){ return m.weights[a] > m.weights[b]; } );
      SFDGMMModel s = m;
      for(int k=0;k<K;k++)
      {
        const int o = order[k];
        s.weights[k] = m.weights[o];
        std::copy( m.means.begin() + size_t(o)*D, m.means.begin() + size_t(o+1)*D, s.means.begin() + size_t(k)*D );
        std::copy( m.covariances.begin() + size_t(o)*D*D, m.covariances.begin() + size_t(o+1)*D*D, s.covariances.begin() + size_t(k)*D*D );
      }
      update_precisions( s );
      m = std::move( s );
    }
  }

  // ===================================================================================
  // sfd_gmm_fit
  // ===================================================================================
  template<class GridT>
  class SFDGMMFit : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;
    using StringVector = std::vector<std::string>;
    using LongVector = std::vector<long>;

    ADD_SLOT( MPI_Comm     , mpi            , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT        , grid           , INPUT , REQUIRED );
    ADD_SLOT( DoubleVector , sfd_invariants , INPUT , REQUIRED , DocString{"Per-particle strain functional descriptors (compute_strain_functional)"} );
    ADD_SLOT( StringVector , features       , INPUT , StringVector{} , DocString{"Descriptor columns used as features; empty for the 23 invariants P1_I0 ... P4_I11 (no w, no orientation factors)"} );
    ADD_SLOT( LongVector   , n_components   , INPUT , LongVector{2} , DocString{"Number of mixture components; a list fits each value and keeps the best one according to select_by"} );
    ADD_SLOT( std::string  , select_by      , INPUT , std::string("dbi") , DocString{"Model selection when several n_components are given: dbi (lowest Davies-Bouldin), chi (highest Calinski-Harabasz) or bic (lowest)"} );
    ADD_SLOT( bool         , standardize    , INPUT , false , DocString{"Center and scale each feature by its global mean and standard deviation before fitting"} );
    ADD_SLOT( long         , max_iter       , INPUT , 100 , DocString{"Maximum number of EM iterations"} );
    ADD_SLOT( double       , tol            , INPUT , 1.0e-3 , DocString{"Convergence threshold on the change of the mean log-likelihood"} );
    ADD_SLOT( double       , reg_covar      , INPUT , 1.0e-6 , DocString{"Added to the covariance diagonals"} );
    ADD_SLOT( long         , init_sample    , INPUT , 20000 , DocString{"Number of particles (chosen from their ids) used by the k-means initialization"} );
    ADD_SLOT( long         , n_init         , INPUT , 1 , DocString{"Number of initializations; the fit with the highest likelihood is kept"} );
    ADD_SLOT( long         , seed           , INPUT , 0 , DocString{"Seed of the initialization"} );
    ADD_SLOT( std::string  , model_file     , INPUT , std::string("") , DocString{"If not empty, the selected model is written to this YAML file (input of sfd_gmm_classify model_file)"} );
    ADD_SLOT( std::string  , init_model_file , INPUT , std::string("") , DocString{"If not empty, the initial model (before EM) of the selected fit is written to this YAML file (reproducibility checks)"} );
    ADD_SLOT( std::string  , scores_file    , INPUT , std::string("") , DocString{"If not empty, the scores of all fits are written to this CSV file"} );
    ADD_SLOT( SFDGMMModel  , sfd_gmm_model  , OUTPUT , DocString{"Fitted Gaussian mixture model"} );

  public:
    inline void execute () override final
    {
      using namespace sfd_gmm_detail;
      MPI_Comm comm = *mpi;
      int rank = 0;
      MPI_Comm_rank( comm, &rank );

      const StringVector feat = features->empty() ? default_features() : *features;
      const int D = static_cast<int>( feat.size() );
      if( D < 1 || D > sfd::GMM_MAX_FEATURES ) fatal_error() << "sfd_gmm_fit: 1 to " << sfd::GMM_MAX_FEATURES << " features expected" << std::endl;
      const std::vector<int> cols = resolve_columns( feat );
      if( *select_by != "bic" && *select_by != "chi" && *select_by != "dbi" ) fatal_error() << "sfd_gmm_fit: select_by must be bic, chi or dbi" << std::endl;
      if( n_components->empty() ) fatal_error() << "sfd_gmm_fit: n_components is empty" << std::endl;
      for(long K : *n_components) if( K < 1 || K > sfd::GMM_MAX_COMPONENTS ) fatal_error() << "sfd_gmm_fit: n_components must be in [1," << sfd::GMM_MAX_COMPONENTS << "]" << std::endl;

      // owned particles' feature rows
      LocalData data;
      data.D = D;
      if( sfd_invariants->size() != grid->number_of_particles() * sfd::N_INVARIANTS ) fatal_error() << "sfd_gmm_fit: sfd_invariants does not match the grid (run compute_strain_functional first)" << std::endl;
      {
        auto cells = grid->cells_accessor();
        const size_t * const offset = grid->cell_particle_offset_data();
        const size_t n_cells = grid->number_of_cells();
        for(size_t c=0;c<n_cells;c++)
        {
          if( grid->is_ghost_cell(c) ) continue;
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            const double* row = sfd_invariants->data() + ( offset[c] + p ) * sfd::N_INVARIANTS;
            for(int d=0;d<D;d++) data.x.push_back( row[cols[d]] );
            data.id.push_back( cells[c][field::id][p] );
          }
        }
        data.n = data.id.size();
      }
      double ntot = double( data.n );
      MPI_Allreduce( MPI_IN_PLACE, &ntot, 1, MPI_DOUBLE, MPI_SUM, comm );
      long kmax = 0; for(long K : *n_components) kmax = std::max( kmax, K );
      if( ntot < double(kmax) ) fatal_error() << "sfd_gmm_fit: fewer particles than components" << std::endl;

      // standardization
      std::vector<double> center( D, 0.0 ), scale( D, 1.0 );
      if( *standardize )
      {
        std::vector<double> s( 2*D, 0.0 );
        for(size_t i=0;i<data.n;i++) for(int d=0;d<D;d++) s[d] += data.x[i*D+d];
        MPI_Allreduce( MPI_IN_PLACE, s.data(), D, MPI_DOUBLE, MPI_SUM, comm );
        for(int d=0;d<D;d++) center[d] = s[d] / ntot;
        for(size_t i=0;i<data.n;i++) for(int d=0;d<D;d++) { const double t = data.x[i*D+d] - center[d]; s[D+d] += t*t; }
        MPI_Allreduce( MPI_IN_PLACE, s.data()+D, D, MPI_DOUBLE, MPI_SUM, comm );
        for(int d=0;d<D;d++) { const double sd = std::sqrt( s[D+d] / ntot ); scale[d] = ( sd > 0.0 ) ? sd : 1.0; }
        for(size_t i=0;i<data.n;i++) for(int d=0;d<D;d++) data.x[i*D+d] = ( data.x[i*D+d] - center[d] ) / scale[d];
      }

      // initialization sample: particles selected by a hash of their id (independent of the domain decomposition)
      std::vector<double> sample;
      {
        const double frac = std::min( 1.0 , double( std::max( *init_sample , kmax ) ) / ntot );
        const double thr = frac * 18446744073709551615.0;
        std::vector<uint64_t> sid; std::vector<double> sx;
        for(size_t i=0;i<data.n;i++)
        {
          if( frac >= 1.0 || double( splitmix64( data.id[i] ^ splitmix64( uint64_t(*seed) ) ) ) < thr )
          {
            sid.push_back( data.id[i] );
            sx.insert( sx.end(), data.x.begin() + i*D, data.x.begin() + (i+1)*D );
          }
        }
        int nloc = int( sid.size() ), nranks = 1;
        MPI_Comm_size( comm, &nranks );
        std::vector<int> counts( nranks ), displs( nranks, 0 );
        MPI_Gather( &nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm );
        int nall = 0;
        if( rank == 0 ) { for(int r=0;r<nranks;r++) { displs[r] = nall; nall += counts[r]; } }
        std::vector<uint64_t> all_id( nall );
        MPI_Gatherv( sid.data(), nloc, MPI_UINT64_T, all_id.data(), counts.data(), displs.data(), MPI_UINT64_T, 0, comm );
        std::vector<int> countsD( nranks ), displsD( nranks );
        for(int r=0;r<nranks;r++) { countsD[r] = counts[r]*D; displsD[r] = displs[r]*D; }
        std::vector<double> all_x( size_t(nall)*D );
        MPI_Gatherv( sx.data(), nloc*D, MPI_DOUBLE, all_x.data(), countsD.data(), displsD.data(), MPI_DOUBLE, 0, comm );
        if( rank == 0 )
        {
          std::vector<int> order( nall );
          std::iota( order.begin(), order.end(), 0 );
          std::sort( order.begin(), order.end(), [&all_id](int a, int b){ return all_id[a] < all_id[b]; } );
          sample.resize( size_t(nall)*D );
          for(int i=0;i<nall;i++) std::copy( all_x.begin() + size_t(order[i])*D, all_x.begin() + size_t(order[i]+1)*D, sample.begin() + size_t(i)*D );
          if( nall < kmax ) fatal_error() << "sfd_gmm_fit: initialization sample smaller than n_components, increase init_sample" << std::endl;
        }
      }

      std::vector<FitScores> all_scores;
      SFDGMMModel best_model, best_init;
      FitScores best_sc;
      bool have_best = false;
      std::vector<double> resp;
      std::vector<int> labels( data.n );

      for(long Kl : *n_components)
      {
        const int K = static_cast<int>( Kl );
        SFDGMMModel fit_best, fit_init_best;
        FitScores fit_sc;
        double fit_ll = -std::numeric_limits<double>::infinity();
        for(long init=0; init<std::max(1L,*n_init); init++)
        {
          SFDGMMModel m;
          m.features = feat;
          m.columns.assign( cols.begin(), cols.end() );
          m.center.assign( center.begin(), center.end() );
          m.scale.assign( scale.begin(), scale.end() );

          // k-means centers on rank 0, broadcast, hard assignment of all particles -> initial M-step
          std::vector<double> c( size_t(K)*D );
          if( rank == 0 ) c = kmeans( sample, D, K, splitmix64( uint64_t(*seed) * 1000003ULL + uint64_t(init) * 7919ULL + uint64_t(K) ) );
          MPI_Bcast( c.data(), K*D, MPI_DOUBLE, 0, comm );
          resp.assign( data.n * K, 0.0 );
#         pragma omp parallel for schedule(static)
          for(size_t i=0;i<data.n;i++)
          {
            int best = 0; double bd = std::numeric_limits<double>::max();
            for(int k=0;k<K;k++) { double dd=0; for(int d=0;d<D;d++){ const double t=data.x[i*D+d]-c[size_t(k)*D+d]; dd+=t*t; } if( dd < bd ) { bd = dd; best = k; } }
            resp[i*K+best] = 1.0;
          }
          m_step( data, resp, K, *reg_covar, comm, m );
          SFDGMMModel init_model = m;

          // EM (scikit-learn convention: E-step with the current parameters, then M-step)
          double lower = -std::numeric_limits<double>::infinity();
          FitScores sc; sc.K = K;
          for(long it=1; it<=*max_iter; it++)
          {
            const double prev = lower;
            lower = e_step( data, m, &resp, nullptr, ntot, comm );
            m_step( data, resp, K, *reg_covar, comm, m );
            sc.n_iter = int(it);
            if( std::fabs( lower - prev ) < *tol ) { sc.converged = true; break; }
          }
          if( lower > fit_ll || init == 0 )
          {
            fit_ll = lower;
            fit_best = m;
            fit_init_best = init_model;
            fit_sc = sc;
          }
        }

        sort_components( fit_best );
        sort_components( fit_init_best );
        fit_sc.log_likelihood = e_step( data, fit_best, nullptr, &labels, ntot, comm );
        const double n_params = K * ( D*(D+1)/2.0 ) + K*D + K - 1;
        fit_sc.bic = -2.0 * fit_sc.log_likelihood * ntot + n_params * std::log( ntot );
        cluster_scores( data, labels, K, ntot, comm, fit_sc );
        all_scores.push_back( fit_sc );

        bool better = !have_best;
        if( have_best )
        {
          if( *select_by == "bic" ) better = fit_sc.bic < best_sc.bic;
          else if( *select_by == "chi" ) better = fit_sc.chi > best_sc.chi;
          else better = fit_sc.dbi < best_sc.dbi;
        }
        if( better ) { best_model = fit_best; best_init = fit_init_best; best_sc = fit_sc; have_best = true; }
      }

      lout << "sfd_gmm_fit: " << size_t(ntot) << " particles, " << D << " features" << ( *standardize ? " (standardized)" : "" ) << std::endl;
      lout << "  n_components  iterations  converged  log-likelihood            BIC        CHI        DBI" << std::endl;
      for(const auto& s : all_scores)
      {
        lout << "  " << std::setw(12) << s.K << std::setw(12) << s.n_iter << std::setw(11) << ( s.converged ? "yes" : "no" )
             << std::setw(16) << std::setprecision(8) << s.log_likelihood << std::setw(15) << std::setprecision(8) << s.bic
             << std::setw(11) << std::setprecision(5) << s.chi << std::setw(11) << std::setprecision(5) << s.dbi
             << ( s.K == best_sc.K ? "  <- selected" : "" ) << std::endl;
      }
      lout << "  weights:";
      for(int k=0;k<best_model.n_components();k++) lout << " " << std::setprecision(5) << best_model.weights[k];
      lout << std::endl;

      if( rank == 0 )
      {
        if( ! model_file->empty() ) write_model( best_model, *model_file );
        if( ! init_model_file->empty() ) write_model( best_init, *init_model_file );
        if( ! scores_file->empty() )
        {
          std::ofstream f( *scores_file );
          f << "n_components,n_iter,converged,log_likelihood,bic,chi,dbi\n" << std::setprecision(17);
          for(const auto& s : all_scores) f << s.K << "," << s.n_iter << "," << int(s.converged) << "," << s.log_likelihood << "," << s.bic << "," << s.chi << "," << s.dbi << "\n";
        }
      }
      *sfd_gmm_model = std::move( best_model );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Fits a Gaussian mixture model (full covariances, expectation-maximization) to the strain
functional descriptors of all particles, as in the unsupervised defect analysis of Kober et
al. (arXiv:2402.04191, section 3.6). Same algorithm and defaults as scikit-learn's
GaussianMixture: k-means initialization (here on a sample of particles chosen from their ids,
so the fit does not depend on the number of MPI ranks), reg_covar 1e-6, tol 1e-3, at most 100
iterations. Components are sorted by decreasing weight: class 0 is the most populated one.

With several n_components, each one is fitted and scored (BIC, Calinski-Harabasz CHI,
Davies-Bouldin DBI, computed on the hard assignment) and the best according to select_by
is kept (default dbi). The paper chose the number of components where CHI and DBI both
reach a local extremum. BIC tends to split the thermally broadened bulk into several
components.

Usage example:

compute_strain_functional: { rcut: 5.1 ang, sigma: 0.91 ang }
sfd_gmm_fit: { n_components: [ 2, 3, 4, 5, 6 ], select_by: dbi, model_file: gmm_model.yaml }
sfd_gmm_classify: { class_field: sfd_class }

)EOF";
    }
  };

  // ===================================================================================
  // sfd_gmm_classify
  // ===================================================================================
  struct SFDGMMClassifyFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const double * const __restrict__ m_invariants = nullptr;
    const sfd::GMMView m_gmm = {};
    const double m_min_probability = 0.0;

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t cell, unsigned int part, double& cls, double& prob ) const
    {
      const double* row = m_invariants + ( m_cell_particle_offset[cell] + part ) * sfd::N_INVARIANTS;
      double post = 0.0;
      const int k = sfd::gmm_classify( m_gmm, row, post );
      cls = ( post >= m_min_probability ) ? double(k) : -1.0;
      prob = post;
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits< exaStamp::SFDGMMClassifyFunctor > { static inline constexpr bool CudaCompatible = true; };
}

namespace exaStamp
{
  template<class GridT>
  class SFDGMMClassify : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;
    using LongVector = std::vector<long>;

    ADD_SLOT( MPI_Comm     , mpi               , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT        , grid              , INPUT_OUTPUT );
    ADD_SLOT( DoubleVector , sfd_invariants    , INPUT , REQUIRED , DocString{"Per-particle strain functional descriptors (compute_strain_functional)"} );
    ADD_SLOT( SFDGMMModel  , sfd_gmm_model     , INPUT , OPTIONAL , DocString{"Model from sfd_gmm_fit, used when model_file is empty"} );
    ADD_SLOT( std::string  , model_file        , INPUT , std::string("") , DocString{"YAML model written by sfd_gmm_fit or scripts/python/sfd_gmm_train.py"} );
    ADD_SLOT( std::string  , class_field       , INPUT , std::string("sfd_class") , DocString{"Output per-particle class (component index, -1 if below min_probability)"} );
    ADD_SLOT( std::string  , probability_field , INPUT , std::string("sfd_class_prob") , DocString{"Output per-particle posterior probability of the assigned component"} );
    ADD_SLOT( double       , min_probability   , INPUT , 0.0 , DocString{"Particles whose highest posterior probability is below this value get class -1"} );
    ADD_SLOT( LongVector   , sfd_class_counts  , OUTPUT , DocString{"Number of particles of each class over all ranks (last entry: class -1)"} );

    sfd::SFDGMMModel m_file_model;
    std::string m_loaded_file;

  public:
    inline void execute () override final
    {
      const sfd::SFDGMMModel* model = nullptr;
      if( ! model_file->empty() )
      {
        if( *model_file != m_loaded_file ) { m_file_model = sfd_gmm_detail::read_model( *model_file ); m_loaded_file = *model_file; }
        model = &m_file_model;
      }
      else if( sfd_gmm_model.has_value() ) { model = sfd_gmm_model.get_pointer(); }
      else { fatal_error() << "sfd_gmm_classify: no model (run sfd_gmm_fit or set model_file)" << std::endl; }

      if( sfd_invariants->size() != grid->number_of_particles() * sfd::N_INVARIANTS ) fatal_error() << "sfd_gmm_classify: sfd_invariants does not match the grid (run compute_strain_functional first)" << std::endl;

      for(const std::string* f : { class_field.get_pointer() , probability_field.get_pointer() })
      {
        if( f->size() > 15 ) fatal_error() << "sfd_gmm_classify: field name '" << *f << "' longer than 15 characters" << std::endl;
      }
      auto cls_acc = grid->field_accessor( field::mk_generic_real( *class_field ) );
      auto prob_acc = grid->field_accessor( field::mk_generic_real( *probability_field ) );
      SFDGMMClassifyFunctor func = { grid->cell_particle_offset_data() , sfd_invariants->data() , sfd::gmm_view( *model ) , *min_probability };
      compute_cell_particles( *grid , false , func , onika::make_flat_tuple( cls_acc , prob_acc ) , parallel_execution_context() );

      const int K = model->n_components();
      std::vector<long> counts( K+1, 0 );
      {
        auto cells = grid->cells_accessor();
        const size_t n_cells = grid->number_of_cells();
        for(size_t c=0;c<n_cells;c++)
        {
          if( grid->is_ghost_cell(c) ) continue;
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            const int k = static_cast<int>( cells[c][cls_acc][p] );
            ++ counts[ ( k >= 0 && k < K ) ? k : K ];
          }
        }
      }
      MPI_Allreduce( MPI_IN_PLACE, counts.data(), K+1, MPI_LONG, MPI_SUM, *mpi );
      lout << "sfd_gmm_classify:";
      for(int k=0;k<K;k++) lout << " class " << k << ": " << counts[k];
      if( counts[K] > 0 ) lout << ", unassigned: " << counts[K];
      lout << std::endl;
      *sfd_class_counts = counts;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Assigns each particle to the most probable component of a Gaussian mixture model of its
strain functional descriptors (model from sfd_gmm_fit, or a YAML file from sfd_gmm_fit
model_file or scripts/python/sfd_gmm_train.py). Writes the class and its posterior
probability to per-particle fields. GPU compatible. To keep the same class numbering over a
run, fit once (or offline) and classify every snapshot with the saved model_file.

Usage example:

compute_strain_functional: { rcut: 5.1 ang, sigma: 0.91 ang }
sfd_gmm_classify: { model_file: gmm_model.yaml, min_probability: 0.5 }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(sfd_gmm)
  {
    OperatorNodeFactory::instance()->register_factory( "sfd_gmm_fit", make_grid_variant_operator< SFDGMMFit > );
    OperatorNodeFactory::instance()->register_factory( "sfd_gmm_classify", make_grid_variant_operator< SFDGMMClassify > );
  }
}
