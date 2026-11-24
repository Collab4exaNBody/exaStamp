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
#include <onika/math/basic_types.h>
#include <onika/memory/allocator.h>
#include <onika/string_utils.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <algorithm>
#include <cmath>
#include <vector>
#include <string>

namespace exaStamp
{
  using namespace exanb;

  // maximum degree l : bounds the per-particle q_lm accumulators (a per-thread local)
  static inline constexpr unsigned int STEINHARDT_LMAX = 12;
  // maximum number of nearest neighbors kept in nnn mode (as compute_local_csp)
  static inline constexpr unsigned int STEINHARDT_MAX_NNN = 32;
  // number of (l,m) pairs with 0 <= m <= l <= LMAX, stored at triangular index l(l+1)/2+m
  static inline constexpr unsigned int STEINHARDT_NLM = ( STEINHARDT_LMAX + 1 ) * ( STEINHARDT_LMAX + 2 ) / 2;

  ONIKA_HOST_DEVICE_FUNC inline constexpr unsigned int steinhardt_lm_index(unsigned int l, unsigned int m) { return l * ( l + 1 ) / 2 + m; }

  // per-particle context
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) LocalSteinhardtExtStorage
  {
    // sum over neighbors of Y_lm, m >= 0 only (q_l,-m = (-1)^m conj(q_lm))
    // shortcut: accumulates every l <= max requested degree, not only the requested ones; trim if GPU local memory hurts
    double m_re[STEINHARDT_NLM];
    double m_im[STEINHARDT_NLM];
    // nnn mode : the (at most nnn) nearest neighbors seen so far, sorted by increasing distance
    double m_d2[STEINHARDT_MAX_NNN];
    double m_rx[STEINHARDT_MAX_NNN];
    double m_ry[STEINHARDT_MAX_NNN];
    double m_rz[STEINHARDT_MAX_NNN];
    unsigned int m_count = 0;
  };

  // Per-atom Steinhardt bond-orientational order parameters (Steinhardt, Nelson, Ronchetti, PRB 28, 1983),
  // same definitions as LAMMPS compute orientorder/atom :
  //   q_lm = 1/N sum_j Y_lm(r_ij/|r_ij|) , Q_l = sqrt( 4pi/(2l+1) sum_m |q_lm|^2 )
  //   W_l = sum_{m1+m2+m3=0} (l l l ; m1 m2 m3) q_lm1 q_lm2 q_lm3 , What_l = W_l / ( sum_m |q_lm|^2 )^(3/2)
  // Spherical harmonics without trigonometry nor complex type : Y_lm(u) = P_lm(z) (x+iy)^m, where P_lm is the
  // orthonormalized associated Legendre function divided by sin^m(theta) (a polynomial in z), obtained by the
  // standard 3-term recurrence with precomputed coefficients, and (x+iy)^m by a complex product recurrence.
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) LocalSteinhardtFunctor
  {
    const double * __restrict__ m_rec_a = nullptr;   // [NLM] recurrence coefficient on z*P_{l-1,m} (diagonal factor when m==l)
    const double * __restrict__ m_rec_c = nullptr;   // [NLM] recurrence coefficient on P_{l-2,m}
    const double * __restrict__ m_w3j = nullptr;     // (l l l ; m1 m2 -m1-m2) at m_w3j_offset[k] + (m1+l)*(2l+1) + (m2+l)
    unsigned int m_nnn = 12;                         // 0 = all neighbors within rcut
    unsigned int m_lmax = 6;                         // max requested degree
    unsigned int m_ndeg = 0;
    unsigned int m_deg[STEINHARDT_LMAX+1] = {};
    unsigned int m_w3j_offset[STEINHARDT_LMAX+1] = {};
    bool m_wl = false;
    bool m_wl_hat = false;
    // output flat arrays, indexed by cell_particle_offset[cell]+p (raw pointers : 39 field accessors exceed the functor size limit)
    const size_t * __restrict__ m_cell_particle_offset = nullptr;
    double * m_q[STEINHARDT_LMAX+1] = {};
    double * m_w[STEINHARDT_LMAX+1] = {};
    double * m_what[STEINHARDT_LMAX+1] = {};

    template<class ExtT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void accumulate_ylm(ExtT & ext, double rx, double ry, double rz, double d2) const
    {
      const double inv_r = 1.0 / sqrt( d2 );
      const double x = rx * inv_r, y = ry * inv_r, z = rz * inv_r;
      double cm = 1.0, sm = 0.0;                         // (x+iy)^m
      double pmm = 0.28209479177387814347403972578;      // P_00 = 1/sqrt(4pi)
      for(unsigned int m=0;m<=m_lmax;m++)
      {
        if( m > 0 )
        {
          pmm *= m_rec_a[ steinhardt_lm_index(m,m) ];
          const double c = x * cm - y * sm;
          sm = x * sm + y * cm;
          cm = c;
        }
        unsigned int i = steinhardt_lm_index(m,m);
        ext.m_re[i] += pmm * cm;
        ext.m_im[i] += pmm * sm;
        double p2 = 0.0, p1 = pmm;
        for(unsigned int l=m+1;l<=m_lmax;l++)
        {
          i = steinhardt_lm_index(l,m);
          const double p = m_rec_a[i] * z * p1 - m_rec_c[i] * p2;
          ext.m_re[i] += p * cm;
          ext.m_im[i] += p * sm;
          p2 = p1; p1 = p;
        }
      }
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, LocalCellsT, size_t, size_t, exanb::ComputePairParticleContextStart) const
    {
      ctx.ext.m_count = 0;
      const unsigned int nlm = steinhardt_lm_index(m_lmax+1,0);
      for(unsigned int i=0;i<nlm;i++) { ctx.ext.m_re[i] = 0.0; ctx.ext.m_im[i] = 0.0; }
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
        ComputeBufferT& ctx
      , const Vec3d& dr , double d2
      , LocalCellsT , size_t , size_t
      , double ) const
    {
      auto & ext = ctx.ext;
      if( ! ( d2 > 0.0 ) ) return; // no direction for a superimposed atom
      if( m_nnn == 0 )
      {
        accumulate_ylm( ext, dr.x, dr.y, dr.z, d2 );
        ++ ext.m_count;
        return;
      }
      // nnn mode : keep the nnn nearest by insertion (same selection as compute_local_csp)
      unsigned int pos = 0;
      if( ext.m_count < m_nnn ) { pos = ext.m_count ++ ; }
      else if( d2 < ext.m_d2[m_nnn-1] ) { pos = m_nnn - 1; } // replaces the farthest kept neighbor
      else return;                                            // ties keep the first neighbor seen
      while( pos > 0 && ext.m_d2[pos-1] > d2 )
      {
        ext.m_d2[pos] = ext.m_d2[pos-1];
        ext.m_rx[pos] = ext.m_rx[pos-1];
        ext.m_ry[pos] = ext.m_ry[pos-1];
        ext.m_rz[pos] = ext.m_rz[pos-1];
        -- pos;
      }
      ext.m_d2[pos] = d2;
      ext.m_rx[pos] = dr.x;
      ext.m_ry[pos] = dr.y;
      ext.m_rz[pos] = dr.z;
    }

    template<class ComputeBufferT, class LocalCellsT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, LocalCellsT, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop) const
    {
      auto & ext = ctx.ext;
      const size_t ip = m_cell_particle_offset[cell_a] + p_a;
      // nnn mode with fewer than nnn neighbors within rcut, or no neighbor at all : every output is 0 (as LAMMPS)
      const bool valid = ( m_nnn == 0 ) ? ( ext.m_count > 0 ) : ( ext.m_count >= m_nnn );
      if( valid && m_nnn > 0 )
      {
        for(unsigned int j=0;j<m_nnn;j++) accumulate_ylm( ext, ext.m_rx[j], ext.m_ry[j], ext.m_rz[j], ext.m_d2[j] );
      }
      const double inv_n = valid ? 1.0 / ext.m_count : 0.0;

      for(unsigned int k=0;k<m_ndeg;k++)
      {
        const unsigned int l = m_deg[k];
        const int il = l;
        // full q_lm, m = -l..l at index m+l
        double qr[2*STEINHARDT_LMAX+1];
        double qi[2*STEINHARDT_LMAX+1];
        double sum2 = 0.0;
        for(int m=0;m<=il;m++)
        {
          const unsigned int i = steinhardt_lm_index(l,m);
          const double re = ext.m_re[i] * inv_n;
          const double im = ext.m_im[i] * inv_n;
          const double sign = ( m % 2 == 0 ) ? 1.0 : -1.0;
          qr[il+m] = re;         qi[il+m] = im;
          qr[il-m] = sign * re;  qi[il-m] = - sign * im;
          sum2 += ( m == 0 ? 1.0 : 2.0 ) * ( re * re + im * im );
        }
        m_q[k][ip] = sqrt( 4.0 * M_PI / ( 2 * l + 1 ) * sum2 );

        if( m_wl || m_wl_hat )
        {
          const double * __restrict__ w3j = m_w3j + m_w3j_offset[k];
          double wl = 0.0;
          for(int m1=-il;m1<=il;m1++)
          {
            const int m2_min = ( -il > -il-m1 ) ? -il : -il-m1;
            const int m2_max = (  il <  il-m1 ) ?  il :  il-m1;
            for(int m2=m2_min;m2<=m2_max;m2++)
            {
              const int m3 = - m1 - m2;
              // real part of q_m1 q_m2 q_m3
              const double ar = qr[il+m1] * qr[il+m2] - qi[il+m1] * qi[il+m2];
              const double ai = qr[il+m1] * qi[il+m2] + qi[il+m1] * qr[il+m2];
              wl += w3j[ (m1+il) * (2*il+1) + (m2+il) ] * ( ar * qr[il+m3] - ai * qi[il+m3] );
            }
          }
          if( m_wl ) m_w[k][ip] = wl;
          if( m_wl_hat ) m_what[k][ip] = ( sum2 > 0.0 ) ? wl / ( sum2 * sqrt( sum2 ) ) : 0.0;
        }
      }
    }
  };

  // Wigner 3j symbol (l l l ; m1 m2 m3) by the Racah formula, host only (factorials up to (3l+1)! = 37!, fine in double)
  static inline double steinhardt_wigner3j_lll(int l, int m1, int m2, int m3)
  {
    if( m1 + m2 + m3 != 0 || std::abs(m1) > l || std::abs(m2) > l || std::abs(m3) > l ) return 0.0;
    auto f = [](int n) { return std::tgamma( n + 1.0 ); };
    // (j1 j2 j3 ; m1 m2 m3) with j1=j2=j3=l
    const double delta = f(l) * f(l) * f(l) / f(3*l+1);
    const double pref = std::sqrt( delta * f(l+m1) * f(l-m1) * f(l+m2) * f(l-m2) * f(l+m3) * f(l-m3) );
    // k range where every factorial argument below is >= 0
    const int kmin = std::max( { 0 , -m1 , m2 } );
    const int kmax = std::min( { l , l - m1 , l + m2 } );
    double s = 0.0;
    for(int k=kmin;k<=kmax;k++)
    {
      const double t = 1.0 / ( f(k) * f(k+m1) * f(k-m2) * f(l-k) * f(l-k-m1) * f(l-k+m2) );
      s += ( k % 2 == 0 ) ? t : -t;
    }
    const double sign = ( ( ( -m3 ) % 2 ) == 0 ) ? 1.0 : -1.0; // (-1)^(j1-j2-m3) = (-1)^(-m3)
    return sign * pref * s;
  }

  template<class GridT>
  class ComputeLocalSteinhardt : public OperatorNode
  {
    using DoubleVector = onika::memory::CudaMMVector<double>;

    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT , DocString{"Local sub-domain particles grid"} );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED , DocString{"Simulation domain"} );
    ADD_SLOT( double                    , rcut            , INPUT        , REQUIRED , DocString{"Neighbor search radius (must enclose the nnn nearest neighbors in nnn mode)"} );
    ADD_SLOT( size_t                    , nnn             , INPUT        , 12 , DocString{"Number of nearest neighbors, in [1;32] (12 for FCC, 8 or 14 for BCC). 0 = all neighbors within rcut"} );
    ADD_SLOT( std::vector<long>         , degrees         , INPUT        , std::vector<long>{4,6} , DocString{"Degrees l to compute, distinct, in [0;12]"} );
    ADD_SLOT( bool                      , wl              , INPUT        , false , DocString{"Also compute the third-order invariants W_l"} );
    ADD_SLOT( bool                      , wl_hat          , INPUT        , false , DocString{"Also compute the normalized third-order invariants What_l"} );
    ADD_SLOT( std::string               , name            , INPUT        , std::string("") , DocString{"Prefix (max 11 characters) of the per-particle scalar fields <name>q<l>, <name>w<l>, <name>wh<l>"} );
    ADD_SLOT( bool                      , ghost           , INPUT        , false , DocString{"Also compute for ghost particles"} );
    ADD_SLOT( bool                      , verbose         , INPUT        , false , DocString{"Print parameters"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT        , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double                    , rcut_max        , INPUT        , 0.0 , DocString{"Neighbor list cutoff (global rcut_max), must be >= rcut"} );

    DoubleVector m_rec_a;
    DoubleVector m_rec_c;
    DoubleVector m_w3j;

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;
      if( chunk_neighbors->number_of_cells() != grid->number_of_cells() )
      {
        fatal_error() << "compute_local_steinhardt: neighbor list does not match the grid (missing chunk_neighbors ?)" << std::endl;
      }
      if( ! ( *rcut > 0.0 ) ) { fatal_error() << "compute_local_steinhardt: rcut must be > 0 (got " << *rcut << ")" << std::endl; }
      // the neighbor list is sized from the global rcut_max at initialization, before analysis operators run
      if( *rcut > *rcut_max )
      {
        fatal_error() << "compute_local_steinhardt: rcut=" << *rcut << " exceeds the neighbor list cutoff rcut_max=" << *rcut_max
                      << ", neighbors would be missing. Set 'rcut_max' >= rcut in the global section of the input file" << std::endl;
      }
      if( *nnn > STEINHARDT_MAX_NNN )
      {
        fatal_error() << "compute_local_steinhardt: nnn=" << *nnn << " must be in [1;" << STEINHARDT_MAX_NNN << "] (or 0 = all neighbors within rcut)" << std::endl;
      }
      if( degrees->empty() || degrees->size() > STEINHARDT_LMAX + 1 )
      {
        fatal_error() << "compute_local_steinhardt: degrees must hold between 1 and " << STEINHARDT_LMAX + 1 << " values" << std::endl;
      }
      // dynamic field names are truncated to NAME_MAX_LEN-1 characters : longest one is <name>wh12
      if( name->size() + 4 >= field::generic_real::NAME_MAX_LEN )
      {
        fatal_error() << "compute_local_steinhardt: name prefix '" << *name << "' is too long (max " << field::generic_real::NAME_MAX_LEN - 5 << " characters)" << std::endl;
      }
      bool seen[STEINHARDT_LMAX+1] = {};
      for(long l : *degrees)
      {
        if( l < 0 || l > long(STEINHARDT_LMAX) || seen[l] )
        {
          fatal_error() << "compute_local_steinhardt: degree " << l << " must be distinct and in [0;" << STEINHARDT_LMAX << "]" << std::endl;
        }
        seen[l] = true;
      }

      LocalSteinhardtFunctor compute_op = {};
      compute_op.m_cell_particle_offset = grid->cell_particle_offset_data();
      compute_op.m_nnn = *nnn;
      compute_op.m_ndeg = degrees->size();
      compute_op.m_wl = *wl;
      compute_op.m_wl_hat = *wl_hat;
      compute_op.m_lmax = 0;
      for(size_t k=0;k<degrees->size();k++)
      {
        const long l = degrees->at(k);
        compute_op.m_deg[k] = l;
        if( unsigned(l) > compute_op.m_lmax ) compute_op.m_lmax = l;
        const std::string ls = std::to_string(l);
        compute_op.m_q[k] = grid->flat_array_data( field::mk_generic_real( *name + "q" + ls ) );
        if( *wl ) compute_op.m_w[k] = grid->flat_array_data( field::mk_generic_real( *name + "w" + ls ) );
        if( *wl_hat ) compute_op.m_what[k] = grid->flat_array_data( field::mk_generic_real( *name + "wh" + ls ) );
      }

      // recurrence coefficients of the orthonormalized P_lm(z) / sin^m(theta) :
      //   P_mm = a_mm P_{m-1,m-1} , a_mm = -sqrt((2m+1)/(2m)) (Condon-Shortley phase)
      //   P_lm = a_lm z P_{l-1,m} - c_lm P_{l-2,m} , a_lm = sqrt((4l^2-1)/(l^2-m^2)) , c_lm = sqrt((2l+1)/(2l-3) ((l-1)^2-m^2)/(l^2-m^2))
      m_rec_a.assign( STEINHARDT_NLM , 0.0 );
      m_rec_c.assign( STEINHARDT_NLM , 0.0 );
      for(unsigned int l=1;l<=STEINHARDT_LMAX;l++)
      {
        for(unsigned int m=0;m<=l;m++)
        {
          const unsigned int i = steinhardt_lm_index(l,m);
          const double ld = l, md = m;
          if( m == l ) { m_rec_a[i] = - std::sqrt( ( 2.0 * md + 1.0 ) / ( 2.0 * md ) ); continue; }
          m_rec_a[i] = std::sqrt( ( 4.0 * ld * ld - 1.0 ) / ( ld * ld - md * md ) );
          if( l >= m + 2 ) m_rec_c[i] = std::sqrt( ( 2.0 * ld + 1.0 ) / ( 2.0 * ld - 3.0 ) * ( ( ld - 1.0 ) * ( ld - 1.0 ) - md * md ) / ( ld * ld - md * md ) );
        }
      }

      // Wigner 3j tables, one (2l+1)^2 block per requested degree
      m_w3j.clear();
      if( *wl || *wl_hat )
      {
        for(size_t k=0;k<degrees->size();k++)
        {
          const int l = degrees->at(k);
          compute_op.m_w3j_offset[k] = m_w3j.size();
          for(int m1=-l;m1<=l;m1++) for(int m2=-l;m2<=l;m2++) m_w3j.push_back( steinhardt_wigner3j_lll( l, m1, m2, -m1-m2 ) );
        }
      }
      compute_op.m_rec_a = m_rec_a.data();
      compute_op.m_rec_c = m_rec_c.data();
      compute_op.m_w3j = m_w3j.data();

      if( *verbose )
      {
        std::string dl;
        for(long l : *degrees) dl += " " + std::to_string(l);
        lout << "compute_local_steinhardt: rcut=" << *rcut << " nnn=" << *nnn << " degrees=" << dl
             << " wl=" << std::boolalpha << *wl << " wl_hat=" << *wl_hat << std::endl;
      }

      using ComputeBuffer = ComputePairBuffer2<false,false,LocalSteinhardtExtStorage>;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr std::integral_constant<bool,true> force_use_cells_accessor = {};
      static constexpr DefaultPositionFields posfields = {};
      compute_cell_particle_pairs2( *grid, *rcut, *ghost, optional, compute_buf, compute_op, compute_field_set
                                   , posfields, parallel_execution_context(), force_use_cells_accessor );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Computes per-atom Steinhardt bond-orientational order parameters (Steinhardt, Nelson, Ronchetti,
Phys. Rev. B 28, 1983; same definitions as LAMMPS compute orientorder/atom) as per-particle
scalar fields, GPU-compatibility to be tested:

  q_lm(i) = 1/N sum_j Y_lm(r_ij/|r_ij|)      over the N neighbors j of atom i
  Q_l     = sqrt( 4 pi/(2l+1) sum_m |q_lm|^2 )                         -> field <name>q<l>
  W_l     = sum_{m1+m2+m3=0} (l l l ; m1 m2 m3) q_lm1 q_lm2 q_lm3     -> field <name>w<l>  (wl: true)
  What_l  = W_l / ( sum_m |q_lm|^2 )^(3/2)                             -> field <name>wh<l> (wl_hat: true)

Neighbors: the nnn nearest within rcut, or all neighbors within rcut when nnn = 0.
Atoms with fewer than nnn neighbors within rcut (or none when nnn = 0) get 0 for every output.
Reference values (perfect lattices): FCC (nnn 12) Q4=0.19094 Q6=0.57452 What6=-0.01316,
BCC (nnn 8) Q4=0.50918 Q6=0.62854, BCC (nnn 14) Q4=0.03637 Q6=0.51069.

Requirements / limits:
  - degrees distinct, in [0;12]; nnn in [1;32] or 0
  - name prefix at most 11 characters
  - rcut must not exceed the neighbor list cutoff: set 'rcut_max: <value>' (>= rcut) in the
    global section when rcut is larger than the potential cutoff

Usage example:

compute_local_steinhardt:
  rcut: 3.0 ang
  nnn: 12
  degrees: [ 4, 6 ]
  wl_hat: true
  name: st         # -> fields stq4, stq6, stwh4, stwh6

)EOF";
    }
  };

  template<class GridT> using ComputeLocalSteinhardtTmpl = ComputeLocalSteinhardt<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_local_steinhardt)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_local_steinhardt", make_grid_variant_operator< ComputeLocalSteinhardtTmpl > );
  }

}

namespace exanb
{
  // specialize functor traits to allow Cuda execution space
  template<>
  struct ComputePairTraits< exaStamp::LocalSteinhardtFunctor >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };
}
