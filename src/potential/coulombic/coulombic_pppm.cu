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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/flat_tuple.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_context.h>
#include <onika/parallel/parallel_for.h>
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <climits>
#include <limits>

#include <exaStamp/potential/coulombic/pppm.h>
#include <exaStamp/potential/coulombic/pppm_fft.h>

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using namespace exanb;

  template<bool PerAtomCharge, class ChargeOrTypeT>
  ONIKA_HOST_DEVICE_FUNC static inline double pppm_particle_charge( const ParticleSpecie* __restrict__ species, ChargeOrTypeT ct )
  {
    if constexpr ( PerAtomCharge ) return ct;
    else return species[ct].m_charge;
  }

  // LAMMPS PPPM::make_rho : charge density on the (global, replicated) mesh.
  // GPU : atomic adds into m_density. CPU with several OpenMP threads (m_thread_density set) : each thread adds into
  // its own mesh (m_nthreads meshes, nfft apart), summed afterwards by PPPMSumThreadMeshesFunc. CPU with one thread :
  // plain adds into m_density.
  template<bool PerAtomCharge>
  struct PPPMSpreadFunc
  {
    ReadOnlyPPPMParameters p;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    double * __restrict__ m_density = nullptr;
    double * __restrict__ m_thread_density = nullptr;
    int m_nthreads = 0;
    size_t m_nfft = 0;

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      const double q = pppm_particle_charge<PerAtomCharge>( m_species , ct );
      if( q == 0.0 ) return;
      PPPMStencil st;
      pppm_stencil( p , Vec3d{rx,ry,rz} , st );
      const double z0 = p.delvolinv * q;
#     ifndef ONIKA_GPU_DEVICE_COMPILE
      if( m_nthreads == 1 || m_thread_density != nullptr )
      {
        double * __restrict__ mesh = m_density;
        if( m_thread_density != nullptr )
        {
          const int t = omp_get_thread_num();
          if( t >= m_nthreads ) onika::fatal_error() << "PPPMSpreadFunc : thread "<<t<<" out of "<<m_nthreads<<" thread meshes"<<std::endl;
          mesh = m_thread_density + size_t(t) * m_nfft;
        }
        for( int n = 0 ; n < p.order ; n++ )
        {
          int mz = st.iz + n; if( mz >= p.mnz ) mz -= p.mnz;
          const double y0 = z0 * st.wz[n];
          for( int m = 0 ; m < p.order ; m++ )
          {
            int my = st.iy + m; if( my >= p.mny ) my -= p.mny;
            const double x0 = y0 * st.wy[m];
            double * __restrict__ row = mesh + ( size_t(mz) * p.mny + my ) * p.mnx;
            for( int l = 0 ; l < p.order ; l++ )
            {
              int mx = st.ix + l; if( mx >= p.mnx ) mx -= p.mnx;
              row[mx] += x0 * st.wx[l];
            }
          }
        }
        return;
      }
#     endif
      for( int n = 0 ; n < p.order ; n++ )
      {
        int mz = st.iz + n; if( mz >= p.mnz ) mz -= p.mnz;
        const double y0 = z0 * st.wz[n];
        for( int m = 0 ; m < p.order ; m++ )
        {
          int my = st.iy + m; if( my >= p.mny ) my -= p.mny;
          const double x0 = y0 * st.wy[m];
          double * __restrict__ row = m_density + ( size_t(mz) * p.mny + my ) * p.mnx;
          for( int l = 0 ; l < p.order ; l++ )
          {
            int mx = st.ix + l; if( mx >= p.mnx ) mx -= p.mnx;
            ONIKA_CU_ATOMIC_ADD( row[mx] , x0 * st.wx[l] );
          }
        }
      }
    }
  };

  // slab correction (LAMMPS slabcorr) : sum of q.z and q.z^2 over local particles, z = real height
  template<bool PerAtomCharge>
  struct PPPMSlabDipoleFunc
  {
    const ParticleSpecie * __restrict__ m_species = nullptr;
    double * __restrict__ m_sum = nullptr; // [ sum q.z , sum q.z^2 ]
    double m_zscale = 1.0;
    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double rz, ChargeOrTypeT ct ) const
    {
      const double q = pppm_particle_charge<PerAtomCharge>( m_species , ct );
      const double z = m_zscale * rz;
      ONIKA_CU_ATOMIC_ADD( m_sum[0] , q * z );
      ONIKA_CU_ATOMIC_ADD( m_sum[1] , q * z * z );
    }
  };

  // real space meshes interpolated back to particles. Two real meshes are stored in one complex mesh (real and
  // imaginary parts), as they come out of one backward FFT (see PPPMPairFactorFunc).
  struct PPPMMeshes
  {
    const Complexd * __restrict__ exy = nullptr; // ik : gradient of the potential x (real) , y (imaginary) ; ad : potential (real)
    const Complexd * __restrict__ ezu = nullptr; // ik : gradient z (real) , potential for per particle energy (imaginary)
    const Complexd * __restrict__ v01 = nullptr; // virial meshes xx , yy
    const Complexd * __restrict__ v23 = nullptr; // virial meshes zz , xy
    const Complexd * __restrict__ v45 = nullptr; // virial meshes xz , yz
  };

  // LAMMPS PPPM::fieldforce_ik or fieldforce_ad (+ fieldforce_peratom and per atom self/background correction)
  template<bool PerAtomCharge, bool ComputeEnergy, bool ComputeVirial>
  struct PPPMForceFunc
  {
    ReadOnlyPPPMParameters p;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    PPPMMeshes m;

    ONIKA_HOST_DEVICE_FUNC inline void compute( double q, double rx, double ry, double rz, Vec3d& f, double& ep, Mat3d& vir ) const
    {
      if( q == 0.0 ) return;
      PPPMStencil st;
      pppm_stencil( p , Vec3d{rx,ry,rz} , st );
      double ekx = 0.0, eky = 0.0, ekz = 0.0;
      double u = 0.0;
      double v0 = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0, v4 = 0.0, v5 = 0.0;
      double sfx = 0.0, sfy = 0.0, sfz = 0.0;
      if( p.diff_ad )
      {
        // field = - gradient of the interpolated potential, through the derivative of the weights
        double dwx[pppm_constants::MAXORDER], dwy[pppm_constants::MAXORDER], dwz[pppm_constants::MAXORDER];
        pppm_stencil_derivative( p , st , dwx , dwy , dwz );
        for( int n = 0 ; n < p.order ; n++ )
        {
          int mz = st.iz + n; if( mz >= p.mnz ) mz -= p.mnz;
          for( int mm = 0 ; mm < p.order ; mm++ )
          {
            int my = st.iy + mm; if( my >= p.mny ) my -= p.mny;
            const size_t row = ( size_t(mz) * p.mny + my ) * p.mnx;
            for( int l = 0 ; l < p.order ; l++ )
            {
              int mx = st.ix + l; if( mx >= p.mnx ) mx -= p.mnx;
              const size_t idx = row + mx;
              const double uval = m.exy[idx].r;
              ekx += dwx[l] * st.wy[mm] * st.wz[n] * uval;
              eky += st.wx[l] * dwy[mm] * st.wz[n] * uval;
              ekz += st.wx[l] * st.wy[mm] * dwz[n] * uval;
              const double x0 = st.wx[l] * st.wy[mm] * st.wz[n];
              if constexpr ( ComputeEnergy ) u += x0 * uval;
              if constexpr ( ComputeVirial )
              {
                const Complexd a = m.v01[idx], b = m.v23[idx], c = m.v45[idx];
                v0 += x0 * a.r;
                v1 += x0 * a.i;
                v2 += x0 * b.r;
                v3 += x0 * b.i;
                v4 += x0 * c.r;
                v5 += x0 * c.i;
              }
            }
          }
        }
        ekx *= p.hinv.x; eky *= p.hinv.y; ekz *= p.hinv.z;
        // self force correction, as LAMMPS on the absolute position x[i]*hx_inv
        const double s1 = rx * p.xf.x * p.hinv.x;
        const double s2 = ry * p.xf.y * p.hinv.y;
        const double s3 = rz * p.xf.z * p.hinv.z;
        sfx = 2.0*q*q * ( p.sf_coeff[0]*sin(2.0*M_PI*s1) + p.sf_coeff[1]*sin(4.0*M_PI*s1) );
        sfy = 2.0*q*q * ( p.sf_coeff[2]*sin(2.0*M_PI*s2) + p.sf_coeff[3]*sin(4.0*M_PI*s2) );
        sfz = 2.0*q*q * ( p.sf_coeff[4]*sin(2.0*M_PI*s3) + p.sf_coeff[5]*sin(4.0*M_PI*s3) );
      }
      else
      for( int n = 0 ; n < p.order ; n++ )
      {
        int mz = st.iz + n; if( mz >= p.mnz ) mz -= p.mnz;
        const double z0 = st.wz[n];
        for( int mm = 0 ; mm < p.order ; mm++ )
        {
          int my = st.iy + mm; if( my >= p.mny ) my -= p.mny;
          const double y0 = z0 * st.wy[mm];
          const size_t row = ( size_t(mz) * p.mny + my ) * p.mnx;
          for( int l = 0 ; l < p.order ; l++ )
          {
            int mx = st.ix + l; if( mx >= p.mnx ) mx -= p.mnx;
            const size_t idx = row + mx;
            const double x0 = y0 * st.wx[l];
            const Complexd exy = m.exy[idx];
            const Complexd ezu = m.ezu[idx];
            ekx -= x0 * exy.r;
            eky -= x0 * exy.i;
            ekz -= x0 * ezu.r;
            if constexpr ( ComputeEnergy ) u += x0 * ezu.i;
            if constexpr ( ComputeVirial )
            {
              const Complexd a = m.v01[idx], b = m.v23[idx], c = m.v45[idx];
              v0 += x0 * a.r;
              v1 += x0 * a.i;
              v2 += x0 * b.r;
              v3 += x0 * b.i;
              v4 += x0 * c.r;
              v5 += x0 * c.i;
            }
          }
        }
      }
      const double qfactor = COULOMB_CONSTANT * q;
      f.x += qfactor * ekx - COULOMB_CONSTANT * sfx;
      f.y += qfactor * eky - COULOMB_CONSTANT * sfy;
      f.z += qfactor * ekz - COULOMB_CONSTANT * sfz;
      if( p.slab )
      {
        const double z = p.zscale * rz;
        f.z += pppm_slab_force_z( p , q , z );
        if constexpr ( ComputeEnergy ) ep += pppm_slab_energy( p , q , z );
      }
      if constexpr ( ComputeEnergy )
      {
        ep += 0.5 * qfactor * u + pppm_self_energy( p , q );
        if constexpr ( ComputeVirial )
        {
          const double h = 0.5 * qfactor;
          vir.m11 += h * v0; vir.m22 += h * v1; vir.m33 += h * v2;
          vir.m12 += h * v3; vir.m21 += h * v3;
          vir.m13 += h * v4; vir.m31 += h * v4;
          vir.m23 += h * v5; vir.m32 += h * v5;
        }
      }
    }

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double & fx, double & fy, double & fz, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      static_assert( ! ComputeEnergy && ! ComputeVirial );
      Vec3d f = {0.,0.,0.}; double ep = 0.0; Mat3d vir;
      compute( pppm_particle_charge<PerAtomCharge>( m_species , ct ) , rx, ry, rz, f, ep, vir );
      fx += f.x; fy += f.y; fz += f.z;
    }

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double & fx, double & fy, double & fz, double & ep, Mat3d & virial, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      static_assert( ComputeEnergy && ComputeVirial );
      Vec3d f = {0.,0.,0.}; double e = 0.0; Mat3d vir = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      compute( pppm_particle_charge<PerAtomCharge>( m_species , ct ) , rx, ry, rz, f, e, vir );
      fx += f.x; fy += f.y; fz += f.z; ep += e; virial += vir;
    }
  };

  // ------------- per mesh point operations (parallel_for, one thread per mesh point on GPU) -------------

  struct PPPMSumThreadMeshesFunc
  {
    const double * __restrict__ m_thread_density = nullptr;
    double * __restrict__ m_density = nullptr;
    int m_nthreads = 0;
    size_t m_nfft = 0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      double s = 0.0;
      for( int t = 0 ; t < m_nthreads ; t++ ) s += m_thread_density[ t * m_nfft + i ];
      m_density[i] = s;
    }
  };

  struct PPPMZeroFunc
  {
    double * __restrict__ a = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { a[i] = 0.0; }
  };

  // w = density (real) , then w *= scale * G(k) once transformed
  struct PPPMLoadDensityFunc
  {
    const double * __restrict__ density = nullptr;
    Complexd * __restrict__ w = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { w[i] = Complexd{ density[i] , 0.0 }; }
  };

  struct PPPMApplyGreenFunc
  {
    Complexd * __restrict__ w = nullptr;
    const double * __restrict__ greensfn = nullptr;
    double scaleinv = 1.0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      const double s = scaleinv * greensfn[i];
      w[i].r *= s;
      w[i].i *= s;
    }
  };

  // Two real meshes A and B from one backward FFT : w2 = H(A(k)) + i.H(B(k)), where H(X)(k) = ( X(k) + conj(X(-k)) )/2
  // is the hermitian part. The backward FFT of a hermitian array is real, so the result holds A in its real part and
  // B in its imaginary part. LAMMPS keeps the real part of each backward FFT, which is exactly the transform of the
  // hermitian part : results are the same, including at Nyquist frequencies where i.k.V(k) is not hermitian.
  // FIELD  : (A,B) = (i.kx.V , i.ky.V) into w2[0] and (i.kz.V , V or 0) into w2[1] (LAMMPS poisson_ik, poisson_peratom u)
  // VIRIAL : (vg_xx.V , vg_yy.V) , (vg_zz.V , vg_xy.V) , (vg_xz.V , vg_yz.V) into w2[0..2] (LAMMPS poisson_peratom v)
  ONIKA_HOST_DEVICE_FUNC static inline Complexd pppm_hermitian_pair( const Complexd& a , const Complexd& am , const Complexd& b , const Complexd& bm )
  {
    const Complexd ha = { 0.5 * ( a.r + am.r ) , 0.5 * ( a.i - am.i ) };
    const Complexd hb = { 0.5 * ( b.r + bm.r ) , 0.5 * ( b.i - bm.i ) };
    return Complexd{ ha.r - hb.i , ha.i + hb.r };
  }
  ONIKA_HOST_DEVICE_FUNC static inline Complexd pppm_ik( double k , const Complexd& c ) { return Complexd{ - k * c.i , k * c.r }; }
  ONIKA_HOST_DEVICE_FUNC static inline Complexd pppm_scale( double s , const Complexd& c ) { return Complexd{ s * c.r , s * c.i }; }

  // POTENTIAL : V into the real part of w2[0] (LAMMPS poisson_ad u_brick)
  enum class PPPMBackwardKind { FIELD , VIRIAL , POTENTIAL };

  // Work item = one mesh point (m_rows false, GPU), or one x row of the mesh (m_rows true, CPU : no index division).
  struct PPPMPairFactorFunc
  {
    const Complexd * __restrict__ w1 = nullptr;     // V(k)
    Complexd * __restrict__ w2[3] = { nullptr , nullptr , nullptr };
    const double * __restrict__ fkx = nullptr;
    const double * __restrict__ fky = nullptr;
    const double * __restrict__ fkz = nullptr;
    double g_ewald = 0.0;
    PPPMBackwardKind kind = PPPMBackwardKind::FIELD;
    bool with_u = false; // FIELD : potential mesh in the imaginary part of w2[1]
    int nx = 0, ny = 0, nz = 0;
    bool m_rows = false;
    const int * __restrict__ krow_partner = nullptr; // distributed layout (l*nx+ix)*nz+iz : local row of -y for each local row

    ONIKA_HOST_DEVICE_FUNC inline void point( size_t i , size_t im ) const
    {
      const Complexd c = w1[i], cm = w1[im];
      if( kind == PPPMBackwardKind::POTENTIAL )
      {
        const Complexd zero = { 0.0 , 0.0 };
        w2[0][i] = pppm_hermitian_pair( c , cm , zero , zero );
      }
      else if( kind == PPPMBackwardKind::FIELD )
      {
        w2[0][i] = pppm_hermitian_pair( pppm_ik( fkx[i] , c ) , pppm_ik( fkx[im] , cm ) , pppm_ik( fky[i] , c ) , pppm_ik( fky[im] , cm ) );
        const Complexd zero = { 0.0 , 0.0 };
        w2[1][i] = pppm_hermitian_pair( pppm_ik( fkz[i] , c ) , pppm_ik( fkz[im] , cm ) , with_u ? c : zero , with_u ? cm : zero );
      }
      else
      {
        double v[6], vm[6];
        pppm_virial_coeffs( fkx[i] , fky[i] , fkz[i] , g_ewald , v );
        pppm_virial_coeffs( fkx[im] , fky[im] , fkz[im] , g_ewald , vm );
        for( int q = 0 ; q < 3 ; q++ )
        {
          w2[q][i] = pppm_hermitian_pair( pppm_scale( v[2*q] , c ) , pppm_scale( vm[2*q] , cm ) , pppm_scale( v[2*q+1] , c ) , pppm_scale( vm[2*q+1] , cm ) );
        }
      }
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      if( krow_partner != nullptr )
      {
        // distributed layout : work item = one z column (CPU, m_rows) or one point (GPU)
        if( m_rows )
        {
          const int li = i / nx;
          const int ix = i % nx;
          const size_t col = size_t(i) * nz;
          const size_t colm = ( size_t( krow_partner[li] ) * nx + ( ix ? nx - ix : 0 ) ) * nz;
          point( col , colm );
          for( int iz = 1 ; iz < nz ; iz++ ) point( col + iz , colm + nz - iz );
        }
        else
        {
          const unsigned int ui = i;
          const unsigned int iz = ui % nz;
          const unsigned int c = ui / nz;
          const unsigned int ix = c % nx;
          const unsigned int li = c / nx;
          point( i , ( size_t( krow_partner[li] ) * nx + ( ix ? nx - ix : 0 ) ) * nz + ( iz ? nz - iz : 0 ) );
        }
        return;
      }
      if( m_rows )
      {
        const int iy = i % ny;
        const int iz = i / ny;
        const size_t row = size_t(i) * nx;
        const size_t rowm = ( size_t( iz ? nz - iz : 0 ) * ny + ( iy ? ny - iy : 0 ) ) * nx; // row of -k
        point( row , rowm );
        for( int ix = 1 ; ix < nx ; ix++ ) point( row + ix , rowm + nx - ix );
      }
      else
      {
        const unsigned int ui = i;
        const unsigned int ix = ui % nx;
        const unsigned int iy = ( ui / nx ) % ny;
        const unsigned int iz = ui / ( unsigned(nx) * ny );
        point( i , ( size_t( iz ? nz - iz : 0 ) * ny + ( iy ? ny - iy : 0 ) ) * nx + ( ix ? nx - ix : 0 ) ); // -k
      }
    }
  };

  // ------------- distributed mesh : brick bounds, pack / unpack through index maps -------------

  // bounds of the local particles' stencils, unwrapped mesh indices : [xmin,ymin,zmin,xmax,ymax,zmax] (inclusive)
  struct PPPMBrickBoundsFunc
  {
    ReadOnlyPPPMParameters p;
    int * __restrict__ m_bounds = nullptr;
    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double rx, double ry, double rz, ChargeOrTypeT ) const
    {
      PPPMStencil st;
      pppm_stencil( p , Vec3d{rx,ry,rz} , st );
      ONIKA_CU_ATOMIC_MIN( m_bounds[0] , st.gx0 );
      ONIKA_CU_ATOMIC_MIN( m_bounds[1] , st.gy0 );
      ONIKA_CU_ATOMIC_MIN( m_bounds[2] , st.gz0 );
      ONIKA_CU_ATOMIC_MAX( m_bounds[3] , st.gx0 + p.order - 1 );
      ONIKA_CU_ATOMIC_MAX( m_bounds[4] , st.gy0 + p.order - 1 );
      ONIKA_CU_ATOMIC_MAX( m_bounds[5] , st.gz0 + p.order - 1 );
    }
  };

  // pack / unpack loops skip this rank's own segment [skip_lo, skip_lo+skip_n) of the index map, copied directly
  // (PPPMCopyAddRealFunc, PPPMCopyComplexFunc) : loop index i < total - skip_n
  ONIKA_HOST_DEVICE_FUNC inline size_t pppm_skip_own( size_t i, size_t skip_lo, size_t skip_n ) { return i < skip_lo ? i : i + skip_n; }

  struct PPPMPackRealFunc
  {
    const double * __restrict__ src = nullptr;
    const size_t * __restrict__ idx = nullptr;
    double * __restrict__ buf = nullptr;
    size_t skip_lo = 0, skip_n = 0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { const size_t j = pppm_skip_own( i , skip_lo , skip_n ); buf[j] = src[ idx[j] ]; }
  };

  // several sources may map to the same destination point (several bricks) : atomic adds
  struct PPPMUnpackAddRealFunc
  {
    double * __restrict__ dst = nullptr;
    const size_t * __restrict__ idx = nullptr;
    const double * __restrict__ buf = nullptr;
    size_t skip_lo = 0, skip_n = 0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { const size_t j = pppm_skip_own( i , skip_lo , skip_n ); ONIKA_CU_ATOMIC_ADD( dst[ idx[j] ] , buf[j] ); }
  };

  // own points, from the send map to the receive map (both list them in the same order)
  struct PPPMCopyAddRealFunc
  {
    const double * __restrict__ src = nullptr;
    const size_t * __restrict__ sidx = nullptr;
    double * __restrict__ dst = nullptr;
    const size_t * __restrict__ didx = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { ONIKA_CU_ATOMIC_ADD( dst[ didx[i] ] , src[ sidx[i] ] ); }
  };

  static constexpr int PPPM_MAX_PACK_MESHES = 3;

  // K complex meshes packed point by point : buf[i*K+k] = src[k][idx[i]]
  struct PPPMPackComplexFunc
  {
    const Complexd * src[PPPM_MAX_PACK_MESHES] = { nullptr , nullptr , nullptr };
    int K = 1;
    const size_t * __restrict__ idx = nullptr;
    Complexd * __restrict__ buf = nullptr;
    size_t skip_lo = 0, skip_n = 0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      const size_t ii = pppm_skip_own( i , skip_lo , skip_n );
      const size_t j = idx[ii];
      for( int k = 0 ; k < K ; k++ ) buf[ ii*K + k ] = src[k][j];
    }
  };

  struct PPPMUnpackComplexFunc
  {
    Complexd * dst[PPPM_MAX_PACK_MESHES] = { nullptr , nullptr , nullptr };
    int K = 1;
    const size_t * __restrict__ idx = nullptr;
    const Complexd * __restrict__ buf = nullptr;
    size_t skip_lo = 0, skip_n = 0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      const size_t ii = pppm_skip_own( i , skip_lo , skip_n );
      const size_t j = idx[ii];
      for( int k = 0 ; k < K ; k++ ) dst[k][j] = buf[ ii*K + k ];
    }
  };

  struct PPPMCopyComplexFunc
  {
    const Complexd * src[PPPM_MAX_PACK_MESHES] = { nullptr , nullptr , nullptr };
    Complexd * dst[PPPM_MAX_PACK_MESHES] = { nullptr , nullptr , nullptr };
    int K = 1;
    const size_t * __restrict__ sidx = nullptr;
    const size_t * __restrict__ didx = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      const size_t a = sidx[i], b = didx[i];
      for( int k = 0 ; k < K ; k++ ) dst[k][b] = src[k][a];
    }
  };

}
}

namespace exanb
{
  template<bool PerAtomCharge> struct ComputeCellParticlesTraits< exaStamp::PPPMSpreadFunc<PerAtomCharge> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<> struct ComputeCellParticlesTraits< exaStamp::PPPMBrickBoundsFunc >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<bool PerAtomCharge> struct ComputeCellParticlesTraits< exaStamp::PPPMSlabDipoleFunc<PerAtomCharge> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<bool PerAtomCharge, bool ComputeEnergy, bool ComputeVirial>
  struct ComputeCellParticlesTraits< exaStamp::PPPMForceFunc<PerAtomCharge,ComputeEnergy,ComputeVirial> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace onika
{
  namespace parallel
  {
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMZeroFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMSumThreadMeshesFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMLoadDensityFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMApplyGreenFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMPairFactorFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMPackRealFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMUnpackAddRealFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMPackComplexFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMUnpackComplexFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMCopyAddRealFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMCopyComplexFunc > { static inline constexpr bool CudaCompatible = true; };
  }
}

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class PPPMLongRangePC : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( PPPMParameters  , pppm_config          , INPUT , OPTIONAL );
    ADD_SLOT( GridT           , grid                 , INPUT_OUTPUT );
    ADD_SLOT( Domain          , domain               , INPUT , REQUIRED );
    ADD_SLOT( double          , rcut_max             , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species              , INPUT , REQUIRED );
    ADD_SLOT( bool            , per_atom_charge      , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( MPI_Comm        , mpi                  , INPUT );
    ADD_SLOT( bool            , trigger_thermo_state , INPUT , OPTIONAL );

    // mesh work buffers, kept between time steps (replicated mesh)
    onika::memory::CudaMMVector<double> m_density;
    onika::memory::CudaMMVector<double> m_thread_density; // CPU, several OpenMP threads : one mesh (or brick) per thread
    onika::memory::CudaMMVector<Complexd> m_work1; // rho(k), then V(k)
    onika::memory::CudaMMVector<Complexd> m_field; // 2 complex meshes : (vdx,vdy) , (vdz,u)
    onika::memory::CudaMMVector<Complexd> m_vir;   // 3 complex meshes : (xx,yy) , (zz,xy) , (xz,yz)
    onika::memory::CudaMMVector<double> m_slab_sum; // slab correction : sum q.z , sum q.z^2
    PPPMFFT m_fft;
    bool m_gpu = false;

    // distributed mesh (see PPPMDecomposition) : local brick of the particles' stencils, z slab of the real space mesh,
    // local rows of the reciprocal space mesh, and index maps of the exchanges between them
    static constexpr int BRICK_PAD = 1; // brick margin (mesh points) : maps are rebuilt only when particles leave it
    PPPMDistFFT m_dfft;
    onika::memory::CudaMMVector<int> m_bounds;
    std::vector<int> m_all_bricks;      // lo[3],dims[3] of every rank's brick, as the brick maps were built for
    int m_brick_lo[3] = {0,0,0}, m_brick_dims[3] = {0,0,0};
    size_t m_brick_n = 0;
    onika::memory::CudaMMVector<double> m_brick_density;
    onika::memory::CudaMMVector<Complexd> m_brick_mesh;  // 5 bricks : (vdx,vdy) or u , (vdz,u) , virial pairs
    onika::memory::CudaMMVector<double> m_slab_density;
    onika::memory::CudaMMVector<Complexd> m_slab_mesh;   // 3 slabs : FFT work, then backward outputs
    onika::memory::CudaMMVector<Complexd> m_kwork;       // rho(k), then V(k), local rows
    onika::memory::CudaMMVector<Complexd> m_kout;        // 3 local reciprocal meshes, backward inputs
    onika::memory::CudaMMVector<size_t> m_t_slab_idx, m_t_k_idx;     // slab <-> rows transpose, grouped by peer
    std::vector<size_t> m_t_slab_cnt, m_t_k_cnt;
    onika::memory::CudaMMVector<size_t> m_b_brick_idx, m_b_slab_idx; // brick <-> slab exchanges, grouped by peer
    std::vector<size_t> m_b_brick_cnt, m_b_slab_cnt;
    int m_tmaps_key[5] = { -1, -1, -1, -1, -1 }; // nx,ny,nz,nprocs,rank the transpose maps were built for
    onika::memory::CudaMMVector<double> m_sendbuf, m_recvbuf;

    // per mesh point loop : GPU (one thread per point) when the FFT runs there, OpenMP otherwise
    template<class FuncT>
    inline void mesh_for( size_t n, const FuncT& func )
    {
      if( n == 0 ) return;
      onika::parallel::ParallelForOptions opts;
      opts.enable_gpu = m_gpu;
      onika::parallel::parallel_for( n , func , parallel_execution_context() , opts );
    }

    inline PPPMPairFactorFunc pair_factor( const PPPMParameters& p, PPPMBackwardKind kind, bool with_u, const Complexd* w1, Complexd* const out[] ) const
    {
      const int nout = ( kind == PPPMBackwardKind::POTENTIAL ) ? 1 : ( ( kind == PPPMBackwardKind::FIELD ) ? 2 : 3 );
      PPPMPairFactorFunc func = {};
      func.w1 = w1;
      for( int q = 0 ; q < nout ; q++ ) func.w2[q] = out[q];
      func.fkx = p.fkx.data(); func.fky = p.fky.data(); func.fkz = p.fkz.data();
      func.g_ewald = p.g_ewald;
      func.kind = kind;
      func.with_u = with_u;
      func.nx = p.nx; func.ny = p.ny; func.nz = p.nz;
      func.m_rows = ! m_gpu;
      if( p.dec.distributed ) func.krow_partner = p.krow_partner.data();
      return func;
    }

    static inline int backward_count( PPPMBackwardKind kind ) { return ( kind == PPPMBackwardKind::POTENTIAL ) ? 1 : ( ( kind == PPPMBackwardKind::FIELD ) ? 2 : 3 ); }

    // backward FFTs of the FIELD (2 meshes out), POTENTIAL (1) or VIRIAL (3) pairs, see PPPMPairFactorFunc
    inline void backward_pairs( const PPPMParameters& p, PPPMBackwardKind kind, bool with_u, Complexd* const out[] )
    {
      const int nout = backward_count( kind );
      const PPPMPairFactorFunc func = pair_factor( p , kind , with_u , m_work1.data() , out );
      mesh_for( func.m_rows ? size_t(p.ny) * p.nz : p.nfft() , func );
      for( int q = 0 ; q < nout ; q++ ) m_fft.backward( out[q] );
      m_fft.sync();
    }

    // ------------------------------- distributed mesh -------------------------------

    // MPI_Alltoallv of points made of dpp doubles each, counts given in points per peer. This rank's own points are
    // not exchanged (copied directly, see move_real_add / move_complex) : their place in the buffers is left unused.
    inline void exchange( const double* send, const std::vector<size_t>& scnt, double* recv, const std::vector<size_t>& rcnt, int dpp, int rank )
    {
      const int P = scnt.size();
      if( P == 1 ) return;
      std::vector<int> sc(P), sd(P), rc(P), rd(P);
      size_t so = 0, ro = 0;
      for( int r = 0 ; r < P ; r++ )
      {
        const size_t s = scnt[r] * dpp, q = rcnt[r] * dpp;
        if( s > size_t(INT_MAX) || q > size_t(INT_MAX) || so > size_t(INT_MAX) || ro > size_t(INT_MAX) ) fatal_error() << "coulombic_pppm : MPI message too large" << std::endl;
        sc[r] = ( r == rank ) ? 0 : s; sd[r] = so; so += s;
        rc[r] = ( r == rank ) ? 0 : q; rd[r] = ro; ro += q;
      }
      MPI_Alltoallv( send , sc.data() , sd.data() , MPI_DOUBLE , recv , rc.data() , rd.data() , MPI_DOUBLE , *mpi );
    }

    static inline size_t total( const std::vector<size_t>& c ) { size_t t = 0; for( auto x : c ) t += x; return t; }
    static inline size_t offset( const std::vector<size_t>& c, int rank ) { size_t t = 0; for( int r = 0 ; r < rank ; r++ ) t += c[r]; return t; }

    // real mesh points sidx (grouped by destination rank, counts scnt) of src added to points didx (grouped by source
    // rank, counts rcnt) of dst. Own points are added directly ; the send and receive maps list them in the same order.
    inline void move_real_add( const double* src, const size_t* sidx, const std::vector<size_t>& scnt,
                               double* dst, const size_t* didx, const std::vector<size_t>& rcnt, int rank )
    {
      const size_t ns = total( scnt ), nr = total( rcnt ), own = scnt[rank];
      const size_t so = offset( scnt , rank ), ro = offset( rcnt , rank );
      mesh_for( own , PPPMCopyAddRealFunc{ src , sidx + so , dst , didx + ro } );
      if( ns == own && nr == own ) return;
      resize_buffers( ns , nr );
      mesh_for( ns - own , PPPMPackRealFunc{ src , sidx , m_sendbuf.data() , so , own } );
      exchange( m_sendbuf.data() , scnt , m_recvbuf.data() , rcnt , 1 , rank );
      mesh_for( nr - own , PPPMUnpackAddRealFunc{ dst , didx , m_recvbuf.data() , ro , own } );
    }

    // same for K complex meshes, destination points overwritten
    inline void move_complex( int K, const Complexd* const src[], const size_t* sidx, const std::vector<size_t>& scnt,
                              Complexd* const dst[], const size_t* didx, const std::vector<size_t>& rcnt, int rank )
    {
      const size_t ns = total( scnt ), nr = total( rcnt ), own = scnt[rank];
      const size_t so = offset( scnt , rank ), ro = offset( rcnt , rank );
      PPPMCopyComplexFunc copy = {}; copy.K = K; copy.sidx = sidx + so; copy.didx = didx + ro;
      for( int q = 0 ; q < K ; q++ ) { copy.src[q] = src[q]; copy.dst[q] = dst[q]; }
      mesh_for( own , copy );
      if( ns == own && nr == own ) return;
      resize_buffers( 2*K*ns , 2*K*nr );
      PPPMPackComplexFunc pack = {}; pack.K = K; pack.idx = sidx; pack.buf = reinterpret_cast<Complexd*>( m_sendbuf.data() ); pack.skip_lo = so; pack.skip_n = own;
      for( int q = 0 ; q < K ; q++ ) pack.src[q] = src[q];
      mesh_for( ns - own , pack );
      exchange( m_sendbuf.data() , scnt , m_recvbuf.data() , rcnt , 2*K , rank );
      PPPMUnpackComplexFunc unpack = {}; unpack.K = K; unpack.idx = didx; unpack.buf = reinterpret_cast<const Complexd*>( m_recvbuf.data() ); unpack.skip_lo = ro; unpack.skip_n = own;
      for( int q = 0 ; q < K ; q++ ) unpack.dst[q] = dst[q];
      mesh_for( nr - own , unpack );
    }

    inline void resize_buffers( size_t send_doubles, size_t recv_doubles )
    {
      if( m_sendbuf.size() < send_doubles ) m_sendbuf.resize( send_doubles );
      if( m_recvbuf.size() < recv_doubles ) m_recvbuf.resize( recv_doubles );
    }

    // z slab (index ((z-z0)*ny+y)*nx+x) <-> local rows (index (l*nx+x)*nz+z) : points this rank sends to / receives from each peer
    inline void build_transpose_maps( const PPPMParameters& p )
    {
      const auto& d = p.dec;
      const int key[5] = { p.nx , p.ny , p.nz , d.nprocs , d.rank };
      if( std::equal( key , key+5 , m_tmaps_key ) ) return;
      std::copy( key , key+5 , m_tmaps_key );
      const int nx = p.nx, ny = p.ny, nz = p.nz, P = d.nprocs;
      const int z0 = d.z0(), nzl = d.nzl(), nyl = d.nyl();
      m_t_slab_cnt.assign( P , 0 ); m_t_k_cnt.assign( P , 0 );
      std::vector<size_t> sidx, kidx;
      for( int s = 0 ; s < P ; s++ )
      {
        const int* rows = d.rows(s);
        for( int i = 0 ; i < d.nyl(s) ; i++ )
          for( int z = 0 ; z < nzl ; z++ )
            for( int x = 0 ; x < nx ; x++ ) sidx.push_back( ( size_t(z) * ny + rows[i] ) * nx + x );
        m_t_slab_cnt[s] = size_t( d.nyl(s) ) * nzl * nx;
      }
      for( int t = 0 ; t < P ; t++ )
      {
        for( int l = 0 ; l < nyl ; l++ )
          for( int z = d.zlo[t] ; z < d.zlo[t+1] ; z++ )
            for( int x = 0 ; x < nx ; x++ ) kidx.push_back( ( size_t(l) * nx + x ) * nz + z );
        m_t_k_cnt[t] = size_t( nyl ) * d.nzl(t) * nx;
      }
      m_t_slab_idx.assign( sidx.begin() , sidx.end() );
      m_t_k_idx.assign( kidx.begin() , kidx.end() );
      ldbg << "coulombic_pppm : transpose maps, slab "<< z0 <<"+"<< nzl <<" planes, "<< nyl <<" rows" << std::endl;
    }

    // local brick (unwrapped indices) <-> z slabs of the owners of its planes
    inline void build_brick_maps( const PPPMParameters& p )
    {
      const auto& d = p.dec;
      const int nx = p.nx, ny = p.ny, nz = p.nz, P = d.nprocs;
      const int z0 = d.z0(), z1 = d.z0() + d.nzl();
      const int* b = m_all_bricks.data() + 6*d.rank;
      m_b_brick_cnt.assign( P , 0 ); m_b_slab_cnt.assign( P , 0 );
      std::vector<size_t> bidx, sidx;
      for( int s = 0 ; s < P ; s++ )
      {
        const size_t before = bidx.size();
        for( int bz = 0 ; bz < b[5] ; bz++ )
        {
          if( d.z_owner( pppm_wrap( b[2] + bz , nz ) ) != s ) continue;
          for( int by = 0 ; by < b[4] ; by++ )
            for( int bx = 0 ; bx < b[3] ; bx++ ) bidx.push_back( ( size_t(bz) * b[4] + by ) * b[3] + bx );
        }
        m_b_brick_cnt[s] = bidx.size() - before;
      }
      for( int t = 0 ; t < P ; t++ )
      {
        const int* bt = m_all_bricks.data() + 6*t;
        const size_t before = sidx.size();
        for( int bz = 0 ; bz < bt[5] ; bz++ )
        {
          const int wz = pppm_wrap( bt[2] + bz , nz );
          if( wz < z0 || wz >= z1 ) continue;
          for( int by = 0 ; by < bt[4] ; by++ )
          {
            const int wy = pppm_wrap( bt[1] + by , ny );
            for( int bx = 0 ; bx < bt[3] ; bx++ ) sidx.push_back( ( size_t(wz - z0) * ny + wy ) * nx + pppm_wrap( bt[0] + bx , nx ) );
          }
        }
        m_b_slab_cnt[t] = sidx.size() - before;
      }
      m_b_brick_idx.assign( bidx.begin() , bidx.end() );
      m_b_slab_idx.assign( sidx.begin() , sidx.end() );
      ldbg << "coulombic_pppm : brick maps, brick "<<b[0]<<","<<b[1]<<","<<b[2]<<" + "<<b[3]<<"x"<<b[4]<<"x"<<b[5]
           <<", "<< bidx.size() <<" points sent, "<< sidx.size() <<" received" << std::endl;
    }

    // backward chain of the distributed mesh : pair factor on local rows, z FFTs, transpose to slabs, xy FFTs, slab ->
    // bricks. Outputs nout complex bricks.
    inline void backward_pairs_distributed( const PPPMParameters& p, PPPMBackwardKind kind, bool with_u, Complexd* const brick_out[] )
    {
      const auto& d = p.dec;
      const int nout = backward_count( kind );
      const size_t nk = p.nk_local();
      const size_t nslab = size_t( d.nzl() ) * p.ny * p.nx;
      Complexd* kout[3] = { m_kout.data() , m_kout.data() + nk , m_kout.data() + 2*nk };
      Complexd* slab_out[3] = { m_slab_mesh.data() , m_slab_mesh.data() + nslab , m_slab_mesh.data() + 2*nslab };

      const PPPMPairFactorFunc func = pair_factor( p , kind , with_u , m_kwork.data() , kout );
      mesh_for( func.m_rows ? size_t( d.nyl() ) * p.nx : nk , func );
      for( int q = 0 ; q < nout ; q++ ) m_dfft.columns( kout[q] , false );
      m_dfft.sync();

      // rows -> slabs
      move_complex( nout , kout , m_t_k_idx.data() , m_t_k_cnt , slab_out , m_t_slab_idx.data() , m_t_slab_cnt , d.rank );
      for( int q = 0 ; q < nout ; q++ ) m_dfft.planes( slab_out[q] , false );
      m_dfft.sync();

      // slabs -> bricks
      move_complex( nout , slab_out , m_b_slab_idx.data() , m_b_slab_cnt , brick_out , m_b_brick_idx.data() , m_b_brick_cnt , d.rank );
    }

  public:
    inline void execute () override final
    {
      const bool log_energy = trigger_thermo_state.has_value() ? *trigger_thermo_state : false;

      if( ! pppm_config.has_value() )
      {
        ldbg << "pppm_config not set, skip coulombic_pppm" << std::endl;
        return;
      }
      const auto & p = *pppm_config;
      *rcut_max = std::max( *rcut_max , p.radius );

      if( grid->number_of_cells() == 0 ) return;

      const Mat3d cell = ewald_cell_matrix( domain->xform() , domain->bounds_size() );
      if( p.volume == 0.0 || ! ewald_same_cell( cell , p.cell ) )
      {
        fatal_error() << "coulombic_pppm : domain cell "<<cell<<" differs from the one PPPM was set up for "
                      << p.cell << ". Call coulombic_pppm_init before force computation when the cell changes." << std::endl;
      }

      int nprocs = 1, rank = 0;
      MPI_Comm_size( *mpi , &nprocs );
      MPI_Comm_rank( *mpi , &rank );
      const bool distributed = p.dec.distributed;
      if( distributed && ( p.dec.nprocs != nprocs || p.dec.rank != rank ) )
      {
        fatal_error() << "coulombic_pppm : mesh decomposition was built for "<<p.dec.nprocs<<" ranks, running on "<<nprocs << std::endl;
      }

      // GPU path when a device is available : particle kernels, mesh loops and cuFFT all run there (unified memory)
      const bool gpu_available = ( global_cuda_ctx() != nullptr ) && global_cuda_ctx()->has_devices() && PPPMFFT::gpu_support();
      void* stream = nullptr;
#     ifdef EXASTAMP_PPPM_CUFFT
      if( gpu_available ) stream = global_cuda_ctx()->getThreadStream(0);
#     endif
      m_gpu = gpu_available;

      const size_t nfft = p.nfft();
      const size_t nk = p.nk_local();
      const size_t nslab = distributed ? size_t( p.dec.nzl() ) * p.ny * p.nx : 0;
      if( ! distributed )
      {
        if( m_density.size() != nfft )
        {
          m_density.resize( nfft );
          m_work1.resize( nfft );
          m_field.resize( 2*nfft );
          m_vir.clear();
        }
        m_fft.resize( p.nx , p.ny , p.nz , gpu_available , stream );
        if( log_energy && m_vir.size() != 3*nfft ) m_vir.resize( 3*nfft );
      }
      else
      {
        if( m_slab_density.size() != nslab ) { m_slab_density.resize( nslab ); m_slab_mesh.resize( 3*nslab ); }
        if( m_kwork.size() != nk ) { m_kwork.resize( nk ); m_kout.resize( 3*nk ); }
        m_dfft.resize( p.nx , p.ny , p.nz , p.dec.nzl() , p.dec.nyl() * p.nx , gpu_available , stream );
        build_transpose_maps( p );
      }

      ReadOnlyPPPMParameters ro( p , domain->bounds().bmin , domain->bounds_size() );

      auto rx = grid->field_accessor( field::rx );
      auto ry = grid->field_accessor( field::ry );
      auto rz = grid->field_accessor( field::rz );
      auto fx = grid->field_accessor( field::fx );
      auto fy = grid->field_accessor( field::fy );
      auto fz = grid->field_accessor( field::fz );
      auto ep = grid->field_accessor( field::ep );
      auto virial = grid->field_accessor( field::virial );

      auto compute_with_charges = [&]( auto per_atom_charge_tag , auto charge_or_type )
      {
        static constexpr bool PerAtomCharge = decltype(per_atom_charge_tag)::value;

        // 0. slab correction : total dipole along z (LAMMPS slabcorr), needed by the force pass
        if( ro.slab )
        {
          if( m_slab_sum.size() != 2 ) m_slab_sum.resize( 2 );
          mesh_for( 2 , PPPMZeroFunc{ m_slab_sum.data() } );
          PPPMSlabDipoleFunc<PerAtomCharge> dipole_func = { species->data() , m_slab_sum.data() , ro.zscale };
          compute_cell_particles( *grid , false , dipole_func , onika::make_flat_tuple(rz,charge_or_type) , parallel_execution_context() );
          double sums[2] = { m_slab_sum[0] , m_slab_sum[1] };
          if( nprocs > 1 ) MPI_Allreduce( MPI_IN_PLACE , sums , 2 , MPI_DOUBLE , MPI_SUM , *mpi );
          ro.dipole = sums[0];
          ro.dipole_r2 = sums[1];
        }

        // distributed : brick of the local particles' stencils (padded, kept while particles stay inside), shared with all ranks
        if( distributed )
        {
          if( m_bounds.size() != 6 ) m_bounds.resize( 6 );
          for( int i = 0 ; i < 3 ; i++ ) { m_bounds[i] = std::numeric_limits<int>::max(); m_bounds[3+i] = std::numeric_limits<int>::min(); }
          PPPMBrickBoundsFunc bounds_func = { ro , m_bounds.data() };
          compute_cell_particles( *grid , false , bounds_func , onika::make_flat_tuple(rx,ry,rz,charge_or_type) , parallel_execution_context() );
          int brick[6] = { 0, 0, 0, 0, 0, 0 };
          if( m_bounds[0] <= m_bounds[3] )
          {
            // a direction where the padded stencils span the whole mesh is the whole periodic mesh (lo 0, wrapped
            // indices) : no mesh point is held twice
            const int n[3] = { int(p.nx) , int(p.ny) , int(p.nz) };
            bool inside = m_brick_n > 0;
            for( int i = 0 ; i < 3 ; i++ )
              if( m_brick_dims[i] != n[i] ) inside = inside && m_bounds[i] >= m_brick_lo[i] && m_bounds[3+i] < m_brick_lo[i] + m_brick_dims[i];
            for( int i = 0 ; i < 3 ; i++ )
            {
              const int len = m_bounds[3+i] - m_bounds[i] + 1 + 2*BRICK_PAD;
              if( inside )         { brick[i] = m_brick_lo[i]; brick[3+i] = m_brick_dims[i]; }
              else if( len >= n[i] ) { brick[i] = 0; brick[3+i] = n[i]; }
              else                 { brick[i] = m_bounds[i] - BRICK_PAD; brick[3+i] = len; }
            }
          }
          std::vector<int> all( 6*nprocs );
          MPI_Allgather( brick , 6 , MPI_INT , all.data() , 6 , MPI_INT , *mpi );
          for( int i = 0 ; i < 3 ; i++ ) { m_brick_lo[i] = brick[i]; m_brick_dims[i] = brick[3+i]; }
          m_brick_n = size_t( brick[3] ) * brick[4] * brick[5];
          if( all != m_all_bricks )
          {
            m_all_bricks = all;
            build_brick_maps( p );
          }
          if( m_brick_density.size() != m_brick_n ) { m_brick_density.resize( m_brick_n ); m_brick_mesh.resize( 5*m_brick_n ); }
          ro.set_brick( m_brick_lo , m_brick_dims );
        }

        // 1. charge density of local particles : replicated, whole mesh summed over all ranks ; distributed, local brick
        const size_t nmesh = distributed ? m_brick_n : nfft;
        double * __restrict__ density = distributed ? m_brick_density.data() : m_density.data();
        const int nthreads = m_gpu ? 0 : omp_get_max_threads();
        double * thread_density = nullptr;
        if( nthreads > 1 )
        {
          if( m_thread_density.size() != nthreads * nmesh ) m_thread_density.resize( nthreads * nmesh );
          thread_density = m_thread_density.data();
          mesh_for( nthreads * nmesh , PPPMZeroFunc{ thread_density } );
        }
        else mesh_for( nmesh , PPPMZeroFunc{ density } );
        PPPMSpreadFunc<PerAtomCharge> spread_func = { ro , species->data() , density , thread_density , nthreads , nmesh };
        compute_cell_particles( *grid , false , spread_func , onika::make_flat_tuple(rx,ry,rz,charge_or_type) , parallel_execution_context() );
        if( nthreads > 1 ) mesh_for( nmesh , PPPMSumThreadMeshesFunc{ thread_density , density , nthreads , nmesh } );

        PPPMMeshes meshes = {};
        if( ! distributed )
        {
          if( nprocs > 1 ) MPI_Allreduce( MPI_IN_PLACE , density , nfft , MPI_DOUBLE , MPI_SUM , *mpi );

          // 2. rho(k), then V(k) = G(k) rho(k) / N (LAMMPS poisson_ik)
          mesh_for( nfft , PPPMLoadDensityFunc{ density , m_work1.data() } );
          m_fft.forward( m_work1.data() );
          m_fft.sync();
          mesh_for( nfft , PPPMApplyGreenFunc{ m_work1.data() , p.greensfn.data() , 1.0 / double(nfft) } );

          // 3. gradient of the potential, i.k V(k) back to real space (LAMMPS poisson_ik), and on energy steps the
          // potential and virial meshes (LAMMPS poisson_peratom), two real meshes per backward FFT
          // ad (LAMMPS poisson_ad) : the potential only, its gradient is taken on the particles
          Complexd * exy = m_field.data();
          Complexd * ezu = m_field.data() + nfft;
          if( p.diff_ad )
          {
            Complexd* const out[1] = { exy };
            backward_pairs( p , PPPMBackwardKind::POTENTIAL , false , out );
          }
          else
          {
            Complexd* const out[2] = { exy , ezu };
            backward_pairs( p , PPPMBackwardKind::FIELD , log_energy , out );
          }
          meshes = { exy , ezu , nullptr , nullptr , nullptr };
          if( log_energy )
          {
            Complexd* const out[3] = { m_vir.data() , m_vir.data() + nfft , m_vir.data() + 2*nfft };
            meshes.v01 = out[0];
            meshes.v23 = out[1];
            meshes.v45 = out[2];
            backward_pairs( p , PPPMBackwardKind::VIRIAL , false , out );
          }
        }
        else
        {
          // 2. bricks -> z slabs (sum), rho(k) : xy FFTs, transpose to rows, z FFTs, V(k) = G(k) rho(k) / N
          mesh_for( nslab , PPPMZeroFunc{ m_slab_density.data() } );
          move_real_add( density , m_b_brick_idx.data() , m_b_brick_cnt , m_slab_density.data() , m_b_slab_idx.data() , m_b_slab_cnt , rank );
          Complexd* slab_work = m_slab_mesh.data();
          mesh_for( nslab , PPPMLoadDensityFunc{ m_slab_density.data() , slab_work } );
          m_dfft.planes( slab_work , true );
          m_dfft.sync();
          {
            const Complexd* const src[1] = { slab_work };
            Complexd* const dst[1] = { m_kwork.data() };
            move_complex( 1 , src , m_t_slab_idx.data() , m_t_slab_cnt , dst , m_t_k_idx.data() , m_t_k_cnt , rank );
          }
          m_dfft.columns( m_kwork.data() , true );
          m_dfft.sync();
          mesh_for( nk , PPPMApplyGreenFunc{ m_kwork.data() , p.greensfn.data() , 1.0 / double(nfft) } );

          // 3. backward chains down to the bricks (same pairs as the replicated mesh)
          Complexd* b[5];
          for( int q = 0 ; q < 5 ; q++ ) b[q] = m_brick_mesh.data() + q*m_brick_n;
          if( p.diff_ad )
          {
            Complexd* const out[1] = { b[0] };
            backward_pairs_distributed( p , PPPMBackwardKind::POTENTIAL , false , out );
          }
          else
          {
            Complexd* const out[2] = { b[0] , b[1] };
            backward_pairs_distributed( p , PPPMBackwardKind::FIELD , log_energy , out );
          }
          meshes = { b[0] , b[1] , nullptr , nullptr , nullptr };
          if( log_energy )
          {
            Complexd* const out[3] = { b[2] , b[3] , b[4] };
            meshes.v01 = b[2]; meshes.v23 = b[3]; meshes.v45 = b[4];
            backward_pairs_distributed( p , PPPMBackwardKind::VIRIAL , false , out );
          }
        }

        // 4. interpolate field (and energy, virial) back to local particles
        if( log_energy )
        {
          PPPMForceFunc<PerAtomCharge,true,true> force_func = { ro , species->data() , meshes };
          compute_cell_particles( *grid , false , force_func , onika::make_flat_tuple(fx,fy,fz,ep,virial,rx,ry,rz,charge_or_type) , parallel_execution_context() );
        }
        else
        {
          PPPMForceFunc<PerAtomCharge,false,false> force_func = { ro , species->data() , meshes };
          compute_cell_particles( *grid , false , force_func , onika::make_flat_tuple(fx,fy,fz,rx,ry,rz,charge_or_type) , parallel_execution_context() );
        }
      };

      if( *per_atom_charge ) compute_with_charges( std::true_type{}  , grid->field_accessor( field::charge ) );
      else                   compute_with_charges( std::false_type{} , grid->field_accessor( field::type ) );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Reciprocal space part of PPPM long range coulomb (ik or ad differentiation), set up by coulombic_pppm_init.
Same algorithm as LAMMPS kspace_style pppm, orthogonal and triclinic cells (ad : orthogonal only, as LAMMPS), with
optional slab correction (EW3DC, z non periodic). Use with coulombic_ewald_short_range for
the real space part. Computes forces, and when trigger_thermo_state is true, per particle energy (reciprocal + self +
neutralizing background) and per particle reciprocal virial.
Mesh decomposition (coulombic_pppm_init mesh_decomposition) : distributed (default), each rank spreads its particles on
a local brick, the mesh is split in z slabs for the xy FFTs and in y rows for the z FFTs (MPI_Alltoallv exchanges, the
rank's own points copied directly) ; replicated, the whole mesh on every MPI rank (density summed with MPI_Allreduce,
every rank runs the full FFTs). When a GPU
is available, particle kernels, mesh loops and FFTs (cuFFT) run on the GPU ; otherwise on CPU with pocketfft.
)EOF";
    }
  };

  template<class GridT> using PPPMLongRangePCTmpl = PPPMLongRangePC<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(coulombic_pppm)
  {
    OperatorNodeFactory::instance()->register_factory( "coulombic_pppm" , make_grid_variant_operator<PPPMLongRangePCTmpl> );
  }

}
}
