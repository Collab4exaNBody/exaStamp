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

// PPPM (particle-particle particle-mesh) long range coulomb, ik differentiation, orthogonal and triclinic cells.
// Follows LAMMPS KSPACE/pppm.cpp (kspace_style pppm) : same grid selection, g_ewald adjustment, charge assignment
// and optimal (Hockney-Eastwood) influence function, so results match LAMMPS up to round-off.
// All mesh quantities are expressed in fractional (lamda) coordinates s = (r - bmin) / bounds_size, which is what
// LAMMPS does for triclinic cells ; for orthogonal cells it is the same algorithm.

#include <exaStamp/potential/coulombic/ewald.h>
#include <onika/memory/allocator.h>
#include <onika/cuda/cuda.h>
#include <onika/log.h>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <vector>
#include <mpi.h>

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using namespace exanb;

  namespace pppm_constants
  {
    static constexpr int MAXORDER = 7;
    static constexpr int OFFSET = 16384;       // avoids int(-0.75)=0 when mapping particles to the mesh
    static constexpr double EPS_HOC = 1.0e-7;  // aliasing sum cutoff of the influence function
    static constexpr double SMALL = 0.00001;   // g_ewald Newton-Raphson tolerance (eV/ang, absolute as in LAMMPS)
    static constexpr int LARGE = 10000;
    static constexpr double FOUR_PI_LMP = 12.5663706; // truncated 4.pi used by LAMMPS compute_gf_ik, kept for parity
  }

  // Distributed mesh (MPI ranks > 1) : real space mesh split in z slabs (rank r owns planes [zlo[r],zlo[r+1]) ), reciprocal
  // space split in y rows, each rank's rows closed under y -> -y so that the -k partner of a local k point is local.
  // Local reciprocal layout : index (l*nx + ix)*nz + iz for local row l (global y = rows_all[row_first[rank]+l]).
  struct PPPMDecomposition
  {
    bool distributed = false;
    int nprocs = 1, rank = 0;
    std::vector<int> zlo;       // nprocs+1
    std::vector<int> row_first; // nprocs+1
    std::vector<int> rows_all;  // rows grouped by owner rank
    std::vector<int> row_owner; // ny
    std::vector<int> row_local; // ny : index of a row in its owner's list

    inline int z0() const { return zlo[rank]; }
    inline int nzl() const { return zlo[rank+1] - zlo[rank]; }
    inline int nzl( int r ) const { return zlo[r+1] - zlo[r]; }
    inline int nyl() const { return row_first[rank+1] - row_first[rank]; }
    inline int nyl( int r ) const { return row_first[r+1] - row_first[r]; }
    inline const int* rows() const { return rows_all.data() + row_first[rank]; }
    inline const int* rows( int r ) const { return rows_all.data() + row_first[r]; }
    inline int z_owner( int z ) const { int r = 0; while( z >= zlo[r+1] ) ++r; return r; }
  };

  inline void pppm_make_decomposition( PPPMDecomposition& d, bool distributed, int nprocs, int rank, int nx, int ny, int nz )
  {
    d = PPPMDecomposition{};
    d.distributed = distributed;
    d.nprocs = distributed ? nprocs : 1;
    d.rank = distributed ? rank : 0;
    const int P = d.nprocs;
    d.zlo.resize( P+1 );
    for( int r = 0 ; r <= P ; r++ ) d.zlo[r] = int( ( int64_t(r) * nz ) / P );
    // row groups {0}, {y, ny-y}, {ny/2} if ny even : group g holds rows g and ny-g
    const int ngroups = ny/2 + 1;
    d.row_first.resize( P+1 );
    d.row_owner.assign( ny , 0 );
    d.row_local.assign( ny , 0 );
    d.rows_all.clear();
    for( int r = 0 ; r < P ; r++ )
    {
      d.row_first[r] = d.rows_all.size();
      const int g0 = int( ( int64_t(r) * ngroups ) / P ), g1 = int( ( int64_t(r+1) * ngroups ) / P );
      for( int g = g0 ; g < g1 ; g++ )
      {
        const int ys[2] = { g , ( ny - g ) % ny };
        const int n = ( ys[1] == ys[0] ) ? 1 : 2;
        for( int i = 0 ; i < n ; i++ )
        {
          d.row_owner[ ys[i] ] = r;
          d.row_local[ ys[i] ] = d.rows_all.size() - d.row_first[r];
          d.rows_all.push_back( ys[i] );
        }
      }
    }
    d.row_first[P] = d.rows_all.size();
  }

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PPPMParameters
  {
    // user settings
    double accuracy_relative = 0.0;
    double g_ewald_user = 0.0;  // 0 = automatic
    double radius = 0.0;
    long order = 5;
    long mesh_user[3] = { 0, 0, 0 }; // 0 = automatic
    bool diff_ad = false;            // LAMMPS kspace_modify diff ad (analytic differentiation), otherwise diff ik
    double slab_user = 0.0;          // LAMMPS kspace_modify slab <volfactor> (EW3DC), 0 = no slab correction
    bool slab_auto = false;          // LAMMPS kspace_modify slab auto : volfactor computed from accuracy and g_ewald

    // derived
    double g_ewald = 0.0;
    int nx = 0, ny = 0, nz = 0;
    Mat3d cell = { 0.,0.,0., 0.,0.,0., 0.,0.,0. }; // cell matrix H the influence function was built for
    double volume = 0.0;             // volume of the (z extended with slab) cell
    double slab_volfactor = 1.0;     // z extension of the cell for the slab correction (1 = no slab)
    double qsum = 0.0;
    double qsqsum = 0.0;
    double estimated_accuracy = 0.0; // absolute rms force error estimate, eV/ang

    // charge assignment polynomial coefficients rho_coeff[l][k-nlower] and influence function denominator expansion.
    // ad uses their derivative, LAMMPS drho_coeff[l-1][k] = l * rho_coeff[l][k] (see pppm_stencil_derivative)
    double rho_coeff[pppm_constants::MAXORDER][pppm_constants::MAXORDER] = {};
    double gf_b[pppm_constants::MAXORDER] = {};
    double sf_coeff[6] = {}; // ad : self force correction coefficients (LAMMPS sf_coeff)

    // per local k point (replicated : whole mesh, x fastest ; distributed : local rows, see PPPMDecomposition) :
    // influence function, k vectors (actual frame). Virial coefficients (LAMMPS vg) are computed from k when needed,
    // see pppm_virial_coeffs
    onika::memory::CudaMMVector<double> greensfn;
    onika::memory::CudaMMVector<double> fkx, fky, fkz;

    PPPMDecomposition dec;
    bool mesh_distributed_user = false; // distributed mesh in use (several ranks, mesh_decomposition not replicated)
    onika::memory::CudaMMVector<int> krow_partner; // distributed : local row of -y for each local row

    inline size_t nfft() const { return size_t(nx) * size_t(ny) * size_t(nz); }
    inline size_t nk_local() const { return dec.distributed ? size_t(dec.nyl()) * nx * nz : nfft(); }
    inline bool slab() const { return slab_auto || slab_user > 0.0; }
    inline int nlower() const { return -(order-1)/2; }
    inline int nupper() const { return order/2; }
  };

  // ------------------- LAMMPS lamda <-> box transposed transforms, on restricted cell parameters -------------------
  // lamda2xT uses absolute tilt values, as LAMMPS KSpace::lamda2xT
  inline void pppm_lamda2xT( const RestrictedCell& c, const double* l, double* v )
  {
    const double v0 = c.lx*l[0];
    const double v1 = std::fabs(c.xy)*l[0] + c.ly*l[1];
    const double v2 = std::fabs(c.xz)*l[0] + std::fabs(c.yz)*l[1] + c.lz*l[2];
    v[0]=v0; v[1]=v1; v[2]=v2;
  }

  // x2lamdaT = h_inv^T . v (LAMMPS Domain h_inv of the restricted cell). With the slab correction, LAMMPS divides the
  // z component by slab_volfactor (z extended cell, slab normal along z : xz = yz = 0)
  inline void pppm_x2lamdaT( const RestrictedCell& c, const double* v, double* l, double slab_volfactor = 1.0 )
  {
    const double hi0 = 1.0/c.lx, hi1 = 1.0/c.ly, hi2 = 1.0/c.lz;
    const double hi3 = -c.yz/(c.ly*c.lz);
    const double hi4 = (c.yz*c.xy - c.ly*c.xz)/(c.lx*c.ly*c.lz);
    const double hi5 = -c.xy/(c.lx*c.ly);
    const double l0 = hi0*v[0];
    const double l1 = hi5*v[0] + hi1*v[1];
    double l2 = hi4*v[0] + hi3*v[1] + hi2*v[2];
    if( slab_volfactor != 1.0 ) l2 /= slab_volfactor;
    l[0]=l0; l[1]=l1; l[2]=l2;
  }

  inline bool pppm_factorable( int n )
  {
    static constexpr int factors[3] = { 2, 3, 5 };
    while( n > 1 )
    {
      int i = 0;
      for( ; i < 3 ; i++ ) { if( n % factors[i] == 0 ) { n /= factors[i]; break; } }
      if( i == 3 ) return false;
    }
    return true;
  }

  // (sin(x)/x)^n, LAMMPS MathSpecial::powsinxx
  inline double pppm_powsinxx( double x, int n )
  {
    if( x == 0.0 ) return 1.0;
    double ww = std::sin(x) / x;
    double yy = 1.0;
    for( ; n != 0 ; n >>= 1, ww *= ww ) if( n & 1 ) yy *= ww;
    return yy;
  }

  // ------------------- accuracy estimates (LAMMPS metal units : q2 in eV.ang, forces in eV/ang) -------------------
  struct PPPMEstimator
  {
    RestrictedCell c;
    int order = 5;
    double g_ewald = 0.0;
    double cutoff = 0.0;
    double q2 = 0.0;         // sum q^2 * qqr2e
    double natoms = 1.0;
    double h_x = 0.0, h_y = 0.0, h_z = 0.0;
    bool diff_ad = false;
    int nx = 0, ny = 0, nz = 0; // mesh, used by the ad error estimate
    double slab_volfactor = 1.0;
    inline double zprd_slab() const { return c.lz*slab_volfactor; }

    // LAMMPS PPPM::compute_qopt (ad differentiation, orthogonal cell). qopt is a difference of nearly equal sums and
    // g_ewald comes from a finite difference derivative of it : terms are computed in parallel but summed in the LAMMPS
    // order, with the same operation order, so that g_ewald matches LAMMPS.
    inline double compute_qopt() const
    {
      static constexpr double MY_4PI = 4.0*M_PI;
      const double xprd = c.lx, yprd = c.ly, zprd = zprd_slab();
      const double unitkx = 2.0*M_PI/xprd, unitky = 2.0*M_PI/yprd, unitkz = 2.0*M_PI/zprd;
      const int twoorder = 2*order;
      const int64_t nxy = int64_t(nx) * ny;
      const int64_t ngrid = nxy * nz;
      std::vector<double> terms( ngrid , 0.0 );
#     pragma omp parallel for schedule(static)
      for( int64_t i = 0 ; i < ngrid ; i++ )
      {
        const int k = i % nx;
        const int l = (i/nx) % ny;
        const int m = i / nxy;
        const int kper = k - nx*(2*k/nx);
        const int lper = l - ny*(2*l/ny);
        const int mper = m - nz*(2*m/nz);
        const double sqk = (unitkx*kper)*(unitkx*kper) + (unitky*lper)*(unitky*lper) + (unitkz*mper)*(unitkz*mper);
        if( sqk == 0.0 ) continue;
        double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
        for( int ax = -2 ; ax <= 2 ; ax++ )
        {
          double qx = unitkx*(kper+nx*ax);
          const double sx = std::exp(-0.25*(qx/g_ewald)*(qx/g_ewald));
          const double wx = pppm_powsinxx( 0.5*qx*xprd/nx , twoorder );
          qx *= qx;
          for( int ay = -2 ; ay <= 2 ; ay++ )
          {
            double qy = unitky*(lper+ny*ay);
            const double sy = std::exp(-0.25*(qy/g_ewald)*(qy/g_ewald));
            const double wy = pppm_powsinxx( 0.5*qy*yprd/ny , twoorder );
            qy *= qy;
            for( int az = -2 ; az <= 2 ; az++ )
            {
              double qz = unitkz*(mper+nz*az);
              const double sz = std::exp(-0.25*(qz/g_ewald)*(qz/g_ewald));
              const double wz = pppm_powsinxx( 0.5*qz*zprd/nz , twoorder );
              qz *= qz;
              const double dot2 = qx+qy+qz;
              const double u1 = sx*sy*sz;
              const double u2 = wx*wy*wz;
              sum1 += u1*u1/dot2*MY_4PI*MY_4PI;
              sum2 += u1*u2*MY_4PI;
              sum3 += u2;
              sum4 += dot2*u2;
            }
          }
        }
        sum2 *= sum2;
        terms[i] = sum1 - sum2/(sum3*sum4);
      }
      double qopt = 0.0;
      for( int64_t i = 0 ; i < ngrid ; i++ ) qopt += terms[i];
      return qopt;
    }

    inline double estimate_ik_error( double h, double prd ) const
    {
      static constexpr double acons[8][7] = {
        { 0 },
        { 2.0/3.0 },
        { 1.0/50.0, 5.0/294.0 },
        { 1.0/588.0, 7.0/1440.0, 21.0/3872.0 },
        { 1.0/4320.0, 3.0/1936.0, 7601.0/2271360.0, 143.0/28800.0 },
        { 1.0/23232.0, 7601.0/13628160.0, 143.0/69120.0, 517231.0/106536960.0, 106640677.0/11737571328.0 },
        { 691.0/68140800.0, 13.0/57600.0, 47021.0/35512320.0, 9694607.0/2095994880.0, 733191589.0/59609088000.0, 326190917.0/11700633600.0 },
        { 1.0/345600.0, 3617.0/35512320.0, 745739.0/838397952.0, 56399353.0/12773376000.0, 25091609.0/1560084480.0, 1755948832039.0/36229939200000.0, 4887769399.0/37838389248.0 } };
      if( natoms == 0 ) return 0.0;
      double sum = 0.0;
      for( int m = 0 ; m < order ; m++ ) sum += acons[order][m] * std::pow( h*g_ewald , 2.0*m );
      return q2 * std::pow( h*g_ewald , (double)order ) * std::sqrt( g_ewald*prd*std::sqrt(2.0*M_PI)*sum/natoms ) / (prd*prd);
    }

    inline double df_kspace() const
    {
      if( diff_ad ) return std::sqrt( compute_qopt()/natoms ) * q2 / ( c.lx*c.ly*zprd_slab() );
      const double lprx = estimate_ik_error( h_x , c.lx );
      const double lpry = estimate_ik_error( h_y , c.ly );
      const double lprz = estimate_ik_error( h_z , zprd_slab() );
      return std::sqrt( lprx*lprx + lpry*lpry + lprz*lprz ) / std::sqrt(3.0);
    }

    inline double df_rspace() const
    {
      return 2.0*q2*std::exp(-g_ewald*g_ewald*cutoff*cutoff) / std::sqrt(natoms*cutoff*c.lx*c.ly*c.lz);
    }

    inline double newton_raphson_f() const { return df_rspace() - df_kspace(); }

    inline double derivf()
    {
      const double h = 0.000001;
      const double f1 = newton_raphson_f();
      const double g_old = g_ewald;
      g_ewald += h;
      const double f2 = newton_raphson_f();
      g_ewald = g_old;
      return (f2 - f1) / h;
    }

    // LAMMPS PPPM::adjust_gewald
    inline bool adjust_gewald()
    {
      using namespace pppm_constants;
      for( int i = 0 ; i < LARGE ; i++ )
      {
        const double dfx = derivf();
        if( dfx == 0.0 || dfx != dfx ) break;
        double dx = newton_raphson_f() / dfx;
        while( g_ewald - dx <= 0.0 ) dx *= 0.5;
        g_ewald -= dx;
        if( std::fabs( newton_raphson_f() ) < SMALL ) return true;
      }
      return false;
    }
  };

  // LAMMPS PPPM::set_grid_global (ik or ad differentiation, z extended by the slab volfactor) : chooses nx,ny,nz when not user defined, sets h_x,h_y,h_z
  inline void pppm_set_grid_global( PPPMEstimator& e, bool triclinic, const long mesh_user[3], int& nx, int& ny, int& nz, double accuracy )
  {
    const double xprd = e.c.lx, yprd = e.c.ly, zprd = e.c.lz;
    const double zprd_slab = e.zprd_slab();
    const bool gridflag = mesh_user[0] > 0 && mesh_user[1] > 0 && mesh_user[2] > 0;
    if( gridflag )
    {
      nx = mesh_user[0]; ny = mesh_user[1]; nz = mesh_user[2];
    }
    else if( e.diff_ad )
    {
      // LAMMPS : shrink the spacing until the qopt error estimate meets the accuracy (orthogonal cells only)
      double h = 4.0/e.g_ewald;
      int count = 0;
      while( true )
      {
        nx = static_cast<int>( xprd/h ); ny = static_cast<int>( yprd/h ); nz = static_cast<int>( zprd_slab/h );
        if( nx <= 1 ) nx = 2;
        if( ny <= 1 ) ny = 2;
        if( nz <= 1 ) nz = 2;
        e.nx = nx; e.ny = ny; e.nz = nz;
        const double df = e.df_kspace();
        count++;
        if( df <= accuracy ) break;
        if( count > 500 ) ::onika::fatal_error() << "PPPM : could not compute grid size" << std::endl;
        h *= 0.95;
      }
    }
    else
    {
      double err;
      e.h_x = e.h_y = e.h_z = 1.0/e.g_ewald;
      nx = static_cast<int>( xprd/e.h_x ) + 1;
      ny = static_cast<int>( yprd/e.h_y ) + 1;
      nz = static_cast<int>( zprd_slab/e.h_z ) + 1;

      // as in LAMMPS, err is evaluated before the increment (the final mesh is one point beyond the first passing one)
      err = e.estimate_ik_error( e.h_x , xprd );
      while( err > accuracy ) { err = e.estimate_ik_error( e.h_x , xprd ); nx++; e.h_x = xprd/nx; }
      err = e.estimate_ik_error( e.h_y , yprd );
      while( err > accuracy ) { err = e.estimate_ik_error( e.h_y , yprd ); ny++; e.h_y = yprd/ny; }
      err = e.estimate_ik_error( e.h_z , zprd_slab );
      while( err > accuracy ) { err = e.estimate_ik_error( e.h_z , zprd_slab ); nz++; e.h_z = zprd_slab/nz; }

      if( triclinic )
      {
        double tmp[3] = { nx/xprd , ny/yprd , nz/zprd };
        pppm_lamda2xT( e.c , tmp , tmp );
        nx = static_cast<int>(tmp[0]) + 1;
        ny = static_cast<int>(tmp[1]) + 1;
        nz = static_cast<int>(tmp[2]) + 1;
      }
    }

    while( ! pppm_factorable(nx) ) nx++;
    while( ! pppm_factorable(ny) ) ny++;
    while( ! pppm_factorable(nz) ) nz++;
    e.nx = nx; e.ny = ny; e.nz = nz;

    if( ! triclinic )
    {
      e.h_x = xprd/nx; e.h_y = yprd/ny; e.h_z = zprd_slab/nz;
    }
    else
    {
      double tmp[3] = { double(nx) , double(ny) , double(nz) };
      pppm_x2lamdaT( e.c , tmp , tmp , e.slab_volfactor );
      e.h_x = 1.0/tmp[0]; e.h_y = 1.0/tmp[1]; e.h_z = 1.0/tmp[2];
    }

    if( nx >= pppm_constants::OFFSET || ny >= pppm_constants::OFFSET || nz >= pppm_constants::OFFSET )
    {
      ::onika::fatal_error() << "PPPM mesh is too large : "<<nx<<"x"<<ny<<"x"<<nz << std::endl;
    }
  }

  // LAMMPS PPPM::compute_gf_denom
  inline void pppm_compute_gf_denom( PPPMParameters& p )
  {
    const int order = p.order;
    double* gf_b = p.gf_b;
    for( int l = 1 ; l < order ; l++ ) gf_b[l] = 0.0;
    gf_b[0] = 1.0;
    for( int m = 1 ; m < order ; m++ )
    {
      int l = m;
      for( ; l > 0 ; l-- ) gf_b[l] = 4.0 * ( gf_b[l]*(l-m)*(l-m-0.5) - gf_b[l-1]*(l-m-1)*(l-m-1) );
      gf_b[0] = 4.0 * ( gf_b[0]*(l-m)*(l-m-0.5) );
    }
    int64_t ifact = 1;
    for( int k = 1 ; k < 2*order ; k++ ) ifact *= k;
    const double gaminv = 1.0/ifact;
    for( int l = 0 ; l < order ; l++ ) gf_b[l] *= gaminv;
  }

  inline double pppm_gf_denom( const PPPMParameters& p, double x, double y, double z )
  {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    for( int l = p.order-1 ; l >= 0 ; l-- )
    {
      sx = p.gf_b[l] + sx*x;
      sy = p.gf_b[l] + sy*y;
      sz = p.gf_b[l] + sz*z;
    }
    const double s = sx*sy*sz;
    return s*s;
  }

  // LAMMPS PPPM::compute_rho_coeff : rho_coeff[l][k-nlower] for k = nlower..nupper
  inline void pppm_compute_rho_coeff( PPPMParameters& p )
  {
    const int order = p.order;
    // a[l][k+order], k = -order..order
    double a[pppm_constants::MAXORDER][2*pppm_constants::MAXORDER+1] = {};
    auto A = [&a,order](int l, int k) -> double& { return a[l][k+order]; };
    A(0,0) = 1.0;
    for( int j = 1 ; j < order ; j++ )
    {
      for( int k = -j ; k <= j ; k += 2 )
      {
        double s = 0.0;
        for( int l = 0 ; l < j ; l++ )
        {
          A(l+1,k) = ( A(l,k+1) - A(l,k-1) ) / (l+1);
          s += std::pow(0.5,(double)l+1) * ( A(l,k-1) + std::pow(-1.0,(double)l) * A(l,k+1) ) / (l+1);
        }
        A(0,k) = s;
      }
    }
    int m = 0;
    for( int k = -(order-1) ; k < order ; k += 2 )
    {
      for( int l = 0 ; l < order ; l++ ) p.rho_coeff[l][m] = A(l,k);
      m++;
    }
  }

  // H^-T . v : reciprocal vector in the actual frame for a lamda-space wave vector v
  inline Vec3d pppm_reciprocal( const Mat3d& Hinv, double v0, double v1, double v2 )
  {
    return Vec3d{ Hinv.m11*v0 + Hinv.m21*v1 + Hinv.m31*v2 ,
                  Hinv.m12*v0 + Hinv.m22*v1 + Hinv.m32*v2 ,
                  Hinv.m13*v0 + Hinv.m23*v1 + Hinv.m33*v2 };
  }

  // calls f(n, k, l, m) for each local k point : n local index, (k,l,m) = global mesh indices along x,y,z
  template<class FuncT>
  inline void pppm_for_local_kpoints( const PPPMParameters& p, const FuncT& f )
  {
    const int nx = p.nx, ny = p.ny, nz = p.nz;
    if( ! p.dec.distributed )
    {
#     pragma omp parallel for schedule(static)
      for( int m = 0 ; m < nz ; m++ )
        for( int l = 0 ; l < ny ; l++ )
          for( int k = 0 ; k < nx ; k++ ) f( ( size_t(m)*ny + l ) * nx + k , k , l , m );
    }
    else
    {
      const int nyl = p.dec.nyl();
      const int* rows = p.dec.rows();
#     pragma omp parallel for schedule(static)
      for( int li = 0 ; li < nyl ; li++ )
        for( int k = 0 ; k < nx ; k++ )
          for( int m = 0 ; m < nz ; m++ ) f( ( size_t(li)*nx + k ) * nz + m , k , rows[li] , m );
    }
  }

  // LAMMPS PPPM::setup + compute_sf_precoeff + compute_gf_ad (ad differentiation, orthogonal cell only)
  inline void pppm_setup_ad( PPPMParameters& p, const RestrictedCell& c, MPI_Comm comm )
  {
    const int nx = p.nx, ny = p.ny, nz = p.nz;
    const double xprd = c.lx, yprd = c.ly, zprd = c.lz*p.slab_volfactor; // zprd_slab
    const double unitkx = 2.0*M_PI/xprd, unitky = 2.0*M_PI/yprd, unitkz = 2.0*M_PI/zprd;
    const double g = p.g_ewald;
    const int order = p.order;
    const int twoorder = 2*order;
    const size_t nk = p.nk_local();
    std::vector<double> sfterms( 6*nk , 0.0 ); // per point self force terms, summed in point order afterwards (as LAMMPS)

    pppm_for_local_kpoints( p , [&]( size_t n, int k, int l, int m )
    {
      const int mper = m - nz*(2*m/nz);
      const double qz = unitkz*mper;
      const double snz = std::pow( std::sin(0.5*qz*zprd/nz) , 2 );
      const double sz = std::exp(-0.25*(qz/g)*(qz/g));
      const double wz = pppm_powsinxx( 0.5*qz*zprd/nz , twoorder );
      const int lper = l - ny*(2*l/ny);
      const double qy = unitky*lper;
      const double sny = std::pow( std::sin(0.5*qy*yprd/ny) , 2 );
      const double sy = std::exp(-0.25*(qy/g)*(qy/g));
      const double wy = pppm_powsinxx( 0.5*qy*yprd/ny , twoorder );
      const int kper = k - nx*(2*k/nx);
      const double qx = unitkx*kper;
      const double snx = std::pow( std::sin(0.5*qx*xprd/nx) , 2 );
      const double sx = std::exp(-0.25*(qx/g)*(qx/g));
      const double wx = pppm_powsinxx( 0.5*qx*xprd/nx , twoorder );
      p.fkx[n] = qx; p.fky[n] = qy; p.fkz[n] = qz;
      const double sqk = qx*qx + qy*qy + qz*qz;
      if( sqk == 0.0 ) { p.greensfn[n] = 0.0; return; }
      const double numerator = (4.0*M_PI)/sqk;
      const double gf = numerator*sx*sy*sz*wx*wy*wz/pppm_gf_denom( p , snx , sny , snz );
      p.greensfn[n] = gf;

      // self force pre-coefficients of this mesh point (LAMMPS compute_sf_precoeff)
      double wx0[5], wy0[5], wz0[5], wx1[5], wy1[5], wz1[5], wx2[5], wy2[5], wz2[5];
      for( int i = 0 ; i < 5 ; i++ )
      {
        wx0[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(kper+nx*(i-2))/nx , order );
        wx1[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(kper+nx*(i-1))/nx , order );
        wx2[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(kper+nx*i)/nx , order );
        wy0[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(lper+ny*(i-2))/ny , order );
        wy1[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(lper+ny*(i-1))/ny , order );
        wy2[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(lper+ny*i)/ny , order );
        wz0[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(mper+nz*(i-2))/nz , order );
        wz1[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(mper+nz*(i-1))/nz , order );
        wz2[i] = pppm_powsinxx( 0.5*(2.0*M_PI)*(mper+nz*i)/nz , order );
      }
      double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0, s6 = 0.0;
      for( int ax = 0 ; ax < 5 ; ax++ )
        for( int ay = 0 ; ay < 5 ; ay++ )
          for( int az = 0 ; az < 5 ; az++ )
          {
            const double u0 = wx0[ax]*wy0[ay]*wz0[az];
            s1 += u0 * wx1[ax]*wy0[ay]*wz0[az];
            s2 += u0 * wx2[ax]*wy0[ay]*wz0[az];
            s3 += u0 * wx0[ax]*wy1[ay]*wz0[az];
            s4 += u0 * wx0[ax]*wy2[ay]*wz0[az];
            s5 += u0 * wx0[ax]*wy0[ay]*wz1[az];
            s6 += u0 * wx0[ax]*wy0[ay]*wz2[az];
          }
      double* t = sfterms.data() + 6*n;
      t[0] = s1*gf; t[1] = s2*gf; t[2] = s3*gf; t[3] = s4*gf; t[4] = s5*gf; t[5] = s6*gf;
    } );

    double sf[6] = { 0., 0., 0., 0., 0., 0. };
    for( size_t n = 0 ; n < nk ; n++ ) for( int i = 0 ; i < 6 ; i++ ) sf[i] += sfterms[6*n+i];
    if( p.dec.distributed ) MPI_Allreduce( MPI_IN_PLACE , sf , 6 , MPI_DOUBLE , MPI_SUM , comm );

    const double pre = M_PI / p.volume;
    const double prex = pre * nx / xprd, prey = pre * ny / yprd, prez = pre * nz / zprd;
    p.sf_coeff[0] = sf[0] * prex; p.sf_coeff[1] = sf[1] * prex * 2;
    p.sf_coeff[2] = sf[2] * prey; p.sf_coeff[3] = sf[3] * prey * 2;
    p.sf_coeff[4] = sf[4] * prez; p.sf_coeff[5] = sf[5] * prez * 2;
  }

  // LAMMPS PPPM::setup_triclinic + compute_gf_ik_triclinic (the orthogonal formulas are the same algebra) :
  // volume dependent quantities, called at init and whenever the cell changes (mesh and g_ewald are kept)
  inline void pppm_setup( PPPMParameters& p, const Mat3d& H, MPI_Comm comm )
  {
    using namespace pppm_constants;
    const RestrictedCell c = restricted_cell( H );
    const Mat3d Hinv = inverse( H );
    const int nx = p.nx, ny = p.ny, nz = p.nz;
    const size_t nk = p.nk_local();
    const double g = p.g_ewald;

    const double vf = p.slab_volfactor;
    p.cell = H;
    p.volume = c.lx * c.ly * ( c.lz * vf );

    p.greensfn.resize( nk );
    p.fkx.resize( nk ); p.fky.resize( nk ); p.fkz.resize( nk );
    if( p.dec.distributed )
    {
      const int nyl = p.dec.nyl();
      p.krow_partner.resize( nyl );
      for( int li = 0 ; li < nyl ; li++ ) p.krow_partner[li] = p.dec.row_local[ ( ny - p.dec.rows()[li] ) % ny ];
    }
    else p.krow_partner.clear();

    if( p.diff_ad )
    {
      if( ! is_diagonal( H ) ) ::onika::fatal_error() << "PPPM : diff ad requires an orthogonal cell (as in LAMMPS)" << std::endl;
      pppm_setup_ad( p , c , comm );
      return;
    }

    double tmp[3] = { (g/(M_PI*nx)) * std::pow(-std::log(EPS_HOC),0.25) ,
                      (g/(M_PI*ny)) * std::pow(-std::log(EPS_HOC),0.25) ,
                      (g/(M_PI*nz)) * std::pow(-std::log(EPS_HOC),0.25) };
    pppm_lamda2xT( c , tmp , tmp );
    tmp[2] *= vf; // slab : alias sum bound on the extended z length, as LAMMPS
    const int nbx = static_cast<int>( tmp[0] );
    const int nby = static_cast<int>( tmp[1] );
    const int nbz = static_cast<int>( tmp[2] );
    const int twoorder = 2 * p.order;

    pppm_for_local_kpoints( p , [&]( size_t n, int k, int l, int m )
    {
      const int mper = m - nz*(2*m/nz);
      const double snz = std::pow( std::sin(M_PI*mper/nz) , 2 );
      const int lper = l - ny*(2*l/ny);
      const double sny = std::pow( std::sin(M_PI*lper/ny) , 2 );
      const int kper = k - nx*(2*k/nx);
      const double snx = std::pow( std::sin(M_PI*kper/nx) , 2 );

      Vec3d fk = pppm_reciprocal( Hinv , 2.0*M_PI*kper , 2.0*M_PI*lper , 2.0*M_PI*mper );
      if( vf != 1.0 ) fk.z /= vf; // slab : z extended cell (c along z, see pppm_init_parameters)
      p.fkx[n] = fk.x; p.fky[n] = fk.y; p.fkz[n] = fk.z;
      const double sqk = fk.x*fk.x + fk.y*fk.y + fk.z*fk.z;

      if( sqk == 0.0 )
      {
        p.greensfn[n] = 0.0;
        return;
      }

      const double numerator = FOUR_PI_LMP / sqk;
      const double denominator = pppm_gf_denom( p , snx , sny , snz );
      double sum1 = 0.0;
      for( int ax = -nbx ; ax <= nbx ; ax++ )
      {
        const double wx = pppm_powsinxx( M_PI*kper/nx + M_PI*ax , twoorder );
        for( int ay = -nby ; ay <= nby ; ay++ )
        {
          const double wy = pppm_powsinxx( M_PI*lper/ny + M_PI*ay , twoorder );
          for( int az = -nbz ; az <= nbz ; az++ )
          {
            const double wz = pppm_powsinxx( M_PI*mper/nz + M_PI*az , twoorder );
            Vec3d b = pppm_reciprocal( Hinv , 2.0*M_PI*nx*ax , 2.0*M_PI*ny*ay , 2.0*M_PI*nz*az );
            if( vf != 1.0 ) b.z /= vf;
            const double qx = fk.x + b.x;
            const double qy = fk.y + b.y;
            const double qz = fk.z + b.z;
            const double sx = std::exp( -0.25 * (qx/g)*(qx/g) );
            const double sy = std::exp( -0.25 * (qy/g)*(qy/g) );
            const double sz = std::exp( -0.25 * (qz/g)*(qz/g) );
            const double dot1 = fk.x*qx + fk.y*qy + fk.z*qz;
            const double dot2 = qx*qx + qy*qy + qz*qz;
            sum1 += ( dot1/dot2 ) * sx*sy*sz * wx*wy*wz;
          }
        }
      }
      p.greensfn[n] = numerator * sum1 / denominator;
    } );
  }

  // LAMMPS auto_slab_volfactor (kspace_modify slab auto) : vacuum large enough for the lateral and reciprocal decay lengths
  inline double pppm_auto_slab_volfactor( double force_tolerance, double alpha, double xprd, double yprd, double zprd )
  {
    if( alpha <= 0.0 ) ::onika::fatal_error() << "PPPM slab_auto requires a positive g_ewald" << std::endl;
    if( !( force_tolerance > 0.0 && force_tolerance < 1.0 ) ) ::onika::fatal_error() << "PPPM slab_auto requires accuracy_relative between 0 and 1" << std::endl;
    const double logeps = std::log( 1.0 / force_tolerance );
    const double lateral = std::max( xprd , yprd ) * logeps / (2.0*M_PI);
    const double reciprocal = std::sqrt( logeps ) / alpha;
    return std::max( ( zprd + std::max( lateral , reciprocal ) ) / zprd , 1.0 );
  }

  // LAMMPS PPPM::init (ik or ad, slab EW3DC with fixed or automatic volfactor, no tip4p) followed by setup.
  // qsqsum = sum of q^2, qsum = sum of q (elementary charges). Estimates use LAMMPS metal units (eV, ang).
  inline void pppm_init_parameters( double g_ewald, double radius, double accuracy_relative, long order, const long mesh[3], bool diff_ad,
                                    double slab_user, bool slab_auto, bool mesh_distributed,
                                    const Mat3d& H, uint64_t natoms, double qsqsum, double qsum, MPI_Comm comm, PPPMParameters& p )
  {
    if( order < 2 || order > pppm_constants::MAXORDER )
    {
      ::onika::fatal_error() << "PPPM order must be in [2,"<<pppm_constants::MAXORDER<<"], got "<<order << std::endl;
    }

    p.accuracy_relative = accuracy_relative;
    p.g_ewald_user = g_ewald;
    p.radius = radius;
    p.order = order;
    p.diff_ad = diff_ad;
    p.slab_user = slab_user;
    p.slab_auto = slab_auto;
    for( int i = 0 ; i < 3 ; i++ ) p.mesh_user[i] = mesh[i];
    p.qsum = qsum;
    p.qsqsum = qsqsum;

    const bool triclinic = ! is_diagonal( H );
    if( diff_ad && triclinic ) ::onika::fatal_error() << "PPPM : diff ad requires an orthogonal cell (as in LAMMPS)" << std::endl;
    if( slab_user != 0.0 && slab_user <= 1.0 ) ::onika::fatal_error() << "PPPM : slab volfactor must be > 1, got "<<slab_user << std::endl;
    if( slab_user > 0.0 && slab_auto ) ::onika::fatal_error() << "PPPM : slab and slab_auto are exclusive" << std::endl;
    if( p.slab() && ( H.m13 != 0.0 || H.m23 != 0.0 || H.m31 != 0.0 || H.m32 != 0.0 ) )
    {
      ::onika::fatal_error() << "PPPM : slab correction requires the third cell vector along z and the two others in the xy plane (xz = yz = 0)" << std::endl;
    }
    PPPMEstimator e;
    e.diff_ad = diff_ad;
    e.c = restricted_cell( H );
    e.order = order;
    e.cutoff = radius;
    e.q2 = qsqsum * COULOMB_CONSTANT_EV_ANG;
    e.natoms = natoms == 0 ? 1.0 : double(natoms);
    const double accuracy = accuracy_relative * COULOMB_CONSTANT_EV_ANG; // two_charge_force in metal units

    const bool gewaldflag = g_ewald > 0.0;
    if( gewaldflag ) e.g_ewald = g_ewald;
    else
    {
      if( accuracy <= 0.0 ) ::onika::fatal_error() << "PPPM accuracy_relative must be > 0" << std::endl;
      if( e.q2 == 0.0 ) ::onika::fatal_error() << "PPPM : g_ewald must be given for an uncharged system" << std::endl;
      double g = accuracy*std::sqrt(e.natoms*radius*e.c.lx*e.c.ly*e.c.lz) / (2.0*e.q2);
      if( g >= 1.0 ) g = (1.35 - 0.15*std::log(accuracy))/radius;
      else g = std::sqrt(-std::log(g)) / radius;
      e.g_ewald = g;
    }

    // slab auto : volfactor depends on g_ewald, which depends on the mesh, which depends on volfactor (LAMMPS loop)
    const double force_tolerance = accuracy / COULOMB_CONSTANT_EV_ANG;
    double vf = slab_user > 0.0 ? slab_user : 1.0;
    if( slab_auto ) vf = pppm_auto_slab_volfactor( force_tolerance , e.g_ewald , e.c.lx , e.c.ly , e.c.lz );
    int slab_iterations = 0;
    while( true )
    {
      e.slab_volfactor = vf;
      pppm_set_grid_global( e , triclinic , mesh , p.nx , p.ny , p.nz , accuracy );
      if( p.nx < order || p.ny < order || p.nz < order )
      {
        ::onika::fatal_error() << "PPPM mesh "<<p.nx<<"x"<<p.ny<<"x"<<p.nz<<" must have at least order="<<order<<" points in each direction" << std::endl;
      }
      if( ! gewaldflag )
      {
        if( ! e.adjust_gewald() ) ::onika::fatal_error() << "PPPM : could not compute g_ewald" << std::endl;
      }
      if( ! slab_auto ) break;
      const double new_vf = pppm_auto_slab_volfactor( force_tolerance , e.g_ewald , e.c.lx , e.c.ly , e.c.lz );
      if( std::fabs( new_vf - vf ) <= pppm_constants::SMALL * new_vf ) break;
      vf = new_vf;
      if( ++slab_iterations > 5 ) ::onika::fatal_error() << "PPPM : could not converge slab_auto" << std::endl;
    }
    p.slab_volfactor = vf;
    p.g_ewald = e.g_ewald;

    // final accuracy estimate (no coulomb table)
    const double dfk = e.df_kspace();
    const double dfr = e.df_rspace();
    p.estimated_accuracy = std::sqrt( dfk*dfk + dfr*dfr );

    int nprocs = 1, rank = 0;
    MPI_Comm_size( comm , &nprocs );
    MPI_Comm_rank( comm , &rank );
    p.mesh_distributed_user = mesh_distributed;
    pppm_make_decomposition( p.dec , mesh_distributed , nprocs , rank , p.nx , p.ny , p.nz );

    pppm_compute_gf_denom( p );
    pppm_compute_rho_coeff( p );
    pppm_setup( p , H , comm );
  }

  // trivially copyable view of PPPMParameters for particle <-> mesh functors
  // LAMMPS PPPM vg : virial coefficients (xx,yy,zz,xy,xz,yz) of mesh point with k vector (kx,ky,kz)
  ONIKA_HOST_DEVICE_FUNC inline void pppm_virial_coeffs( double kx, double ky, double kz, double g, double vg[6] )
  {
    const double sqk = kx*kx + ky*ky + kz*kz;
    if( sqk == 0.0 )
    {
      for( int i = 0 ; i < 6 ; i++ ) vg[i] = 0.0;
      return;
    }
    const double vterm = -2.0 * ( 1.0/sqk + 0.25/(g*g) );
    vg[0] = 1.0 + vterm*kx*kx;
    vg[1] = 1.0 + vterm*ky*ky;
    vg[2] = 1.0 + vterm*kz*kz;
    vg[3] = vterm*kx*ky;
    vg[4] = vterm*kx*kz;
    vg[5] = vterm*ky*kz;
  }

  struct ReadOnlyPPPMParameters
  {
    int order = 5;
    int nlower = -2;
    int nx = 0, ny = 0, nz = 0;
    double shift = 0.0;     // OFFSET (+0.5 for odd orders)
    double shiftone = 0.0;  // 0 for odd orders, 0.5 for even
    double g_ewald = 0.0;
    double qsum = 0.0;
    double volume = 0.0;
    double delvolinv = 0.0; // number of mesh points / volume
    Vec3d bmin = {0.,0.,0.};
    Vec3d delinv = {0.,0.,0.}; // mesh points per grid-space length unit, along each fractional axis
    double rho_coeff[pppm_constants::MAXORDER][pppm_constants::MAXORDER] = {};
    // ad differentiation (orthogonal cell) : mesh points per length unit (LAMMPS hx_inv...),
    // diagonal of xform (real position = xf * grid position) and self force coefficients
    bool diff_ad = false;
    Vec3d hinv = {0.,0.,0.};
    Vec3d xf = {1.,1.,1.};
    double sf_coeff[6] = {};
    // slab correction (EW3DC) : real z = zscale * grid z, extended height, total dipole sum(q.z) and sum(q.z^2)
    bool slab = false;
    double zscale = 1.0;
    double zprd_slab = 0.0;
    double dipole = 0.0;
    double dipole_r2 = 0.0;
    // mesh read/written by the particle functors : whole periodic mesh, or a local brick holding the stencils of the
    // local particles (distributed mesh). In each direction, the brick either covers the whole periodic mesh
    // (mwrap, indices wrapped into [0,n)) or a range of unwrapped mesh indices [mo, mo+mn), mn < n.
    int mnx = 0, mny = 0, mnz = 0;
    int mox = 0, moy = 0, moz = 0;
    bool mwrapx = true, mwrapy = true, mwrapz = true;

    // a brick dimension equal to the mesh size is a whole wrapped direction (lo must then be 0)
    inline void set_brick( const int lo[3], const int dims[3] )
    {
      mox = lo[0]; moy = lo[1]; moz = lo[2];
      mnx = dims[0]; mny = dims[1]; mnz = dims[2];
      mwrapx = ( mnx == nx );
      mwrapy = ( mny == ny );
      mwrapz = ( mnz == nz );
    }

    ReadOnlyPPPMParameters() = default;
    inline ReadOnlyPPPMParameters( const PPPMParameters& p, const Vec3d& bounds_min, const Vec3d& bounds_size )
      : order( p.order ), nlower( p.nlower() ), nx( p.nx ), ny( p.ny ), nz( p.nz )
      , shift( pppm_constants::OFFSET + ( (p.order % 2) ? 0.5 : 0.0 ) )
      , shiftone( (p.order % 2) ? 0.0 : 0.5 )
      , g_ewald( p.g_ewald ), qsum( p.qsum ), volume( p.volume )
      , delvolinv( double(p.nfft()) / p.volume )
      , bmin( bounds_min )
      , delinv( Vec3d{ p.nx / bounds_size.x , p.ny / bounds_size.y , p.nz / ( bounds_size.z * p.slab_volfactor ) } )
      , mnx( p.nx ), mny( p.ny ), mnz( p.nz )
    {
      for( int l = 0 ; l < pppm_constants::MAXORDER ; l++ )
        for( int k = 0 ; k < pppm_constants::MAXORDER ; k++ ) rho_coeff[l][k] = p.rho_coeff[l][k];
      if( p.diff_ad )
      {
        diff_ad = true;
        // z : extended height with the slab correction. LAMMPS fieldforce_ad uses nz/zprd here, which scales the
        // z field by slab_volfactor (wrong forces with diff ad + slab) ; the mesh spacing is zprd_slab/nz.
        hinv = Vec3d{ p.nx / p.cell.m11 , p.ny / p.cell.m22 , p.nz / ( p.cell.m33 * p.slab_volfactor ) };
        xf = Vec3d{ p.cell.m11 / bounds_size.x , p.cell.m22 / bounds_size.y , p.cell.m33 / bounds_size.z };
        for( int i = 0 ; i < 6 ; i++ ) sf_coeff[i] = p.sf_coeff[i];
      }
      if( p.slab() )
      {
        slab = true;
        zscale = p.cell.m33 / bounds_size.z;
        zprd_slab = restricted_cell( p.cell ).lz * p.slab_volfactor;
      }
    }
  };

  // stencil of one particle : lower mesh index in the mesh view (see ReadOnlyPPPMParameters::set_brick) and 1D weights
  struct PPPMStencil
  {
    int ix, iy, iz;
    double wx[pppm_constants::MAXORDER];
    double wy[pppm_constants::MAXORDER];
    double wz[pppm_constants::MAXORDER];
    double dx, dy, dz; // distance to the lower left mesh point, in mesh units (LAMMPS dx,dy,dz)
    int gx0, gy0, gz0; // unwrapped global mesh index of the first stencil point
  };

  // LAMMPS compute_drho1d : derivative of the 1D weights with respect to dx,dy,dz, drho_coeff[l][k] = (l+1).rho_coeff[l+1][k]
  ONIKA_HOST_DEVICE_FUNC inline void pppm_stencil_derivative( const ReadOnlyPPPMParameters& p, const PPPMStencil& st,
                                                             double dwx[pppm_constants::MAXORDER], double dwy[pppm_constants::MAXORDER], double dwz[pppm_constants::MAXORDER] )
  {
    for( int k = 0 ; k < p.order ; k++ )
    {
      double r1 = 0.0, r2 = 0.0, r3 = 0.0;
      for( int l = p.order-2 ; l >= 0 ; l-- )
      {
        const double c = (l+1) * p.rho_coeff[l+1][k];
        r1 = c + r1*st.dx;
        r2 = c + r2*st.dy;
        r3 = c + r3*st.dz;
      }
      dwx[k] = r1; dwy[k] = r2; dwz[k] = r3;
    }
  }

  ONIKA_HOST_DEVICE_FUNC inline int pppm_wrap( int i, int n ) { i %= n; return i < 0 ? i + n : i; }

  // LAMMPS particle_map + compute_rho1d. r is the particle position in grid space (real = xform . r)
  ONIKA_HOST_DEVICE_FUNC inline void pppm_stencil( const ReadOnlyPPPMParameters& p, const Vec3d& r, PPPMStencil& st )
  {
    using pppm_constants::OFFSET;
    const double ux = ( r.x - p.bmin.x ) * p.delinv.x;
    const double uy = ( r.y - p.bmin.y ) * p.delinv.y;
    const double uz = ( r.z - p.bmin.z ) * p.delinv.z;
    const int gx = static_cast<int>( ux + p.shift ) - OFFSET;
    const int gy = static_cast<int>( uy + p.shift ) - OFFSET;
    const int gz = static_cast<int>( uz + p.shift ) - OFFSET;
    const double dx = gx + p.shiftone - ux;
    const double dy = gy + p.shiftone - uy;
    const double dz = gz + p.shiftone - uz;
    for( int k = 0 ; k < p.order ; k++ )
    {
      double r1 = 0.0, r2 = 0.0, r3 = 0.0;
      for( int l = p.order-1 ; l >= 0 ; l-- )
      {
        r1 = p.rho_coeff[l][k] + r1*dx;
        r2 = p.rho_coeff[l][k] + r2*dy;
        r3 = p.rho_coeff[l][k] + r3*dz;
      }
      st.wx[k] = r1; st.wy[k] = r2; st.wz[k] = r3;
    }
    st.dx = dx; st.dy = dy; st.dz = dz;
    st.gx0 = gx + p.nlower; st.gy0 = gy + p.nlower; st.gz0 = gz + p.nlower;
    st.ix = p.mwrapx ? pppm_wrap( st.gx0 , p.nx ) : st.gx0 - p.mox;
    st.iy = p.mwrapy ? pppm_wrap( st.gy0 , p.ny ) : st.gy0 - p.moy;
    st.iz = p.mwrapz ? pppm_wrap( st.gz0 , p.nz ) : st.gz0 - p.moz;
  }

  // LAMMPS PPPM::slabcorr (EW3DC) for one particle at real height z : force along z and per particle energy, internal units
  ONIKA_HOST_DEVICE_FUNC inline double pppm_slab_force_z( const ReadOnlyPPPMParameters& p, double q, double z )
  {
    return COULOMB_CONSTANT * ( -4.0*M_PI / p.volume ) * q * ( p.dipole - p.qsum * z );
  }
  ONIKA_HOST_DEVICE_FUNC inline double pppm_slab_energy( const ReadOnlyPPPMParameters& p, double q, double z )
  {
    return COULOMB_CONSTANT * (2.0*M_PI) / p.volume * q * ( z * p.dipole - 0.5 * ( p.dipole_r2 + p.qsum * z * z ) - p.qsum * p.zprd_slab * p.zprd_slab / 12.0 );
  }

  // self energy + neutralizing background energy of one particle with charge q (LAMMPS PPPM per atom correction), internal units
  ONIKA_HOST_DEVICE_FUNC inline double pppm_self_energy( const ReadOnlyPPPMParameters& p, double q )
  {
    return - COULOMB_CONSTANT * ( p.g_ewald / sqrt(M_PI) * q * q + 0.5 * M_PI * q * p.qsum / ( p.g_ewald * p.g_ewald * p.volume ) );
  }

}
}
