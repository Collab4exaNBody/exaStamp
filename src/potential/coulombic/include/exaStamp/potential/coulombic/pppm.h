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
#include <cstdint>

namespace exaStamp
{
inline namespace coulombic_ewald // distinct symbols from the legacy ewald plugin (plugins are loaded RTLD_GLOBAL)
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

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PPPMParameters
  {
    // user settings
    double accuracy_relative = 0.0;
    double g_ewald_user = 0.0;  // 0 = automatic
    double radius = 0.0;
    long order = 5;
    long mesh_user[3] = { 0, 0, 0 }; // 0 = automatic

    // derived
    double g_ewald = 0.0;
    int nx = 0, ny = 0, nz = 0;
    Mat3d cell = { 0.,0.,0., 0.,0.,0., 0.,0.,0. }; // cell matrix H the influence function was built for
    double volume = 0.0;
    double qsum = 0.0;
    double qsqsum = 0.0;
    double estimated_accuracy = 0.0; // absolute rms force error estimate, eV/ang

    // charge assignment polynomial coefficients rho_coeff[l][k-nlower] and influence function denominator expansion
    double rho_coeff[pppm_constants::MAXORDER][pppm_constants::MAXORDER] = {};
    double gf_b[pppm_constants::MAXORDER] = {};

    // per mesh point (x fastest) : influence function, k vectors (actual frame).
    // virial coefficients (LAMMPS vg) are computed from k when needed, see pppm_virial_coeffs
    onika::memory::CudaMMVector<double> greensfn;
    onika::memory::CudaMMVector<double> fkx, fky, fkz;

    inline size_t nfft() const { return size_t(nx) * size_t(ny) * size_t(nz); }
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

  // x2lamdaT = h_inv^T . v (LAMMPS Domain h_inv of the restricted cell)
  inline void pppm_x2lamdaT( const RestrictedCell& c, const double* v, double* l )
  {
    const double hi0 = 1.0/c.lx, hi1 = 1.0/c.ly, hi2 = 1.0/c.lz;
    const double hi3 = -c.yz/(c.ly*c.lz);
    const double hi4 = (c.yz*c.xy - c.ly*c.xz)/(c.lx*c.ly*c.lz);
    const double hi5 = -c.xy/(c.lx*c.ly);
    const double l0 = hi0*v[0];
    const double l1 = hi5*v[0] + hi1*v[1];
    const double l2 = hi4*v[0] + hi3*v[1] + hi2*v[2];
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
      const double lprx = estimate_ik_error( h_x , c.lx );
      const double lpry = estimate_ik_error( h_y , c.ly );
      const double lprz = estimate_ik_error( h_z , c.lz );
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

  // LAMMPS PPPM::set_grid_global (ik differentiation, no slab) : chooses nx,ny,nz when not user defined, sets h_x,h_y,h_z
  inline void pppm_set_grid_global( PPPMEstimator& e, bool triclinic, const long mesh_user[3], int& nx, int& ny, int& nz, double accuracy )
  {
    const double xprd = e.c.lx, yprd = e.c.ly, zprd = e.c.lz;
    const bool gridflag = mesh_user[0] > 0 && mesh_user[1] > 0 && mesh_user[2] > 0;
    if( gridflag )
    {
      nx = mesh_user[0]; ny = mesh_user[1]; nz = mesh_user[2];
    }
    else
    {
      double err;
      e.h_x = e.h_y = e.h_z = 1.0/e.g_ewald;
      nx = static_cast<int>( xprd/e.h_x ) + 1;
      ny = static_cast<int>( yprd/e.h_y ) + 1;
      nz = static_cast<int>( zprd/e.h_z ) + 1;

      // as in LAMMPS, err is evaluated before the increment (the final mesh is one point beyond the first passing one)
      err = e.estimate_ik_error( e.h_x , xprd );
      while( err > accuracy ) { err = e.estimate_ik_error( e.h_x , xprd ); nx++; e.h_x = xprd/nx; }
      err = e.estimate_ik_error( e.h_y , yprd );
      while( err > accuracy ) { err = e.estimate_ik_error( e.h_y , yprd ); ny++; e.h_y = yprd/ny; }
      err = e.estimate_ik_error( e.h_z , zprd );
      while( err > accuracy ) { err = e.estimate_ik_error( e.h_z , zprd ); nz++; e.h_z = zprd/nz; }

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

    if( ! triclinic )
    {
      e.h_x = xprd/nx; e.h_y = yprd/ny; e.h_z = zprd/nz;
    }
    else
    {
      double tmp[3] = { double(nx) , double(ny) , double(nz) };
      pppm_x2lamdaT( e.c , tmp , tmp );
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

  // LAMMPS PPPM::setup_triclinic + compute_gf_ik_triclinic (the orthogonal formulas are the same algebra) :
  // volume dependent quantities, called at init and whenever the cell changes (mesh and g_ewald are kept)
  inline void pppm_setup( PPPMParameters& p, const Mat3d& H )
  {
    using namespace pppm_constants;
    const RestrictedCell c = restricted_cell( H );
    const Mat3d Hinv = inverse( H );
    const int nx = p.nx, ny = p.ny, nz = p.nz;
    const size_t nfft = p.nfft();
    const double g = p.g_ewald;

    p.cell = H;
    p.volume = c.lx * c.ly * c.lz;

    p.greensfn.resize( nfft );
    p.fkx.resize( nfft ); p.fky.resize( nfft ); p.fkz.resize( nfft );

    double tmp[3] = { (g/(M_PI*nx)) * std::pow(-std::log(EPS_HOC),0.25) ,
                      (g/(M_PI*ny)) * std::pow(-std::log(EPS_HOC),0.25) ,
                      (g/(M_PI*nz)) * std::pow(-std::log(EPS_HOC),0.25) };
    pppm_lamda2xT( c , tmp , tmp );
    const int nbx = static_cast<int>( tmp[0] );
    const int nby = static_cast<int>( tmp[1] );
    const int nbz = static_cast<int>( tmp[2] );
    const int twoorder = 2 * p.order;

#   pragma omp parallel for schedule(static)
    for( int m = 0 ; m < nz ; m++ )
    {
      const int mper = m - nz*(2*m/nz);
      const double snz = std::pow( std::sin(M_PI*mper/nz) , 2 );
      for( int l = 0 ; l < ny ; l++ )
      {
        const int lper = l - ny*(2*l/ny);
        const double sny = std::pow( std::sin(M_PI*lper/ny) , 2 );
        for( int k = 0 ; k < nx ; k++ )
        {
          const size_t n = ( size_t(m)*ny + l ) * nx + k;
          const int kper = k - nx*(2*k/nx);
          const double snx = std::pow( std::sin(M_PI*kper/nx) , 2 );

          const Vec3d fk = pppm_reciprocal( Hinv , 2.0*M_PI*kper , 2.0*M_PI*lper , 2.0*M_PI*mper );
          p.fkx[n] = fk.x; p.fky[n] = fk.y; p.fkz[n] = fk.z;
          const double sqk = fk.x*fk.x + fk.y*fk.y + fk.z*fk.z;

          if( sqk == 0.0 )
          {
            p.greensfn[n] = 0.0;
            continue;
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
                const Vec3d b = pppm_reciprocal( Hinv , 2.0*M_PI*nx*ax , 2.0*M_PI*ny*ay , 2.0*M_PI*nz*az );
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
        }
      }
    }
  }

  // LAMMPS PPPM::init (ik, no slab, no tip4p) followed by setup.
  // qsqsum = sum of q^2, qsum = sum of q (elementary charges). Estimates use LAMMPS metal units (eV, ang).
  inline void pppm_init_parameters( double g_ewald, double radius, double accuracy_relative, long order, const long mesh[3],
                                    const Mat3d& H, uint64_t natoms, double qsqsum, double qsum, PPPMParameters& p )
  {
    if( order < 2 || order > pppm_constants::MAXORDER )
    {
      ::onika::fatal_error() << "PPPM order must be in [2,"<<pppm_constants::MAXORDER<<"], got "<<order << std::endl;
    }

    p.accuracy_relative = accuracy_relative;
    p.g_ewald_user = g_ewald;
    p.radius = radius;
    p.order = order;
    for( int i = 0 ; i < 3 ; i++ ) p.mesh_user[i] = mesh[i];
    p.qsum = qsum;
    p.qsqsum = qsqsum;

    const bool triclinic = ! is_diagonal( H );
    PPPMEstimator e;
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

    pppm_set_grid_global( e , triclinic , mesh , p.nx , p.ny , p.nz , accuracy );
    if( p.nx < order || p.ny < order || p.nz < order )
    {
      ::onika::fatal_error() << "PPPM mesh "<<p.nx<<"x"<<p.ny<<"x"<<p.nz<<" must have at least order="<<order<<" points in each direction" << std::endl;
    }

    if( ! gewaldflag )
    {
      if( ! e.adjust_gewald() ) ::onika::fatal_error() << "PPPM : could not compute g_ewald" << std::endl;
    }
    p.g_ewald = e.g_ewald;

    // final accuracy estimate (no coulomb table)
    const double dfk = e.df_kspace();
    const double dfr = e.df_rspace();
    p.estimated_accuracy = std::sqrt( dfk*dfk + dfr*dfr );

    pppm_compute_gf_denom( p );
    pppm_compute_rho_coeff( p );
    pppm_setup( p , H );
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

    ReadOnlyPPPMParameters() = default;
    inline ReadOnlyPPPMParameters( const PPPMParameters& p, const Vec3d& bounds_min, const Vec3d& bounds_size )
      : order( p.order ), nlower( p.nlower() ), nx( p.nx ), ny( p.ny ), nz( p.nz )
      , shift( pppm_constants::OFFSET + ( (p.order % 2) ? 0.5 : 0.0 ) )
      , shiftone( (p.order % 2) ? 0.0 : 0.5 )
      , g_ewald( p.g_ewald ), qsum( p.qsum ), volume( p.volume )
      , delvolinv( double(p.nfft()) / p.volume )
      , bmin( bounds_min )
      , delinv( Vec3d{ p.nx / bounds_size.x , p.ny / bounds_size.y , p.nz / bounds_size.z } )
    {
      for( int l = 0 ; l < pppm_constants::MAXORDER ; l++ )
        for( int k = 0 ; k < pppm_constants::MAXORDER ; k++ ) rho_coeff[l][k] = p.rho_coeff[l][k];
    }
  };

  // stencil of one particle : lower mesh index (wrapped into [0,n)) and 1D weights in each direction
  struct PPPMStencil
  {
    int ix, iy, iz;
    double wx[pppm_constants::MAXORDER];
    double wy[pppm_constants::MAXORDER];
    double wz[pppm_constants::MAXORDER];
  };

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
    st.ix = pppm_wrap( gx + p.nlower , p.nx );
    st.iy = pppm_wrap( gy + p.nlower , p.ny );
    st.iz = pppm_wrap( gz + p.nlower , p.nz );
  }

  // self energy + neutralizing background energy of one particle with charge q (LAMMPS PPPM per atom correction), internal units
  ONIKA_HOST_DEVICE_FUNC inline double pppm_self_energy( const ReadOnlyPPPMParameters& p, double q )
  {
    return - COULOMB_CONSTANT * ( p.g_ewald / sqrt(M_PI) * q * q + 0.5 * M_PI * q * p.qsum / ( p.g_ewald * p.g_ewald * p.volume ) );
  }

}
}
