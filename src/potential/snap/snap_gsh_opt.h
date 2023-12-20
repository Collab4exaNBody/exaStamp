#pragma once

#include <cmath> // sqrt.
#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>
#include <onika/integral_constant.h>
#include <onika/declare_if.h>
#include "snap_3Dtypes.h"
#include "snap_constants.h"

namespace SnapExt
{


ONIKA_HOST_DEVICE_FUNC inline constexpr int snap_gsh_idx(int j, int m1, int m2) 
{	//returns index in gsh vector corresponding to (j,m1,m2)
    return j*(j+1)*(2*j+1)/6 + (m2+j)/2*(j+1) + (m1+j)/2 + 1; //sum of j first square integers+offset
}

//#define DECLARE_IF_CONSTEXPR(c,t,v) std::conditional_t<c,t,onika::BoolConst<false> > v; if constexpr (!c) if(v) (void)0

template<class GSHArray, class DGSHArray, class TwoJMaxT, bool ComputeGSH=true , bool ComputeDGSH=true>
ONIKA_HOST_DEVICE_FUNC static inline SnapExt::SinCosTheta snap_compute_gsh_ext(
  double3d const &ran,
  double r,
  double rcut,
  double rfac0,
  double rmin0,
  TwoJMaxT m_two_jmax,
  GSHArray m_gsh,
  DGSHArray m_dgsh,
  onika::BoolConst<ComputeGSH> = {},
  onika::BoolConst<ComputeDGSH> = {} )
{
  using namespace std;
  // cout.precision(12);
  // Initialize constants.

  constexpr bool jmax3 = std::is_same_v< TwoJMaxT , onika::IntConst<6> >;
  constexpr double gsh_sqrt_frac_jmax3[7][7] = {
   { 0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0 }
  ,{ 0x0p+0,0x1p+0,0x1.6a09e667f3bcdp-1,0x1.279a74590331cp-1,0x1p-1,0x1.c9f25c5bfedd9p-2,0x1.a20bd700c2c3ep-2 }
  ,{ 0x0p+0,0x1.6a09e667f3bcdp+0,0x1p+0,0x1.a20bd700c2c3ep-1,0x1.6a09e667f3bcdp-1,0x1.43d136248490fp-1,0x1.279a74590331cp-1 }
  ,{ 0x0p+0,0x1.bb67ae8584caap+0,0x1.3988e1409212ep+0,0x1p+0,0x1.bb67ae8584caap-1,0x1.8c97ef43f7248p-1,0x1.6a09e667f3bcdp-1 }
  ,{ 0x0p+0,0x1p+1,0x1.6a09e667f3bcdp+0,0x1.279a74590331cp+0,0x1p+0,0x1.c9f25c5bfedd9p-1,0x1.a20bd700c2c3ep-1 }
  ,{ 0x0p+0,0x1.1e3779b97f4a8p+1,0x1.94c583ada5b53p+0,0x1.4a7e9cb8a3491p+0,0x1.1e3779b97f4a8p+0,0x1p+0,0x1.d363d1848dcbfp-1 }
  ,{ 0x0p+0,0x1.3988e1409212ep+1,0x1.bb67ae8584caap+0,0x1.6a09e667f3bcdp+0,0x1.3988e1409212ep+0,0x1.186f174f88472p+0,0x1p+0 }
  };

  double const x = ran.x, y = ran.y, z = ran.z;
  //double const r = sqrt(x*x + y*y + z*z);
  //double3d     dr;
  //double const PI=3.14159265359; // For Fortran exact comparison

//We could do as Lammps sna.cpp l.334
//theta0=(r-rmin0)*rfac0*PI/(rcut-rmin0)
//then
//r0=r/theta0
//Strictly equivalent
//These parameters should be an input of the class in the future

  //constexpr double rmin0 = 0.;
  //constexpr double rfac0 = 1.0;  
  double const inv_r0 = (rfac0*M_PI*(r-rmin0)) / (r*(rcut-rmin0)) ;

  double const theta0 = r * inv_r0; // Bartok PhD (2.57)
  double sin_theta0 = 0.0 , cos_theta0 = 0.0; 
  sincos(theta0,&sin_theta0,&cos_theta0);
/*
  float sin_theta0f, cos_theta0f; onika_fast_sincosf(theta0,&sin_theta0f,&cos_theta0f);//__sincosf(theta0,&sin_theta0f,&cos_theta0f);
  const double sin_theta0 = sin_theta0f;
  const double cos_theta0 = cos_theta0f;  
*/
  /*const double sin_theta0 = sin(theta0);
  const double cos_theta0 = cos(theta0);*/
    
  const double inv_tan_theta0 = cos_theta0 / sin_theta0; 
  double const l0 = r / sin_theta0;
  double const z0 = l0 * cos_theta0;
  double const l0i = sin_theta0 / r;
//  double const dz=1./tan_theta0-theta0 / pow(sin_theta0,2);
  double const dz = inv_tan_theta0 - theta0 / (sin_theta0*sin_theta0);
  double const dil0 = (cos_theta0 * inv_r0 - l0i) / r;
  const Complexd imag(0., 1.);
//
  const double3d dr = ran / r;

  // Initialize generalised spherical harmonics.

  // if (m_jmax <= 0) {cerr << "Error: bad jmax" << endl; return 1;}

  // avoid any potential unitialized value
  constexpr Complexd zeroc(0.0,0.0);
  if constexpr (ComputeGSH)
  {
    constexpr Complexd u000(1., 0.); // j = 0, m1 = m2 = 0 - Bartok PhD (B.9)
    m_gsh[0] = zeroc;
    //m_gsh[1] = zeroc;
    m_gsh[snap_gsh_idx(0, 0, 0)] = u000;
  }
  
  if constexpr (ComputeDGSH)
  {
    constexpr complex3d zeroc3(zeroc,zeroc,zeroc);
    constexpr Complexd du000(0., 0.);
    m_dgsh[0] = zeroc3;
    //m_dgsh[1] = zeroc3;
    m_dgsh[snap_gsh_idx(0, 0, 0)]= du000;
  }

  Complexd z_minus = (z0-imag*z)*l0i;  // u
  Complexd z_plus  = (z0+imag*z)*l0i;  // u*
  Complexd x_minus =(x-imag*y)*l0i;    // v*eiphi
  Complexd x_plus  =(x+imag*y)*l0i;    // v*e-iphi

  DECLARE_IF_CONSTEXPR(ComputeDGSH,complex3d,dx_plus);
  DECLARE_IF_CONSTEXPR(ComputeDGSH,complex3d,dx_minus);
  DECLARE_IF_CONSTEXPR(ComputeDGSH,complex3d,dz_plus);
  DECLARE_IF_CONSTEXPR(ComputeDGSH,complex3d,dz_minus);

  if constexpr (ComputeDGSH)
  {
    dx_plus  = (x+imag*y)*dil0*dr;
    dx_minus = (x-imag*y)*dil0*dr;
    dz_plus  = ((z0+imag*z)*dil0+dz*l0i)*dr;
    dz_minus = ((z0-imag*z)*dil0+dz*l0i)*dr;
    dx_plus.x  += l0i;
    dx_plus.y  += imag*l0i;
    dx_minus.x += l0i;
    dx_minus.y += -imag*l0i;
    dz_plus.z  += imag*l0i;
    dz_minus.z +=-imag*l0i;
  }

  // Compute generalised spherical harmonics.

  for (int j = 1 /*0 OK: recursion initialized*/; j <= m_two_jmax; ++j) 
  {
  
    {
      const int m2 = -j;
      for (int m1 = -j; m1 <= j; m1 += 2) 
      {
        // Compute gsh for j, m1, m2.
        const int idx_j_m1_m2 = snap_gsh_idx(j,m1,m2);
        const int idx_jm1_m1p1_m2p1 = snap_gsh_idx(j-1,m1+1,m2+1);
        const int idx_jm1_m1m1_m2p1 = snap_gsh_idx(j-1,m1-1,m2+1);
        double first, second;
        if constexpr ( jmax3 ) { first  = gsh_sqrt_frac_jmax3[(j-m1)/2][(j-m2)/2]; second = gsh_sqrt_frac_jmax3[(j+m1)/2][(j-m2)/2]; }
        if constexpr (!jmax3 ) { first  = sqrt((double)(j-m1)/(double)(j-m2));     second = sqrt((double)(j+m1)/(double)(j-m2)); }
        if constexpr (ComputeGSH)
        {
          m_gsh[idx_j_m1_m2] =       first*z_plus *m_gsh[idx_jm1_m1p1_m2p1]
                             - imag*second*x_minus*m_gsh[idx_jm1_m1m1_m2p1];
        }     
        if constexpr (ComputeDGSH)
        {
          m_dgsh[idx_j_m1_m2] =      first *(dz_plus *m_gsh[idx_jm1_m1p1_m2p1]+z_plus *m_dgsh[idx_jm1_m1p1_m2p1])
			                        - imag*second*(dx_minus*m_gsh[idx_jm1_m1m1_m2p1]+x_minus*m_dgsh[idx_jm1_m1m1_m2p1]);
        }
      }
    }
  
    for (int m2 = 2-j; m2 <= j; m2 += 2) 
    {
      for (int m1 = -j; m1 <= j; m1 += 2) 
      {
        // Compute gsh for j, m1, m2.
        const int idx_j_m1_m2 = snap_gsh_idx(j,m1,m2);
        const int idx_jm1_m1m1_m2m1 = snap_gsh_idx(j-1,m1-1,m2-1);
        const int idx_jm1_m1p1_m2m1 = snap_gsh_idx(j-1,m1+1,m2-1);
        double first, second;
        if constexpr ( jmax3 ) { first  = gsh_sqrt_frac_jmax3[(j+m1)/2][(j+m2)/2]; second = gsh_sqrt_frac_jmax3[(j-m1)/2][(j+m2)/2]; }
        if constexpr (!jmax3 ) { first  = sqrt((double)(j+m1)/(double)(j+m2));     second = sqrt((double)(j-m1)/(double)(j+m2)); }

        if constexpr (ComputeGSH)
        {
          m_gsh[idx_j_m1_m2] = first*z_minus* m_gsh[idx_jm1_m1m1_m2m1] 
                             - imag*second*x_plus*m_gsh[idx_jm1_m1p1_m2m1];
        }
        if constexpr (ComputeDGSH)
        {
          m_dgsh[idx_j_m1_m2] =       first*(dz_minus*m_gsh[idx_jm1_m1m1_m2m1]+z_minus*m_dgsh[idx_jm1_m1m1_m2m1])
				                      - imag*second*(dx_plus *m_gsh[idx_jm1_m1p1_m2m1]+x_plus *m_dgsh[idx_jm1_m1p1_m2m1]);
				}
      }
    }
    
  }

  return { sin_theta0 , cos_theta0 };
}


}
