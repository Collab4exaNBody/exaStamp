#include <exaStamp/potential/snaplegacy/SnapLegacyGSH.h>
#include <iostream>
#include <iomanip> // setw.
#include <cmath> // sqrt.
                                              
//Member functions
//
//Public member functions 
//
//Computes all the gsh and dgsh coefficients


int SnapLegacyGSH::compute_gsh(double3d const &ran, double rcut) {

  cout.precision(12);
  // Initialize constants.

  double const x = ran.x, y = ran.y, z = ran.z;
  double const r = sqrt(x*x + y*y + z*z);
  double3d     dr;
  double const PI=3.14159265359; // For Fortran exact comparison
#ifdef LAMMPS
//We could do as Lammps sna.cpp l.334
//theta0=(r-rmin0)*rfac0*PI/(rcut-rmin0)
//then
//r0=r/theta0
//Strictly equivalent
//These parameters should be an input of the class in the future
  double const rmin0 = 0.;
  //  double const rfac0 = 1.0;
  double const rfac0 = 0.99363;  
  double const r0 = (r*(rcut-rmin0))/(rfac0*PI*(r-rmin0));
#else
  double const r0 = rcut/(PI-0.020);                             // Bartok PhD (2.57)
#endif
  double const theta0 = r/r0; // Bartok PhD (2.57)
  double const l0 = r/sin(theta0);
  double const z0 = l0*cos(theta0);
  double const l0i = sin(theta0)/r;
  double const dz=1./tan(theta0)-theta0/pow(sin(theta0),2);
  double const dil0=(cos(theta0)/r0 - l0i)/r;
  complex<double> imag(0., 1.);
//
  dr=ran/r;

  // Initialize generalised spherical harmonics.

  if (m_jmax <= 0) {cerr << "Error: bad jmax" << endl; return 1;}
  m_gsh.resize(size());
  m_dgsh.resize(size());

  complex<double> u000(1., 0.); // j = 0, m1 = m2 = 0 - Bartok PhD (B.9)
  complex<double> du000(0., 0.);
  m_gsh[idx(0, 0, 0)] = u000;
  m_dgsh[idx(0, 0, 0)]= du000;

  complex<double> z_minus = (z0-imag*z)*l0i;  // u
  complex<double> z_plus  = (z0+imag*z)*l0i;  // u*
  complex<double> x_minus =(x-imag*y)*l0i;    // v*eiphi
  complex<double> x_plus  =(x+imag*y)*l0i;    // v*e-iphi
  complex3d       dx_plus ; 
  complex3d       dx_minus;
  complex3d       dz_plus ;
  complex3d       dz_minus;

  for (int i=0; i<=2; i++) {
     dx_plus =(x+imag*y)*dil0*dr;
     dx_minus=(x-imag*y)*dil0*dr;
     dz_plus =((z0+imag*z)*dil0+dz*l0i)*dr;
     dz_minus =((z0-imag*z)*dil0+dz*l0i)*dr;
  }
  dx_plus.x += l0i;
  dx_plus.y +=imag*l0i;
  dx_minus.x += l0i;
  dx_minus.y +=-imag*l0i;

  dz_plus.z +=imag*l0i;
  dz_minus.z +=-imag*l0i;


  // Compute generalised spherical harmonics.

  int two_jmax = floor(2*m_jmax);

  for (int j = 1; j <= two_jmax; ++j) {
    for (int m1=-j; m1<=j ; m1 +=2) {
	m_gsh[idx(j,m1,m1)]=u000;
    }
  }

  for (int j = 1 /*0 OK: recursion initialized*/; j <= two_jmax; ++j) {
    for (int m2 = -j; m2 <= j; m2 += 2) {
      for (int m1 = -j; m1 <= j; m1 += 2) {
        // Compute gsh for j, m1, m2.

#ifdef USE_B7
        if ( m2 == j) {// Bartok PhD (B.7)
          complex<double>  first = sqrt((double)(j+m1)/(double)(j+m2));
          complex<double> second = sqrt((double)(j-m1)/(double)(j+m2));

          m_gsh[idx(j, m1, m2)] = first*z_minus* m_gsh[idx(j-1, m1-1, m2-1)] - imag*second*x_plus*m_gsh[idx(j-1, m1+1, m2-1)];

	  m_dgsh[idx(j, m1, m2)] =  first*(dz_minus*m_gsh[idx(j-1,m1-1,m2-1)]+z_minus*m_dgsh[idx(j-1,m1-1,m2-1)])-
					imag*second*(dx_plus*m_gsh[idx(j-1,m1+1,m2-1)]+x_plus*m_dgsh[idx(j-1,m1+1,m2-1)]);
        }
	else { // Bartok PhD B.7
          complex<double>  first = sqrt((double)(j-m1)/(double)(j-m2));
          complex<double> second = sqrt((double)(j+m1)/(double)(j-m2));

          m_gsh[idx(j, m1, m2)] = first*z_plus      m_*gsh[idx(j-1, m1+1, m2+1)] - imag*second*x_minus*m_gsh[idx(j-1, m1-1, m2+1)];

	  m_dgsh[idx(j, m1, m2)] = first*(dz_plus*m_gsh[idx(j-1,m1+1,m2+1)]+z_plus*m_dgsh[idx(j-1,m1+1,m2+1)])-
				     imag*second*(dx_minus*m_gsh[idx(j-1,m1-1,m2+1)]+x_minus*m_dgsh[idx(j-1,m1-1,m2+1)]);
        }
#else        
        if (m2 == -j) { // Bartok PhD (B.8)
          complex<double>  first = sqrt((double)(j-m1)/(double)(j-m2));
          complex<double> second = sqrt((double)(j+m1)/(double)(j-m2));

          m_gsh[idx(j, m1, m2)] = first*z_plus      *m_gsh[idx(j-1, m1+1, m2+1)] - imag*second*x_minus*m_gsh[idx(j-1, m1-1, m2+1)];

	  m_dgsh[idx(j, m1, m2)] = first*(dz_plus*m_gsh[idx(j-1,m1+1,m2+1)]+z_plus*m_dgsh[idx(j-1,m1+1,m2+1)])-
				     imag*second*(dx_minus*m_gsh[idx(j-1,m1-1,m2+1)]+x_minus*m_dgsh[idx(j-1,m1-1,m2+1)]);
        }
        else { // Bartok PhD (B.8)
          complex<double>  first  = sqrt((double)(j+m1)/(double)(j+m2));
          complex<double> second  = sqrt((double)(j-m1)/(double)(j+m2));

          m_gsh[idx(j, m1, m2)] = first*z_minus* m_gsh[idx(j-1, m1-1, m2-1)] - imag*second*x_plus*m_gsh[idx(j-1, m1+1, m2-1)];

	  m_dgsh[idx(j, m1, m2)] =  first*(dz_minus*m_gsh[idx(j-1,m1-1,m2-1)]+z_minus*m_dgsh[idx(j-1,m1-1,m2-1)])-
					imag*second*(dx_plus*m_gsh[idx(j-1,m1+1,m2-1)]+x_plus*m_dgsh[idx(j-1,m1+1,m2-1)]);
        }
#endif
      }
    }
  }
  return 0;
}

int SnapLegacyGSH::size() {
      int two_jmax = floor(2*m_jmax);
      size_t nb_j = two_jmax + 1; // Account for j = 0.
      size_t nb_m_per_j = two_jmax + 1; // -j <= m <= j with cnt 2.
      return nb_j*nb_m_per_j*nb_m_per_j; // Upper bound.
    };
