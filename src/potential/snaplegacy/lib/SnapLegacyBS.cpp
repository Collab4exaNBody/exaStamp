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

#include <exaStamp/potential/snaplegacy/SnapLegacyBS.h>

#include <iostream>
#include <iomanip> // setw.
#include <cmath> // sqrt.
#include <strings.h>
#include <algorithm> 

#define SNAP_VERBOSE_DBG 1

using namespace std;

static inline complex3d conj(complex3d const &c) {
	complex3d r;
	r.x = conj(c.x);
	r.y = conj(c.y);
	r.z = conj(c.z);
	return r;
}

SnapLegacyBS::SnapLegacyBS(double const jmax /*, double const rcut*/ /*, int myspecy*/, double const *coefs, double const *factor)
  : m_nidx(n_idx_bs(floor(2*jmax)))
  , m_jmax(jmax)
  , m_coefs(coefs)
  , m_factor(factor)
  , m_two_jmax(floor(2*jmax))
  , m_gsh_size(pow(floor(2*jmax)+1,3))

{
  m_bs.resize(m_nidx);

  if( m_coefs == nullptr )
  {
    m_default_coefs.assign( m_nidx , 0.0 );
    m_coefs = m_default_coefs.data();
  }
  
  if( m_factor == nullptr )
  {
    m_default_factors.assign( 1 , 1.0 );
    m_factor = m_default_factors.data();
  }  
}

//Compute the number of bispectrum components
int SnapLegacyBS::n_idx_bs(int m_two_jmax)
{
#ifndef NDEBUG
#ifdef LAMMPS
  cout <<"LAMMPS COMPATIBILITY MODE ENABLED" << endl;
#else
  cout <<"FORTRAN COMPATIBILITY MODE" << endl;
#endif
#endif
  int n_idx=0;
  for (int j1=0; j1 <=m_two_jmax; j1+=1) {
#ifdef LAMMPS
    for (int j2=0; j2 <= j1; j2+=1) { //Thompson-Lammps
      for (int j=abs(j1-j2); j<=min(m_two_jmax,j1+j2); j+=2)
#else
   int j2=j1;
      for (int j=j1-j2; j<=min(m_two_jmax,j1+j2); j+=1)
#endif
        {
          if ((j+j1+j2)%2==1)  continue;
#ifdef LAMMPS
	  if (j<j1) continue;
#endif
    	  n_idx+=1;
        }
#ifdef LAMMPS
   } //Thompson-Lammps
#endif
  }
#ifndef NDEBUG
  cout << "NB of bispectrum components: "<< "\t" << n_idx << endl;
#endif
return n_idx;
}

//Update the list of neighbours
int SnapLegacyBS::set_neighbours(double const *rx, double const *ry, double const *rz, int const * _species, double rcut, size_t N_atom)
{
  m_N_atom=N_atom;
  m_r_vec.resize(m_N_atom);
  m_radius.resize(m_N_atom);
  m_species.resize(m_N_atom); 
  for (int i=0; i<m_N_atom; ++i)
  {
    double r =0.;
    r=sqrt(rx[i]*rx[i] + ry[i]*ry[i] + rz[i]*rz[i]);
    m_r_vec[i].x=rx[i]; m_r_vec[i].y=ry[i]; m_r_vec[i].z=rz[i];
    m_radius[i]=r;
#if !defined(NDEBUG) && defined(SNAP_VERBOSE_DBG)
    if (m_radius[i]<rcut) { cout << "neighbour: " << i << ", x,y,z,r : " << rx[i] << ", " << ry[i] << ", " << rz[i] << ", " << m_radius[i] << endl; }
#endif

    if( _species != nullptr ) { m_species[i] = _species[i]; }
    else { m_species[i] = 0; }
  }
  return 0;
}

//Compute the C_{mm] coefficients (linear combination of GSH)
int SnapLegacyBS::compute_cmm( double rcut )
{
  if(m_vec_gsh.size() != size_t(m_N_atom) )
	{
		m_vec_gsh.resize(m_N_atom);
	}
	 for (unsigned int i=0;i<m_vec_gsh.size();i++)
		m_vec_gsh[i].set_jmax(m_jmax);
  double fcut;
  double dfcut;
  //Resize m_cmm and dcmm if needed
  if( ssize_t(m_cmm.size()) < m_gsh_size)
    m_cmm.resize(m_gsh_size);

  if( ssize_t(m_dcmm.size()) < m_N_atom*m_gsh_size)
    m_dcmm.resize(m_N_atom*m_gsh_size);

  //Initialization
  for (int j = 0; j <= m_two_jmax; ++j) {
    for (int m2 = -j; m2 <= j; m2 += 2) {
      for (int m1 = -j; m1 <= j; m1 += 2) {
	m_cmm[idx(j,m1,m2)]=0.;
	for (int i = 0; i<m_N_atom; ++i) {
       	  m_dcmm[didx(i,j,m1,m2)]=0.;
        }
      }
    }
  }

  for (int j = 0; j <= m_two_jmax; ++j) {
    for (int m1 = -j; m1 <= j; m1 += 2) {
      m_cmm[idx(j,m1,m1)]=1.;
    }
  }
  
 double const PI=3.14159265359;
 double const rmin0=0.;
 m_N_neigh=0;
 for (int i=0; i<m_N_atom; ++i)
 { //loop on neighbours
   if ( (m_radius[i]>rcut) || (m_radius[i]<1.e-12) ) 
     {
#if !defined(NDEBUG) && defined(SNAP_VERBOSE_DBG)
      cout << i << " th atom ignored" << endl;
#endif
      continue;}
   m_N_neigh+=1;
   m_vec_gsh[i].compute_gsh(m_r_vec[i],rcut); //Compute the gsh for atom i
#ifdef LAMMPS
   fcut=0.5*(cos(PI*(m_radius[i]-rmin0)/(rcut-rmin0))+1.)*m_factor[m_species[i]];
   dfcut=-0.5*PI/(rcut-rmin0)*sin(PI*(m_radius[i]-rmin0)/(rcut-rmin0))*m_factor[m_species[i]];
#else    
   fcut=0.5*(cos(PI*m_radius[i]/rcut)+1.)*m_factor[m_species[i]];
   dfcut=-0.5*PI/rcut*sin(PI*m_radius[i]/rcut)*m_factor[m_species[i]];
#endif
   for (int j = 0; j <= m_two_jmax; ++j)
   {
     for (int m2 = -j; m2 <= j; m2 += 2)
     {
       for (int m1 = -j; m1 <= j; m1 += 2)
       {
	       m_cmm[idx(j,m1,m2)] += m_vec_gsh[i].gsh_val(j,m1,m2)*fcut;
	       m_dcmm[didx(i,j,m1,m2)] =dfcut*m_r_vec[i]*m_vec_gsh[i].gsh_val(j,m1,m2)/m_radius[i]+fcut*m_vec_gsh[i].dgsh_val(j,m1,m2);
       }
     }
   }
 }
#if !defined(NDEBUG) && defined(SNAP_VERBOSE_DBG)
 cout << "Number of neighbours: " << "\t" << m_N_neigh << endl;
#endif
 return 0;
}

int SnapLegacyBS::compute_bs(int myspecy, double rcut, const SnapLegacyCG &cg)
{
  m_dbs.resize(m_N_atom*m_nidx);
  size_t bs_idx=0;
  size_t bs_didx=0;
  int m21;
  int m22;
  complex<double> lbs; //local bispectrum
  complex3d ldbs; //local bispectrum derivatives


  //initialize bs and dbs to zero
  m_bs.assign(m_nidx,0.);
  m_dbs.assign(m_nidx*m_N_atom,0.);
  bs_didx=-1;
  for (int i=0; i < m_N_atom ; i++)
  {
    bs_idx=-1;
    for (int j1=0; j1 <= m_two_jmax; j1++)
    {
#ifdef LAMMPS
      // Thompson - Lammps
      for (int j2=0; j2 <= j1; j2++)
      { 
        for (int j=abs(j1-j2); j<=min(m_two_jmax,j1+j2); j+=2)
#else
        int j2=j1; //
        for (int j=j1-j2; j<=min(m_two_jmax,j1+j2); j++)
#endif
        {

          if ((j+j1+j2)%2==1) continue;
#ifdef LAMMPS
          if (j<j1) continue; // Thompson
#endif
	        bs_idx+=1;
	        bs_didx+=1;
          if ( (m_radius[i]<1.e-12) || (m_radius[i]>rcut) ) continue;
	        for (int m1=-j; m1 <=j; m1+=2)
	        {
	          for (int m2=-j; m2 <=j; m2+=2)
	          {
              lbs=0.;
              ldbs=0.;
              for (int m11=max(-j1,m1-j2); m11<=min(j1,m1+j2); m11+=2)
              {
                for (int m12=max(-j1,m2-j2); m12<=min(j1,m2+j2); m12+=2)
                {
                  m21=m1-m11; m22=m2-m12;

                  lbs+=cg.val(j,j1,j2,m1,m11,m21) * cg.val(j,j1,j2,m2,m12,m22) * m_cmm[idx(j1,m11,m12)] * m_cmm[idx(j2,m21,m22)];
                  ldbs+=cg.val(j,j1,j2,m1,m11,m21) * cg.val(j,j1,j2,m2,m12,m22) * (m_dcmm[didx(i,j1,m11,m12)]*m_cmm[idx(j2,m21,m22)] + m_cmm[idx(j1,m11,m12)]* m_dcmm[didx(i,j2,m21,m22)]);
                }
              }
              m_bs[bs_idx]+=conj(m_cmm[idx(j,m1,m2)])*lbs;
              m_dbs[bs_didx]+=-(conj(m_dcmm[didx(i,j,m1,m2)]) * lbs + conj(m_cmm[idx(j,m1,m2)])*ldbs);
            }   //m2
          }	//m1
        }	//j
      }		//j2
#ifdef LAMMPS
    }		//j1  // Thompson
#endif
  }		//atoms
  
#if !defined(NDEBUG) && defined(SNAP_VERBOSE_DBG)
  cout << "bs: values:" << endl;
#endif
  bs_idx=0;
  double R_neigh=m_N_neigh;
  if (R_neigh<0.9) return 0;
  for (int j1=0; j1 <= m_two_jmax; j1++) {
#ifdef LAMMPS
    for (int j2=0; j2 <= j1; j2++) { // Thompson - Lammps
#else
    int j2=j1;
#endif

      for (int j=j1-j2; j<=min(m_two_jmax,j1+j2); j++) {
        if ((j+j1+j2)%2==1) continue;
#ifdef LAMMPS
	if (j<j1) continue;
#endif
      m_bs[bs_idx]=m_bs[bs_idx]/R_neigh*m_factor[myspecy];
#if !defined(NDEBUG) && defined(SNAP_VERBOSE_DBG)
      cout << bs_idx << "\t" << m_bs[bs_idx] << endl;
#endif
      bs_idx+=1;
      }
#ifdef LAMMPS
    }
#endif
  }
/* To display derivatives values, uncoment.
  cout << "dbs: values:" << endl;
  bs_didx=-1;
  for (int i=0; i<m_N_atom; ++i) {
  bs_idx=-1;
    cout << "=====Atom " << "\t" << i << "\t" << "=====" << endl;
    for (int j1=0; j1 <= m_two_jmax; j1++) {
#ifdef LAMMPS
      for (int j2=0; j2 <= j1; j2++) { // Thompson - Lammps
#else
      int j2=j1;
#endif
        for (int j=j1-j2; j<=min(m_two_jmax,j1+j2); j++) {
          if ((j+j1+j2)%2==1) continue;
#ifdef LAMMPS
	  if (j<j1) continue;
#endif
        bs_didx+=1;
      cout << bs_idx << "\t" << m_dbs[bs_didx].x << endl;
      cout << bs_idx << "\t" << m_dbs[bs_didx].y << endl;
      cout << bs_idx << "\t" << m_dbs[bs_didx].z << endl;
        }
#ifdef LAMPPS
      }
#endif
    }
  }
*/
  return 0;
}

complex<double> SnapLegacyBS::en_val(int myspecy)
{
   m_energy=m_coefs[myspecy*m_nidx];
   for (int i=0; i<m_nidx; i++)
   {
     m_energy+=m_bs[i]*m_coefs[myspecy*m_nidx+i+1];
   }
  return m_energy;
}

complex3d SnapLegacyBS::force_val(int i)
{
  int k=m_species[i];
  complex3d f_val;
  f_val=0.;
  for (int j=0; j< m_nidx; j++)
  {
     f_val+=m_dbs[i*m_nidx+j]*m_coefs[k*m_nidx+j+1];
  }
  return f_val;
}

