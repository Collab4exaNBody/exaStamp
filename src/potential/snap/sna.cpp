// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

#include "sna.h"

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace SnapInternal;

/* ----------------------------------------------------------------------

   this implementation is based on the method outlined
   in Bartok[1], using formulae from VMK[2].

   for the Clebsch-Gordan coefficients, we
   convert the VMK half-integral labels
   a, b, c, alpha, beta, gamma
   to array offsets j1, j2, j, m1, m2, m
   using the following relations:

   j1 = 2*a
   j2 = 2*b
   j =  2*c

   m1 = alpha+a      2*alpha = 2*m1 - j1
   m2 = beta+b    or 2*beta = 2*m2 - j2
   m =  gamma+c      2*gamma = 2*m - j

   in this way:

   -a <= alpha <= a
   -b <= beta <= b
   -c <= gamma <= c

   becomes:

   0 <= m1 <= j1
   0 <= m2 <= j2
   0 <= m <= j

   and the requirement that
   a+b+c be integral implies that
   j1+j2+j must be even.
   The requirement that:

   gamma = alpha+beta

   becomes:

   2*m - j = 2*m1 - j1 + 2*m2 - j2

   Similarly, for the Wigner U-functions U(J,m,m') we
   convert the half-integral labels J,m,m' to
   array offsets j,ma,mb:

   j = 2*J
   ma = J+m
   mb = J+m'

   so that:

   0 <= j <= 2*Jmax
   0 <= ma, mb <= j.

   For the bispectrum components B(J1,J2,J) we convert to:

   j1 = 2*J1
   j2 = 2*J2
   j = 2*J

   and the requirement:

   |J1-J2| <= J <= J1+J2, for j1+j2+j integral

   becomes:

   |j1-j2| <= j <= j1+j2, for j1+j2+j even integer

   or

   j = |j1-j2|, |j1-j2|+2,...,j1+j2-2,j1+j2

   [1] Albert Bartok-Partay, "Gaussian Approximation..."
   Doctoral Thesis, Cambridge University, (2009)

   [2] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii,
   "Quantum Theory of Angular Momentum," World Scientific (1988)

------------------------------------------------------------------------- */

static inline double factorial(int n)
{
  double r = 1.0;
  for(;n>=2;--n) r *= n;
  return r;
}

SNA::SNA(Memory* mem, double rfac0_in, int twojmax_in,
         double rmin0_in, int switch_flag_in, int bzero_flag_in,
         int chem_flag_in, int bnorm_flag_in, int wselfall_flag_in,
         int nelements_in, int switch_inner_flag_in)
{
  memory = mem;

  wself = 1.0;

  rfac0 = rfac0_in;
  rmin0 = rmin0_in;
  switch_flag = switch_flag_in;
  switch_inner_flag = switch_inner_flag_in;
  bzero_flag = bzero_flag_in;
  chem_flag = chem_flag_in;
  bnorm_flag = bnorm_flag_in;
  wselfall_flag = wselfall_flag_in;

  if (bnorm_flag != chem_flag)
    std::cerr<< "bnormflag and chemflag are not equal, This is probably not what you intended" << std::endl;

  if (chem_flag)
    nelements = nelements_in;
  else
    nelements = 1;

  twojmax = twojmax_in;

  compute_ncoeff();

  build_indexlist();
  create_twojmax_arrays();

  if (bzero_flag) {
    double www = wself*wself*wself;
    for (int j = 0; j <= twojmax; j++)
      if (bnorm_flag)
        BZERO(j) = www;
      else
        BZERO(j) = www*(j+1);
  }

}

/* ---------------------------------------------------------------------- */

SNA::~SNA()
{  
  delete[] idxz;
  delete[] idxb;
  destroy_twojmax_arrays();
}

void SNA::build_indexlist()
{

  // index list for cglist

  int jdim = twojmax + 1;
  memory->create( idxcg_block, jdim*jdim*jdim, "sna:idxcg_block" );

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        IDXCG_BLOCK(j1,j2,j) = idxcg_count;
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idxcg_max = idxcg_count;

  // index list for uarray
  // need to include both halves

//  memory->create(idxu_block, jdim, "sna:idxu_block");

  int idxu_count = 0;
  for (int j = 0; j <= twojmax; j++) {
//    IDXU_BLOCK(j) = idxu_count;
    if( idxu_count != IDXU_BLOCK(j) ) exanb::fatal_error()<<"bad idxu_count: idxu_count="<<idxu_count<<" , IDXU_BLOCK(j)="<<IDXU_BLOCK(j) <<std::endl;
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  idxu_max = idxu_count;
  if( idxu_max != SUM_INT_SQR(twojmax+1) ) exanb::fatal_error() << "idxu_max is not the sum of twojmax+1 first squares as expected"<<std::endl;

  // index list for beta and B

  int idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  idxb_max = idxb_count;
//  idxb = new SNA_BINDICES[idxb_max];
  memory->create( idxb, idxb_max, "sna:idxb" );

  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          IDXB(idxb_count).j1 = j1;
          IDXB(idxb_count).j2 = j2;
          IDXB(idxb_count).j = j;
          idxb_count++;
        }

  // reverse index list for beta and b

  memory->create( idxb_block, jdim*jdim*jdim, "sna:idxb_block" );
  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          IDXB_BLOCK(j1,j2,j) = idxb_count;
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  idxz_max = idxz_count;
  //idxz = new SNA_ZINDICES[idxz_max];
  memory->create(idxz, idxz_max, "sna:idxz" );

  memory->create( idxz_block, jdim*jdim*jdim, "sna:idxz_block" );

  idxz_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        IDXZ_BLOCK(j1,j2,j) = idxz_count;

        // find right beta[jjb] entry
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            IDXZ(idxz_count).j1 = j1;
            IDXZ(idxz_count).j2 = j2;
            IDXZ(idxz_count).j = j;
            IDXZ(idxz_count).ma1min = std::max(0, (2 * ma - j - j2 + j1) / 2);
            IDXZ(idxz_count).ma2max = (2 * ma - j - (2 * IDXZ(idxz_count).ma1min - j1) + j2) / 2;
            IDXZ(idxz_count).na = std::min(j1, (2 * ma - j + j2 + j1) / 2) - IDXZ(idxz_count).ma1min + 1;
            IDXZ(idxz_count).mb1min = std::max(0, (2 * mb - j - j2 + j1) / 2);
            IDXZ(idxz_count).mb2max = (2 * mb - j - (2 * IDXZ(idxz_count).mb1min - j1) + j2) / 2;
            IDXZ(idxz_count).nb = std::min(j1, (2 * mb - j + j2 + j1) / 2) - IDXZ(idxz_count).mb1min + 1;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = IDXU_BLOCK(j) + (j+1)*mb + ma;
            IDXZ(idxz_count).jju = jju;

            idxz_count++;
          }
      }
      
  assert( idxz_count == idxz_max );

  // create alternative version of ZIndices table for better performances
  idxu_max_alt = exaStamp::snap_force_Y_count(twojmax,idxu_max);
  if( idxu_max_alt < idxu_max )
  {
    memory->create(y_jju_map, idxu_max, "sna:y_jju_map" );
    for(int jju=0;jju<idxu_max;jju++)
    {
      y_jju_map[jju] = exaStamp::snap_force_Y_map(twojmax,jju);
    }
  }

  idxz_max_alt = 0;
  memory->create(idxz_alt, idxz_max, "sna:idxz_alt" );
  for(int jjz=0;jjz<idxz_max;jjz++)
  {
    const int jju = IDXZ(jjz).jju;
    if( y_jju_map == nullptr || y_jju_map[jju] != -1 )
    {
      const int j1 = IDXZ(jjz).j1;
      const int j2 = IDXZ(jjz).j2;
      const int j = IDXZ(jjz).j;
      const int ma1min = IDXZ(jjz).ma1min;
      const int ma2max = IDXZ(jjz).ma2max;
      const int na = IDXZ(jjz).na;
      const int mb1min = IDXZ(jjz).mb1min;
      const int mb2max = IDXZ(jjz).mb2max;
      const int nb = IDXZ(jjz).nb;

      idxz_alt[idxz_max_alt].j1 = j1 ;
      idxz_alt[idxz_max_alt].j2 = j2 ;
      idxz_alt[idxz_max_alt].j = j ;
      idxz_alt[idxz_max_alt].ma1min = ma1min ;
      idxz_alt[idxz_max_alt].ma2max = ma2max ;
      idxz_alt[idxz_max_alt].mb1min = mb1min ;
      idxz_alt[idxz_max_alt].mb2max = mb2max ;
      idxz_alt[idxz_max_alt].na = na ;
      idxz_alt[idxz_max_alt].nb = nb ;
      idxz_alt[idxz_max_alt].jju = y_jju_map[jju];

      //const int idxcg_block_j1_j2_j = IDXCG_BLOCK(j1,j2,j);
      idxz_alt[idxz_max_alt].idx_cgblock = IDXCG_BLOCK(j1,j2,j);
      
      idxz_alt[idxz_max_alt].jjb = -1;
      if (j >= j1) {
        idxz_alt[idxz_max_alt].jjb = IDXB_BLOCK(j1,j2,j);
      } else if (j >= j2) {
        idxz_alt[idxz_max_alt].jjb = IDXB_BLOCK(j,j2,j1);
      } else {
        idxz_alt[idxz_max_alt].jjb = IDXB_BLOCK(j2,j,j1);
      }
      
      ++ idxz_max_alt;
    }
  }
  
  assert( idxz_max_alt <= idxz_max );
  for(int jjz=idxz_max_alt;jjz<idxz_max;jjz++)
  {
    idxz_alt[jjz].jju = -1;
    idxz_alt[jjz].na = 0;
    idxz_alt[jjz].nb = 0;
  }
  
/*
  std::stable_sort( idxz_alt , idxz_alt+idxz_max ,
    [this](const SNA_ZINDICES_ALT& idxz_a , const SNA_ZINDICES_ALT& idxz_b)->bool
    {
      return ( idxz_a.jju < idxz_b.jju );
      if( idxz_a.jju > idxz_b.jju ) return false;
    } );
*/

}

/* ---------------------------------------------------------------------- */

void SNA::init()
{
  init_clebsch_gordan();
  init_rootpqarray();
}

void SNA::create_twojmax_arrays()
{
  //int jdimpq = twojmax + 2;

  // config constants
  memory->create( rootpqarray, twojmax*twojmax, "sna:rootpqarray" );
  memory->create(cglist, idxcg_max, "sna:cglist");
  if (bzero_flag) memory->create(bzero, twojmax+1,"sna:bzero");
  else bzero = nullptr;
}

/* ---------------------------------------------------------------------- */

void SNA::destroy_twojmax_arrays()
{
  // configuration constants
  memory->destroy(rootpqarray);
  memory->destroy(cglist);
  memory->destroy(idxcg_block);
//  memory->destroy(idxu_block);
  memory->destroy(idxz);
  memory->destroy(idxz_alt);
  memory->destroy(idxz_block);
  memory->destroy(idxb);
  memory->destroy(idxb_block);
  
  if (bzero_flag)
    memory->destroy(bzero);

}

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

double SNA::deltacg(int j1, int j2, int j)
{
  double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
  return sqrt(factorial((j1 + j2 - j) / 2) *
              factorial((j1 - j2 + j) / 2) *
              factorial((-j1 + j2 + j) / 2) / sfaccg);
}

/* ----------------------------------------------------------------------
   assign Clebsch-Gordan coefficients using
   the quasi-binomial formula VMK 8.2.1(3)
------------------------------------------------------------------------- */

void SNA::init_clebsch_gordan()
{
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if (m < 0 || m > j) {
              CGLIST(idxcg_count) = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = std::max(0, std::max(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= std::min((j1 + j2 - j) / 2,
                          std::min((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial(z) *
                 factorial((j1 + j2 - j) / 2 - z) *
                 factorial((j1 - aa2) / 2 - z) *
                 factorial((j2 + bb2) / 2 - z) *
                 factorial((j - j2 + aa2) / 2 + z) *
                 factorial((j - j1 - bb2) / 2 + z));
            }

            cc2 = 2 * m - j;
            dcg = deltacg(j1, j2, j);
            sfaccg = sqrt(factorial((j1 + aa2) / 2) *
                          factorial((j1 - aa2) / 2) *
                          factorial((j2 + bb2) / 2) *
                          factorial((j2 - bb2) / 2) *
                          factorial((j  + cc2) / 2) *
                          factorial((j  - cc2) / 2) *
                          (j + 1));

            CGLIST(idxcg_count) = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
}

/* ----------------------------------------------------------------------
   print out values of Clebsch-Gordan coefficients
   format and notation follows VMK Table 8.11
------------------------------------------------------------------------- */

void SNA::print_clebsch_gordan()
{
  //if (comm->me) return;

  int aa2, bb2, cc2;
  for (int j = 0; j <= twojmax; j += 1) {
    printf("c = %g\n",j/2.0);
    printf("a alpha b beta C_{a alpha b beta}^{c alpha+beta}\n");
    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        if (j1-j2 <= j && j1+j2 >= j && (j1+j2+j)%2 == 0) {
          int idxcg_count = IDXCG_BLOCK(j1,j2,j);
          for (int m1 = 0; m1 <= j1; m1++) {
            aa2 = 2*m1-j1;
            for (int m2 = 0; m2 <= j2; m2++) {
              bb2 = 2*m2-j2;
              double cgtmp = CGLIST(idxcg_count);
              cc2 = aa2+bb2;
              if (cc2 >= -j && cc2 <= j)
                if (j1 != j2 || (aa2 > bb2 && aa2 >= -bb2) || (aa2 == bb2 && aa2 >= 0))
                  printf("%4g %4g %4g %4g %10.6g\n",
                         j1/2.0,aa2/2.0,j2/2.0,bb2/2.0,cgtmp);
              idxcg_count++;
            }
          }
        }
  }
}

/* ----------------------------------------------------------------------
   pre-compute table of sqrt[p/m2], p, q = 1,twojmax
   the p = 0, q = 0 entries are allocated and skipped for convenience.
------------------------------------------------------------------------- */

void SNA::init_rootpqarray()
{
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      ROOTPQARRAY(p,q) = sqrt(static_cast<double>(p)/q);
}

/* ---------------------------------------------------------------------- */

void SNA::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2;
           j <= std::min(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  ndoubles = nelements*nelements;
  ntriples = nelements*nelements*nelements;
  if (chem_flag)
    ncoeff = ncount*ntriples;
  else
    ncoeff = ncount;
}

