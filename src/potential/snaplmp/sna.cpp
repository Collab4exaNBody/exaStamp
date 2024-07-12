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

using namespace std;
using namespace LAMMPS_NS;

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

SnapScratchBuffers::SnapScratchBuffers( const SNA* conf , Memory* mem )
  : config(conf)
  , memory(mem)
{
  nmax = 0;
//  rij = nullptr;
//  inside = nullptr;
  wj = nullptr;
  rcutij = nullptr;
  sinnerij = nullptr;
  dinnerij = nullptr;
  element = nullptr;
  ulist_r_ij = nullptr;
  ulist_i_ij = nullptr;
  
  create_twojmax_arrays();
}

SnapScratchBuffers::~SnapScratchBuffers()
{
//  memory->destroy(rij);
//  memory->destroy(inside);
  memory->destroy(wj);
  memory->destroy(rcutij);
  memory->destroy(sinnerij);
  memory->destroy(dinnerij);
  if (config->chem_flag) memory->destroy(element);
  memory->destroy(ulist_r_ij);
  memory->destroy(ulist_i_ij);
  
  destroy_twojmax_arrays();
}

void SnapScratchBuffers::create_twojmax_arrays()
{
//  int jdimpq = config->twojmax + 2;
  memory->create(ulisttot_r, config->idxu_max * config->nelements, "sna:ulisttot_r");
  memory->create(ulisttot_i, config->idxu_max * config->nelements, "sna:ulisttot_i");
  memory->create(dulist_r, config->idxu_max, 3, "sna:dulist_r");
  memory->create(dulist_i, config->idxu_max, 3, "sna:dulist_i");
  memory->create(zlist_r, config->idxz_max * config->ndoubles, "sna:zlist_r");
  memory->create(zlist_i, config->idxz_max * config->ndoubles, "sna:zlist_i");
  memory->create(blist, config->idxb_max * config->ntriples, "sna:blist");
  memory->create(ylist_r, config->idxu_max * config->nelements, "sna:ylist_r");
  memory->create(ylist_i, config->idxu_max * config->nelements, "sna:ylist_i");
}

void SnapScratchBuffers::destroy_twojmax_arrays()
{
  // scratch buffers
  memory->destroy(ulisttot_r);
  memory->destroy(ulisttot_i);
  memory->destroy(dulist_r);
  memory->destroy(dulist_i);
  memory->destroy(zlist_r);
  memory->destroy(zlist_i);
  memory->destroy(blist);
  memory->destroy(ylist_r);
  memory->destroy(ylist_i);
}


void SnapScratchBuffers::grow_rij(int newnmax)
{
  if (newnmax <= nmax) return;

  nmax = newnmax;

//  memory->destroy(rij);
//  memory->destroy(inside);
  memory->destroy(wj);
  memory->destroy(rcutij);
  memory->destroy(sinnerij);
  memory->destroy(dinnerij);
  if (config->chem_flag) memory->destroy(element);
  memory->destroy(ulist_r_ij);
  memory->destroy(ulist_i_ij);
  
//  memory->create(rij, nmax, 3, "pair:rij");
//  memory->create(inside, nmax, "pair:inside");
  memory->create(wj, nmax, "pair:wj");
  memory->create(rcutij, nmax, "pair:rcutij");
  memory->create(sinnerij, nmax, "pair:sinnerij");
  memory->create(dinnerij, nmax, "pair:dinnerij");
  if (config->chem_flag) memory->create(element, nmax, "sna:element");
  memory->create(ulist_r_ij, nmax, config->idxu_max, "sna:ulist_r_ij");
  memory->create(ulist_i_ij, nmax, config->idxu_max, "sna:ulist_i_ij");
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
        bzero[j] = www;
      else
        bzero[j] = www*(j+1);
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
  memory->create(idxcg_block, jdim, jdim, jdim,
                 "sna:idxcg_block");

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxcg_block[j1][j2][j] = idxcg_count;
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idxcg_max = idxcg_count;

  // index list for uarray
  // need to include both halves

  memory->create(idxu_block, jdim,
                 "sna:idxu_block");

  int idxu_count = 0;

  for (int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  idxu_max = idxu_count;

  // index list for beta and B

  int idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  idxb_max = idxb_count;
//  idxb = new SNA_BINDICES[idxb_max];
  memory->create( idxb, idxb_max, "sna:idxb" );

  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count].j1 = j1;
          idxb[idxb_count].j2 = j2;
          idxb[idxb_count].j = j;
          idxb_count++;
        }

  // reverse index list for beta and b

  memory->create(idxb_block, jdim, jdim, jdim,
                 "sna:idxb_block");
  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j1][j2][j] = idxb_count;
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  idxz_max = idxz_count;
  //idxz = new SNA_ZINDICES[idxz_max];
  memory->create(idxz, idxz_max, "sna:idxz" );

  memory->create(idxz_block, jdim, jdim, jdim,
                 "sna:idxz_block");

  idxz_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxz_block[j1][j2][j] = idxz_count;

        // find right beta[jjb] entry
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            idxz[idxz_count].j1 = j1;
            idxz[idxz_count].j2 = j2;
            idxz[idxz_count].j = j;
            idxz[idxz_count].ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
            idxz[idxz_count].ma2max = (2 * ma - j - (2 * idxz[idxz_count].ma1min - j1) + j2) / 2;
            idxz[idxz_count].na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count].ma1min + 1;
            idxz[idxz_count].mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
            idxz[idxz_count].mb2max = (2 * mb - j - (2 * idxz[idxz_count].mb1min - j1) + j2) / 2;
            idxz[idxz_count].nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count].mb1min + 1;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count].jju = jju;

            idxz_count++;
          }
      }
}

/* ---------------------------------------------------------------------- */

void SNA::init()
{
  init_clebsch_gordan();
  init_rootpqarray();
}

void SNA::create_twojmax_arrays()
{
  int jdimpq = twojmax + 2;

  // config constants
  memory->create(rootpqarray, jdimpq, jdimpq,"sna:rootpqarray");
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
  memory->destroy(idxu_block);
  memory->destroy(idxz);
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
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if (m < 0 || m > j) {
              cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= MIN((j1 + j2 - j) / 2,
                          MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
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

            cglist[idxcg_count] = sum * dcg * sfaccg;
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
          int idxcg_count = idxcg_block[j1][j2][j];
          for (int m1 = 0; m1 <= j1; m1++) {
            aa2 = 2*m1-j1;
            for (int m2 = 0; m2 <= j2; m2++) {
              bb2 = 2*m2-j2;
              double cgtmp = cglist[idxcg_count];
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
      rootpqarray[p][q] = sqrt(static_cast<double>(p)/q);
}

/* ---------------------------------------------------------------------- */

void SNA::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2;
           j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  ndoubles = nelements*nelements;
  ntriples = nelements*nelements*nelements;
  if (chem_flag)
    ncoeff = ncount*ntriples;
  else
    ncoeff = ncount;
}

