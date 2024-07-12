#pragma once

#include <cmath>
#include "snap_compute_ui.h"

namespace exaStamp
{
  using namespace exanb;

  static inline double snap_compute_dsfac( // READ ONLY
                                           double rmin0, bool switch_flag, bool switch_inner_flag
                                           // ORIGINAL PARAMETERS
                                         , double r, double rcut, double sinner, double dinner)
  {
    double dsfac, sfac_outer, dsfac_outer, sfac_inner, dsfac_inner;
    if (switch_flag == 0) dsfac_outer = 0.0;
    else if (r <= rmin0) dsfac_outer = 0.0;
    else if (r > rcut) dsfac_outer = 0.0;
    else {
      double rcutfac = M_PI / (rcut - rmin0);
      dsfac_outer = -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
    }

    // some duplicated computation, but rarely visited

    if (switch_inner_flag == 1 && r < sinner + dinner) {
      if (r > sinner - dinner) {

        // calculate sfac_outer

        if (switch_flag == 0) sfac_outer = 1.0;
        else if (r <= rmin0) sfac_outer = 1.0;
        else if (r > rcut) sfac_outer = 0.0;
        else {
	        double rcutfac = M_PI / (rcut - rmin0);
	        sfac_outer = 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
        }

        // calculate sfac_inner

        double rcutfac = (M_PI/2) / dinner;
        sfac_inner = 0.5 * (1.0 - cos( (M_PI/2) + (r - sinner) * rcutfac));
        dsfac_inner = 0.5 * rcutfac * sin( (M_PI/2) + (r - sinner) * rcutfac);
        dsfac = dsfac_outer*sfac_inner + sfac_outer*dsfac_inner;
      } else dsfac = 0.0;
    } else dsfac = dsfac_outer;

    return dsfac;
  }

  /* ----------------------------------------------------------------------
     Compute derivatives of Wigner U-functions for one neighbor
     see comments in compute_uarray()
  ------------------------------------------------------------------------- */
                        
  static inline void snap_compute_duarray( // READ ONLY
                             int twojmax, int const * idxu_block
                           , double const * const * ulist_r_ij, double const * const * ulist_i_ij, double const * const * rootpqarray
                           , double const * sinnerij, double const * dinnerij
                           , double rmin0, bool switch_flag, bool switch_inner_flag
                             // WRITE ONLY
                           , double * const * dulist_r, double * const * dulist_i
                             // ORIGINAL PARAMETERS
                           , double x, double y, double z
                           , double z0, double r, double dz0dr
                           , double wj, double rcut, int jj )
  {
    double r0inv;
    double a_r, a_i, b_r, b_i;
    double da_r[3], da_i[3], db_r[3], db_i[3];
    double dz0[3], dr0inv[3], dr0invdr;
    double rootpq;

    double rinv = 1.0 / r;
    double ux = x * rinv;
    double uy = y * rinv;
    double uz = z * rinv;

    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = z0 * r0inv;
    a_i = -z * r0inv;
    b_r = y * r0inv;
    b_i = -x * r0inv;

    dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

    dr0inv[0] = dr0invdr * ux;
    dr0inv[1] = dr0invdr * uy;
    dr0inv[2] = dr0invdr * uz;

    dz0[0] = dz0dr * ux;
    dz0[1] = dz0dr * uy;
    dz0[2] = dz0dr * uz;

    for (int k = 0; k < 3; k++) {
      da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
      da_i[k] = -z * dr0inv[k];
    }

    da_i[2] += -r0inv;

    for (int k = 0; k < 3; k++) {
      db_r[k] = y * dr0inv[k];
      db_i[k] = -x * dr0inv[k];
    }

    db_i[0] += -r0inv;
    db_r[1] += r0inv;

    double const * ulist_r = ulist_r_ij[jj];
    double const * ulist_i = ulist_i_ij[jj];

    dulist_r[0][0] = 0.0;
    dulist_r[0][1] = 0.0;
    dulist_r[0][2] = 0.0;
    dulist_i[0][0] = 0.0;
    dulist_i[0][1] = 0.0;
    dulist_i[0][2] = 0.0;

    for (int j = 1; j <= twojmax; j++) {
      int jju = idxu_block[j];
      int jjup = idxu_block[j-1];
      for (int mb = 0; 2*mb <= j; mb++) {
        dulist_r[jju][0] = 0.0;
        dulist_r[jju][1] = 0.0;
        dulist_r[jju][2] = 0.0;
        dulist_i[jju][0] = 0.0;
        dulist_i[jju][1] = 0.0;
        dulist_i[jju][2] = 0.0;

        for (int ma = 0; ma < j; ma++) {
          rootpq = rootpqarray[j - ma][j - mb];
          for (int k = 0; k < 3; k++) {
            dulist_r[jju][k] +=
              rootpq * (da_r[k] * ulist_r[jjup] +
                        da_i[k] * ulist_i[jjup] +
                        a_r * dulist_r[jjup][k] +
                        a_i * dulist_i[jjup][k]);
            dulist_i[jju][k] +=
              rootpq * (da_r[k] * ulist_i[jjup] -
                        da_i[k] * ulist_r[jjup] +
                        a_r * dulist_i[jjup][k] -
                        a_i * dulist_r[jjup][k]);
          }

          rootpq = rootpqarray[ma + 1][j - mb];
          for (int k = 0; k < 3; k++) {
            dulist_r[jju+1][k] =
              -rootpq * (db_r[k] * ulist_r[jjup] +
                         db_i[k] * ulist_i[jjup] +
                         b_r * dulist_r[jjup][k] +
                         b_i * dulist_i[jjup][k]);
            dulist_i[jju+1][k] =
              -rootpq * (db_r[k] * ulist_i[jjup] -
                         db_i[k] * ulist_r[jjup] +
                         b_r * dulist_i[jjup][k] -
                         b_i * dulist_r[jjup][k]);
          }
          jju++;
          jjup++;
        }
        jju++;
      }

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

      jju = idxu_block[j];
      jjup = jju+(j+1)*(j+1)-1;
      int mbpar = 1;
      for (int mb = 0; 2*mb <= j; mb++) {
        int mapar = mbpar;
        for (int ma = 0; ma <= j; ma++) {
          if (mapar == 1) {
            for (int k = 0; k < 3; k++) {
              dulist_r[jjup][k] = dulist_r[jju][k];
              dulist_i[jjup][k] = -dulist_i[jju][k];
            }
          } else {
            for (int k = 0; k < 3; k++) {
              dulist_r[jjup][k] = -dulist_r[jju][k];
              dulist_i[jjup][k] = dulist_i[jju][k];
            }
          }
          mapar = -mapar;
          jju++;
          jjup--;
        }
        mbpar = -mbpar;
      }
    }

    double sfac  = snap_compute_sfac ( rmin0, switch_flag, switch_inner_flag, r, rcut, sinnerij[jj], dinnerij[jj] );
    double dsfac = snap_compute_dsfac( rmin0, switch_flag, switch_inner_flag, r, rcut, sinnerij[jj], dinnerij[jj] );

    sfac *= wj;
    dsfac *= wj;
    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j];
      for (int mb = 0; 2*mb <= j; mb++)
        for (int ma = 0; ma <= j; ma++) {
          dulist_r[jju][0] = dsfac * ulist_r[jju] * ux +
                                    sfac * dulist_r[jju][0];
          dulist_i[jju][0] = dsfac * ulist_i[jju] * ux +
                                    sfac * dulist_i[jju][0];
          dulist_r[jju][1] = dsfac * ulist_r[jju] * uy +
                                    sfac * dulist_r[jju][1];
          dulist_i[jju][1] = dsfac * ulist_i[jju] * uy +
                                    sfac * dulist_i[jju][1];
          dulist_r[jju][2] = dsfac * ulist_r[jju] * uz +
                                    sfac * dulist_r[jju][2];
          dulist_i[jju][2] = dsfac * ulist_i[jju] * uz +
                                    sfac * dulist_i[jju][2];
          jju++;
        }
    }
  }


  /* ----------------------------------------------------------------------
     calculate derivative of Ui w.r.t. atom j
  ------------------------------------------------------------------------- */

  static inline void snap_compute_duidrj( // READ ONLY
                                         int twojmax, int const * idxu_block
//                                       , int const * element
                                       , double const * drx, double const * dry, double const * drz
                                       , double const * rcutij, double const * wj
                                       , double const * const * ulist_r_ij, double const * const * ulist_i_ij, double const * const * rootpqarray
                                       , double const * sinnerij, double const * dinnerij
                                       , double rmin0, double rfac0, bool switch_flag, bool switch_inner_flag, bool chem_flag                             
                                         // WRITE ONLY
                                       , double * const * dulist_r, double * const * dulist_i
                                         // ORIGINAL PARAMETERS
                                       , int jj)
  {
    double rsq, r, x, y, z, z0, theta0, cs, sn;
    double dz0dr;
    double rcut = rcutij[jj];

    x = drx[jj];
    y = dry[jj];
    z = drz[jj];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);
    double rscale0 = rfac0 * M_PI / (rcut - rmin0);
    theta0 = (r - rmin0) * rscale0;
    cs = cos(theta0);
    sn = sin(theta0);
    z0 = r * cs / sn;
    dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

//    if (chem_flag) elem_duarray = element[jj];
//    else elem_duarray = 0;

    snap_compute_duarray( twojmax, idxu_block, ulist_r_ij, ulist_i_ij, rootpqarray
                        , sinnerij, dinnerij, rmin0, switch_flag, switch_inner_flag
                        , dulist_r, dulist_i
                        , x, y, z, z0, r, dz0dr, wj[jj], rcut, jj);
  }

}


