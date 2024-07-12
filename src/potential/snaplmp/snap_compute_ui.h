#pragma once

#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  static inline void snap_zero_uarraytot( // READ ONLY
                                          int nelements
                                        , int twojmax
                                        , int const * idxu_block
                                        , int idxu_max
                                        , double wself
                                        , bool wselfall_flag
                                          // WRITE ONLY
                                        , double * ulisttot_r, double * ulisttot_i
                                          // ORIGINAL PARAMETERS
                                        , int ielem ) 
  {
    for (int jelem = 0; jelem < nelements; jelem++)
    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j];
      for (int mb = 0; mb <= j; mb++) {
        for (int ma = 0; ma <= j; ma++) {
          ulisttot_r[jelem*idxu_max+jju] = 0.0;
          ulisttot_i[jelem*idxu_max+jju] = 0.0;

          // utot(j,ma,ma) = wself, sometimes
          if (jelem == ielem || wselfall_flag)
            if (ma==mb)
            ulisttot_r[jelem*idxu_max+jju] = wself; ///// double check this
          jju++;
        }
      }
    }
  }

  static inline void snap_compute_uarray( // READ ONLY
                                          int twojmax
                                        , int const * idxu_block
                                        , double const * const * rootpqarray
                                          // WRITE ONLY
                                        , double * const * ulist_r_ij, double * const * ulist_i_ij
                                          // ORIGINAL PARAMETERS
                                        , double x, double y, double z, double z0, double r, int jj )
  {
    double r0inv;
    double a_r, b_r, a_i, b_i;
    double rootpq;

    // compute Cayley-Klein parameters for unit quaternion

    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    // VMK Section 4.8.2
    //double* ulist_r = ulist_r_ij[jj];
    //double* ulist_i = ulist_i_ij[jj];

    ULIST_R_IJ(jj,0) = 1.0;
    ULIST_I_IJ(jj,0) = 0.0;

    for (int j = 1; j <= twojmax; j++) {
      int jju = idxu_block[j];
      int jjup = idxu_block[j-1];

      // fill in left side of matrix layer from previous layer

      for (int mb = 0; 2*mb <= j; mb++) {
        ULIST_R_IJ(jj,jju) = 0.0;
        ULIST_I_IJ(jj,jju) = 0.0;

        for (int ma = 0; ma < j; ma++) {
          rootpq = rootpqarray[j - ma][j - mb];
          ULIST_R_IJ(jj,jju) +=
            rootpq *
            (a_r * ULIST_R_IJ(jj,jjup) +
             a_i * ULIST_I_IJ(jj,jjup));
          ULIST_I_IJ(jj,jju) +=
            rootpq *
            (a_r * ULIST_I_IJ(jj,jjup) -
             a_i * ULIST_R_IJ(jj,jjup));

          rootpq = rootpqarray[ma + 1][j - mb];
          ULIST_R_IJ(jj,jju+1) =
            -rootpq *
            (b_r * ULIST_R_IJ(jj,jjup) +
             b_i * ULIST_I_IJ(jj,jjup));
          ULIST_I_IJ(jj,jju+1) =
            -rootpq *
            (b_r * ULIST_I_IJ(jj,jjup) -
             b_i * ULIST_R_IJ(jj,jjup));
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
            ULIST_R_IJ(jj,jjup) = ULIST_R_IJ(jj,jju);
            ULIST_I_IJ(jj,jjup) = -ULIST_I_IJ(jj,jju);
          } else {
            ULIST_R_IJ(jj,jjup) = -ULIST_R_IJ(jj,jju);
            ULIST_I_IJ(jj,jjup) = ULIST_I_IJ(jj,jju);
          }
          mapar = -mapar;
          jju++;
          jjup--;
        }
        mbpar = -mbpar;
      }
    }
  }

  static inline double snap_compute_sfac( // READ ONLY
                                          double rmin0, bool switch_flag, bool switch_inner_flag
                                          // ORIGINAL PARAMTERS
                                        , double r, double rcut, double sinner, double dinner)
  {
    double sfac = 0.0;

    // calculate sfac = sfac_outer

    if (switch_flag == 0) sfac = 1.0;
    else if (r <= rmin0) sfac = 1.0;
    else if (r > rcut) sfac = 0.0;
    else {
      double rcutfac = M_PI / (rcut - rmin0);
      sfac = 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }

    // calculate sfac *= sfac_inner, rarely visited

    if (switch_inner_flag == 1 && r < sinner + dinner) {
      if (r > sinner - dinner) {
        double rcutfac = (M_PI/2) / dinner;
        sfac *= 0.5 * (1.0 - cos( (M_PI/2) + (r - sinner) * rcutfac));
      } else sfac = 0.0;
    }

    return sfac;
  }

  static inline void snap_add_uarraytot( // READ ONLY
                                         double rmin0, bool switch_flag, bool switch_inner_flag
                                       , int twojmax, bool chem_flag, int const * element
                                       , double const * rcutij, double const * sinnerij, double const * dinnerij, double const * wj
                                       , int const * idxu_block , int idxu_max
                                       , double const * const * ulist_r_ij, double const * const * ulist_i_ij
                                         // WRITE ONLY
                                       , double * ulisttot_r, double * ulisttot_i
                                         // ORIGINAL PARAMETERS
                                       , double r, int jj)
  {
    double sfac;
    int jelem;

    sfac = snap_compute_sfac( rmin0, switch_flag, switch_inner_flag, r, RCUTIJ(jj), sinnerij[jj], dinnerij[jj]);
    sfac *= WJ(jj);

    if (chem_flag) jelem = element[jj];
    else jelem = 0;

    //double const * ulist_r = ulist_r_ij[jj];
    //double const * ulist_i = ulist_i_ij[jj];

    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j];
      for (int mb = 0; mb <= j; mb++)
        for (int ma = 0; ma <= j; ma++) {
          ulisttot_r[jelem*idxu_max+jju] +=
            sfac * ULIST_R_IJ(jj,jju);
          ulisttot_i[jelem*idxu_max+jju] +=
            sfac * ULIST_I_IJ(jj,jju);
          jju++;
        }
    }
  }


  static inline void snap_compute_ui( // READ ONLY
                                      int nelements, int twojmax, int idxu_max, int const * idxu_block, int const * element
                                    , double const * drx, double const * dry, double const * drz
                                    , double const * rcutij, double const * const * rootpqarray
                                    , double const * sinnerij, double const * dinnerij, double const * wj
                                    , bool wselfall_flag, bool switch_flag, bool switch_inner_flag, bool chem_flag
                                    , double wself, double rmin0, double rfac0
                                      // SCRATCH BUFFER
                                    , double * const * ulist_r_ij, double * const * ulist_i_ij
                                      // WRITE ONLY
                                    , double * ulisttot_r, double * ulisttot_i
                                      // ORIGINAL PARAMETERS
                                    , int jnum, int ielem)
  {
    // utot(j,ma,mb) = 0 for all j,ma,ma
    // utot(j,ma,ma) = 1 for all j,ma
    // for j in neighbors of i:
    //   compute r0 = (x,y,z,z0)
    //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

    snap_zero_uarraytot( nelements, twojmax, idxu_block, idxu_max, wself, wselfall_flag, ulisttot_r, ulisttot_i, ielem );
    
    for (int j = 0; j < jnum; j++)
    {
      const double x = drx[j];
      const double y = dry[j];
      const double z = drz[j];
      const double rsq = x * x + y * y + z * z;
      const double r = sqrt(rsq);

      const double theta0 = (r - rmin0) * rfac0 * M_PI / (RCUTIJ(j) - rmin0);
      const double z0 = r / tan(theta0);

      snap_compute_uarray( twojmax, idxu_block, rootpqarray, ulist_r_ij, ulist_i_ij, x, y, z, z0, r, j);
      
      snap_add_uarraytot( rmin0, switch_flag, switch_inner_flag, twojmax, chem_flag, element
                        , rcutij, sinnerij, dinnerij, wj
                        , idxu_block, idxu_max, ulist_r_ij, ulist_i_ij, ulisttot_r, ulisttot_i
                        , r, j);
    }
  }

}


