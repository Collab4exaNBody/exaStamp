#pragma once

#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  /* ----------------------------------------------------------------------
     compute Bi by summing conj(Ui)*Zi
  ------------------------------------------------------------------------- */

  static inline void snap_compute_bi( // READ ONLY
                                      int nelements, int idxz_max, int idxb_max, int idxu_max
                                    , int const * idxu_block, int const * const * const * idxz_block
                                    , LAMMPS_NS::SNA_ZINDICES const * idxz, LAMMPS_NS::SNA_BINDICES const * idxb
                                    , double const * zlist_r, double const * zlist_i
                                    , double const * ulisttot_r, double const * ulisttot_i
                                    , double const * bzero , bool bzero_flag, bool wselfall_flag
                                      // WRITE ONLY
                                    , double * blist
                                      // ORIGINAL PARAMETERS
                                    , int ielem)
  {
    // for j1 = 0,...,twojmax
    //   for j2 = 0,twojmax
    //     for j = |j1-j2|,Min(twojmax,j1+j2),2
    //        b(j1,j2,j) = 0
    //        for mb = 0,...,jmid
    //          for ma = 0,...,j
    //            b(j1,j2,j) +=
    //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

    int itriple = 0;
    int idouble = 0;
    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++) {

        //double const *zptr_r = &ZLIST_R(idouble*idxz_max);
        //double const *zptr_i = &ZLIST_I(idouble*idxz_max);

        for (int elem3 = 0; elem3 < nelements; elem3++) {
          for (int jjb = 0; jjb < idxb_max; jjb++) {
            const int j1 = idxb[jjb].j1;
            const int j2 = idxb[jjb].j2;
            const int j = idxb[jjb].j;

            int jjz = idxz_block[j1][j2][j];
            int jju = idxu_block[j];
            double sumzu = 0.0;
            for (int mb = 0; 2 * mb < j; mb++)
              for (int ma = 0; ma <= j; ma++) {
                sumzu += ULISTTOT_R(elem3*idxu_max+jju) * ZLIST_R(idouble*idxz_max+jjz) +
                         ULISTTOT_I(elem3*idxu_max+jju) * ZLIST_I(idouble*idxz_max+jjz);
                jjz++;
                jju++;
              } // end loop over ma, mb

            // For j even, handle middle column

            if (j % 2 == 0) {
              int mb = j / 2;
              for (int ma = 0; ma < mb; ma++) {
                sumzu += ULISTTOT_R(elem3*idxu_max+jju) * ZLIST_R(idouble*idxz_max+jjz) +
                         ULISTTOT_I(elem3*idxu_max+jju) * ZLIST_I(idouble*idxz_max+jjz);
                jjz++;
                jju++;
              }

              sumzu += 0.5 * (ULISTTOT_R(elem3*idxu_max+jju) * ZLIST_R(idouble*idxz_max+jjz) +
                              ULISTTOT_I(elem3*idxu_max+jju) * ZLIST_I(idouble*idxz_max+jjz));
            } // end if jeven

            BLIST(itriple*idxb_max+jjb) = 2.0 * sumzu;

          }
          itriple++;
        }
        idouble++;
      }

    // apply bzero shift

    if (bzero_flag) {
      if (!wselfall_flag) {
        itriple = (ielem*nelements+ielem)*nelements+ielem;
        for (int jjb = 0; jjb < idxb_max; jjb++) {
          const int j = idxb[jjb].j;
          BLIST(itriple*idxb_max+jjb) -= bzero[j];
        } // end loop over JJ
      } else {
        int itriple = 0;
        for (int elem1 = 0; elem1 < nelements; elem1++)
          for (int elem2 = 0; elem2 < nelements; elem2++) {
            for (int elem3 = 0; elem3 < nelements; elem3++) {
              for (int jjb = 0; jjb < idxb_max; jjb++) {
                const int j = idxb[jjb].j;
                BLIST(itriple*idxb_max+jjb) -= bzero[j];
              } // end loop over JJ
              itriple++;
            } // end loop over elem3
          } // end loop over elem1,elem2
      }
    }
  }

}


