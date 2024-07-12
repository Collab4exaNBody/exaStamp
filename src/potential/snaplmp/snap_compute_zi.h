#pragma once

#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  /* ----------------------------------------------------------------------
     compute Zi by summing over products of Ui
  ------------------------------------------------------------------------- */

  static inline void snap_compute_zi( // READ ONLY
                                      int nelements, int idxz_max, int idxu_max, int const * idxu_block, int const * const * const * idxcg_block
                                    , LAMMPS_NS::SNA_ZINDICES const * idxz
                                    , const double * cglist, double const * ulisttot_r, double const * ulisttot_i
                                    , bool bnorm_flag
                                    // WRITE ONLY
                                    , double * zlist_r, double * zlist_i )
  {

    int idouble = 0;
    double * zptr_r;
    double * zptr_i;
    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++) {

        zptr_r = &zlist_r[idouble*idxz_max];
        zptr_i = &zlist_i[idouble*idxz_max];

        for (int jjz = 0; jjz < idxz_max; jjz++) {
          const int j1 = idxz[jjz].j1;
          const int j2 = idxz[jjz].j2;
          const int j = idxz[jjz].j;
          const int ma1min = idxz[jjz].ma1min;
          const int ma2max = idxz[jjz].ma2max;
          const int na = idxz[jjz].na;
          const int mb1min = idxz[jjz].mb1min;
          const int mb2max = idxz[jjz].mb2max;
          const int nb = idxz[jjz].nb;

          const double *cgblock = cglist + idxcg_block[j1][j2][j];

          zptr_r[jjz] = 0.0;
          zptr_i[jjz] = 0.0;

          int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
          int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
          int icgb = mb1min * (j2 + 1) + mb2max;
          for (int ib = 0; ib < nb; ib++) {

            double suma1_r = 0.0;
            double suma1_i = 0.0;

            //const double *u1_r = &ULISTTOT_R(elem1*idxu_max+jju1);
            //const double *u1_i = &ULISTTOT_I(elem1*idxu_max+jju1);
            //const double *u2_r = &ULISTTOT_R(elem2*idxu_max+jju2);
            //const double *u2_i = &ULISTTOT_I(elem2*idxu_max+jju2);

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
              suma1_r += cgblock[icga] * (ULISTTOT_R(elem1*idxu_max+jju1+ma1) * ULISTTOT_R(elem2*idxu_max+jju2+ma2) - ULISTTOT_I(elem1*idxu_max+jju1+ma1) * ULISTTOT_I(elem2*idxu_max+jju2+ma2));
              suma1_i += cgblock[icga] * (ULISTTOT_R(elem1*idxu_max+jju1+ma1) * ULISTTOT_I(elem2*idxu_max+jju2+ma2) + ULISTTOT_I(elem1*idxu_max+jju1+ma1) * ULISTTOT_R(elem2*idxu_max+jju2+ma2));
              ma1++;
              ma2--;
              icga += j2;
            } // end loop over ia

            zptr_r[jjz] += cgblock[icgb] * suma1_r;
            zptr_i[jjz] += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
          } // end loop over ib
          if (bnorm_flag) {
            zptr_r[jjz] /= (j+1);
            zptr_i[jjz] /= (j+1);
          }
        } // end loop over jjz
        idouble++;
      }
  }


}


