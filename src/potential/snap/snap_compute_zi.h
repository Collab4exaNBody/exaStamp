#pragma once

#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>
#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  /* ----------------------------------------------------------------------
     compute Zi by summing over products of Ui
  ------------------------------------------------------------------------- */
  ONIKA_HOST_DEVICE_FUNC
  static inline void snap_compute_zi( // READ ONLY
                                      int nelements, int idxz_max, int idxu_max, int twojmax
//                                    , int const * __restrict__ idxu_block
                                    , int const * __restrict__ const idxcg_block
                                    , SnapInternal::SNA_ZINDICES const * __restrict__ idxz
                                    , const double * __restrict__ cglist
                                    , double const * __restrict__ ulisttot_r
                                    , double const * __restrict__ ulisttot_i
                                    , bool bnorm_flag
                                    // WRITE ONLY
                                    , double * __restrict__ zlist_r
                                    , double * __restrict__ zlist_i )
  {

    int idouble = 0;
    //double * zptr_r;
    //double * zptr_i;
    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++) {

        //zptr_r = &ZLIST_R(idouble*idxz_max);
        //zptr_i = &ZLIST_I(idouble*idxz_max);

        for (int jjz = 0; jjz < idxz_max; jjz++) {
          const int j1 = IDXZ(jjz).j1;
          const int j2 = IDXZ(jjz).j2;
          const int j = IDXZ(jjz).j;
          const int ma1min = IDXZ(jjz).ma1min;
          const int ma2max = IDXZ(jjz).ma2max;
          const int na = IDXZ(jjz).na;
          const int mb1min = IDXZ(jjz).mb1min;
          const int mb2max = IDXZ(jjz).mb2max;
          const int nb = IDXZ(jjz).nb;

          const double * const __restrict__ cgblock = cglist + IDXCG_BLOCK(j1,j2,j);

          ZLIST_R(idouble*idxz_max+jjz) = 0.0;
          ZLIST_I(idouble*idxz_max+jjz) = 0.0;

          int jju1 = IDXU_BLOCK(j1) + (j1 + 1) * mb1min;
          int jju2 = IDXU_BLOCK(j2) + (j2 + 1) * mb2max;
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

            ZLIST_R(idouble*idxz_max+jjz) += cgblock[icgb] * suma1_r;
            ZLIST_I(idouble*idxz_max+jjz) += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
          } // end loop over ib
          if (bnorm_flag) {
            ZLIST_R(idouble*idxz_max+jjz) /= (j+1);
            ZLIST_I(idouble*idxz_max+jjz) /= (j+1);
          }
        } // end loop over jjz
        idouble++;
      }
  }


}


