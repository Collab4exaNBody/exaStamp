#pragma once

#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>
#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  ONIKA_HOST_DEVICE_FUNC
  static inline void snap_zero_yi_array( int nelements, int idxu_max, double * __restrict__ ylist_r, double * __restrict__ ylist_i )
  {
    const int N = idxu_max * nelements;
    for(int i=0;i<N;++i)
    {
      YLIST_R(i) = 0.0;
      YLIST_I(i) = 0.0;
    }
    /*
    for (int ielem1 = 0; ielem1 < nelements; ielem1++)
    {
      for (int j = 0; j <= twojmax; j++)
      {
        int jju = IDXU_BLOCK(j);
        for (int mb = 0; 2*mb <= j; mb++)
        {
          for (int ma = 0; ma <= j; ma++)
          {
            YLIST_R(ielem1*idxu_max+jju) = 0.0;
            YLIST_I(ielem1*idxu_max+jju) = 0.0;
            jju++;
          } // end loop over ma, mb
        }
      } // end loop over j
    }
    */
  }

  ONIKA_HOST_DEVICE_FUNC
  static inline void snap_add_yi_contribution_alt( // READ ONLY
                                      int nelements, int twojmax, int idxu_max, int idxz_max
                                    , SnapInternal::SNA_ZINDICES_ALT const * __restrict__ const idxz_alt
                                    , int const * __restrict__ const idxcg_block
                                    , const double * __restrict__ const cglist
                                    , const int * __restrict__ const y_jju_map, int idxu_max_alt
                                    , double const * __restrict__ const ulisttot_r
                                    , double const * __restrict__ const ulisttot_i
                                    , int idxb_max
                                    , int const * __restrict__ idxb_block
                                    , bool bnorm_flag
                                      // WRITE ONLY
                                    , double * __restrict__ const ylist_r
                                    , double * __restrict__ const ylist_i
                                      // ORIGINAL PARAMETERS
                                    , const double * __restrict__ const beta)
  {
//    int jju;
//    double betaj;
//    int itriple;

//    snap_zero_yi_array( nelements, idxu_max, ylist_r, ylist_i );

    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++)
      {
        for (int jjz = 0; jjz < idxz_max; jjz++)
        {
          const int j1 = IDXZ_ALT(jjz).j1;
          const int j2 = IDXZ_ALT(jjz).j2;
          const int j = IDXZ_ALT(jjz).j;
          const int ma1min = IDXZ_ALT(jjz).ma1min;
          const int ma2max = IDXZ_ALT(jjz).ma2max;
          const int na = IDXZ_ALT(jjz).na;
          const int mb1min = IDXZ_ALT(jjz).mb1min;
          const int mb2max = IDXZ_ALT(jjz).mb2max;
          const int nb = IDXZ_ALT(jjz).nb;
          const int jju_map = IDXZ_ALT(jjz).jju;

          assert( jju_map != -1 );

          const double * __restrict__ const cgblock = cglist + IDXZ_ALT(jjz).idx_cgblock; //idxcg_block_j1_j2_j;

          double ztmp_r = 0.0;
          double ztmp_i = 0.0;

          int jju1 = IDXU_BLOCK(j1) + (j1 + 1) * mb1min;
          int jju2 = IDXU_BLOCK(j2) + (j2 + 1) * mb2max;
          int icgb = mb1min * (j2 + 1) + mb2max;
                    
          for (int ib = 0; ib < nb; ib++)
          {
            double suma1_r = 0.0;
            double suma1_i = 0.0;

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;
            
            for (int ia = 0; ia < na; ia++)
            {
              // suma1 += cgblock[icga] * UiTot(elem1,jju1+ma1) * UiTot(elem2,jju2+ma2)
              suma1_r += cgblock[icga] * (ULISTTOT_R(elem1*idxu_max+jju1+ma1) * ULISTTOT_R(elem2*idxu_max+jju2+ma2) \
                       - ULISTTOT_I(elem1*idxu_max+jju1+ma1) * ULISTTOT_I(elem2*idxu_max+jju2+ma2));
              suma1_i += cgblock[icga] * (ULISTTOT_R(elem1*idxu_max+jju1+ma1) * ULISTTOT_I(elem2*idxu_max+jju2+ma2) \
                       + ULISTTOT_I(elem1*idxu_max+jju1+ma1) * ULISTTOT_R(elem2*idxu_max+jju2+ma2));
              ma1++;
              ma2--;
              icga += j2;
            } // end loop over ia

            ztmp_r += cgblock[icgb] * suma1_r;
            ztmp_i += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
          } // end loop over ib
          
          // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
          // find right y_list[jju] and beta[jjb] entries
          // multiply and divide by j+1 factors
          // account for multiplicity of 1, 2, or 3

          if (bnorm_flag) {
            ztmp_i /= j+1;
            ztmp_r /= j+1;
          }

          const int jjb = IDXZ_ALT(jjz).jjb;
          
          for (int elem3 = 0; elem3 < nelements; elem3++) {
            double betaj = 0.0;
            // pick out right beta value
            if (j >= j1) {
              //const int jjb = IDXB_BLOCK(j1,j2,j);
              const int itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb;
              if (j1 == j) {
                if (j2 == j) betaj = 3*beta[itriple];
                else betaj = 2*beta[itriple];
              } else betaj = beta[itriple];
            } else if (j >= j2) {
              //const int jjb = IDXB_BLOCK(j,j2,j1);
              const int itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
              if (j2 == j) betaj = 2*beta[itriple];
              else betaj = beta[itriple];
            } else {
              //const int jjb = IDXB_BLOCK(j2,j,j1);
              const int itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
              betaj = beta[itriple];
            }

            if (!bnorm_flag && j1 > j)
              betaj *= (j1 + 1) / (j + 1.0);
            
            assert( jju_map < idxu_max_alt );
            YLIST_R(elem3 * idxu_max_alt + jju_map) += betaj * ztmp_r;
            YLIST_I(elem3 * idxu_max_alt + jju_map) += betaj * ztmp_i;
          }      
          
        } // end loop over jjz
      }

  }

  ONIKA_HOST_DEVICE_FUNC
  static inline void snap_add_yi_contribution( // READ ONLY
                                      int nelements, int twojmax, int idxu_max, int idxz_max
                                    , SnapInternal::SNA_ZINDICES const * __restrict__ const idxz
                                    , int const * __restrict__ const idxcg_block
                                    , const double * __restrict__ const cglist
                                    , double const * __restrict__ const ulisttot_r
                                    , double const * __restrict__ const ulisttot_i
                                    , int idxb_max
                                    , int const * __restrict__ idxb_block
                                    , bool bnorm_flag
                                      // WRITE ONLY
                                    , double * __restrict__ const ylist_r
                                    , double * __restrict__ const ylist_i
                                      // ORIGINAL PARAMETERS
                                    , const double * __restrict__ const beta)
  {
//    int jju;
//    double betaj;
//    int itriple;

//    snap_zero_yi_array( nelements, idxu_max, ylist_r, ylist_i );

    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++) {
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

            const double * __restrict__ const cgblock = cglist + IDXCG_BLOCK(j1,j2,j);

            double ztmp_r = 0.0;
            double ztmp_i = 0.0;

            int jju1 = IDXU_BLOCK(j1) + (j1 + 1) * mb1min;
            int jju2 = IDXU_BLOCK(j2) + (j2 + 1) * mb2max;
            int icgb = mb1min * (j2 + 1) + mb2max;
            for (int ib = 0; ib < nb; ib++) {

              double suma1_r = 0.0;
              double suma1_i = 0.0;

              int ma1 = ma1min;
              int ma2 = ma2max;
              int icga = ma1min * (j2 + 1) + ma2max;

              for (int ia = 0; ia < na; ia++) {
                // suma1 += cgblock[icga] * UiTot(elem1,jju1+ma1) * UiTot(elem2,jju2+ma2)
                suma1_r += cgblock[icga] * (ULISTTOT_R(elem1*idxu_max+jju1+ma1) * ULISTTOT_R(elem2*idxu_max+jju2+ma2) - ULISTTOT_I(elem1*idxu_max+jju1+ma1) * ULISTTOT_I(elem2*idxu_max+jju2+ma2));
                suma1_i += cgblock[icga] * (ULISTTOT_R(elem1*idxu_max+jju1+ma1) * ULISTTOT_I(elem2*idxu_max+jju2+ma2) + ULISTTOT_I(elem1*idxu_max+jju1+ma1) * ULISTTOT_R(elem2*idxu_max+jju2+ma2));
                ma1++;
                ma2--;
                icga += j2;
              } // end loop over ia

              // ztmp += cgblock[icgb] * suma1
              ztmp_r += cgblock[icgb] * suma1_r;
              ztmp_i += cgblock[icgb] * suma1_i;

              jju1 += j1 + 1;
              jju2 -= j2 + 1;
              icgb += j2;
            } // end loop over ib

            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
            // find right y_list[jju] and beta[jjb] entries
            // multiply and divide by j+1 factors
            // account for multiplicity of 1, 2, or 3

          if (bnorm_flag) {
            ztmp_i /= j+1;
            ztmp_r /= j+1;
          }

          const int jju = IDXZ(jjz).jju;
          for (int elem3 = 0; elem3 < nelements; elem3++) {
            double betaj = 0.0;
            // pick out right beta value
            if (j >= j1) {
              const int jjb = IDXB_BLOCK(j1,j2,j);
              const int itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb;
              if (j1 == j) {
                if (j2 == j) betaj = 3*beta[itriple];
                else betaj = 2*beta[itriple];
              } else betaj = beta[itriple];
            } else if (j >= j2) {
              const int jjb = IDXB_BLOCK(j,j2,j1);
              const int itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
              if (j2 == j) betaj = 2*beta[itriple];
              else betaj = beta[itriple];
            } else {
              const int jjb = IDXB_BLOCK(j2,j,j1);
              const int itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
              betaj = beta[itriple];
            }

            if (!bnorm_flag && j1 > j)
              betaj *= (j1 + 1) / (j + 1.0);

            YLIST_R(elem3 * idxu_max + jju) += betaj * ztmp_r;
            YLIST_I(elem3 * idxu_max + jju) += betaj * ztmp_i;
          }
        } // end loop over jjz
      }

  }

  ONIKA_HOST_DEVICE_FUNC
  static inline void snap_compute_yi( // READ ONLY
                                      int nelements, int twojmax, int idxu_max, int idxz_max
                                    , SnapInternal::SNA_ZINDICES const * __restrict__ idxz
                                    , int const * __restrict__ const idxcg_block
                                    , const double * __restrict__ cglist
                                    , double const * __restrict__ ulisttot_r
                                    , double const * __restrict__ ulisttot_i
                                    , int idxb_max
                                    , int const * __restrict__ idxb_block
                                    , bool bnorm_flag
                                      // WRITE ONLY
                                    , double * __restrict__ ylist_r
                                    , double * __restrict__ ylist_i
                                      // ORIGINAL PARAMETERS
                                    , const double * __restrict__ beta)
  {
    snap_zero_yi_array( nelements, idxu_max, ylist_r, ylist_i );
    snap_add_yi_contribution( nelements, twojmax, idxu_max, idxz_max, idxz, idxcg_block, cglist, ulisttot_r, ulisttot_i, idxb_max, idxb_block, bnorm_flag, ylist_r, ylist_i, beta );
  }


}


