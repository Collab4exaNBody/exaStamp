#pragma once

#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  static inline void snap_compute_yi( // READ ONLY
                                      int nelements, int twojmax, int idxu_max, int const * idxu_block
                                    , int idxz_max, LAMMPS_NS::SNA_ZINDICES const * idxz, int const * const * const * idxcg_block, const double * cglist
                                    , double const * ulisttot_r, double const * ulisttot_i
                                    , int idxb_max, int const * const * const * idxb_block, bool bnorm_flag
                                      // WRITE ONLY
                                    , double * ylist_r, double * ylist_i
                                      // ORIGINAL PARAMETERS
                                    , const double* beta)
  {
    int jju;
    double betaj;
    int itriple;

    for (int ielem1 = 0; ielem1 < nelements; ielem1++)
      for (int j = 0; j <= twojmax; j++) {
        jju = idxu_block[j];
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            ylist_r[ielem1*idxu_max+jju] = 0.0;
            ylist_i[ielem1*idxu_max+jju] = 0.0;
            jju++;
          } // end loop over ma, mb
      } // end loop over j

    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++) {
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

            double ztmp_r = 0.0;
            double ztmp_i = 0.0;

            int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            int icgb = mb1min * (j2 + 1) + mb2max;
            for (int ib = 0; ib < nb; ib++) {

              double suma1_r = 0.0;
              double suma1_i = 0.0;

              const double *u1_r = &ulisttot_r[elem1*idxu_max+jju1];
              const double *u1_i = &ulisttot_i[elem1*idxu_max+jju1];
              const double *u2_r = &ulisttot_r[elem2*idxu_max+jju2];
              const double *u2_i = &ulisttot_i[elem2*idxu_max+jju2];

              int ma1 = ma1min;
              int ma2 = ma2max;
              int icga = ma1min * (j2 + 1) + ma2max;

              for (int ia = 0; ia < na; ia++) {
                suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
                suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
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

          jju = idxz[jjz].jju;
          for (int elem3 = 0; elem3 < nelements; elem3++) {
          // pick out right beta value
            if (j >= j1) {
              const int jjb = idxb_block[j1][j2][j];
              itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb;
              if (j1 == j) {
                if (j2 == j) betaj = 3*beta[itriple];
                else betaj = 2*beta[itriple];
              } else betaj = beta[itriple];
            } else if (j >= j2) {
              const int jjb = idxb_block[j][j2][j1];
              itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
              if (j2 == j) betaj = 2*beta[itriple];
              else betaj = beta[itriple];
            } else {
              const int jjb = idxb_block[j2][j][j1];
              itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
              betaj = beta[itriple];
            }

            if (!bnorm_flag && j1 > j)
              betaj *= (j1 + 1) / (j + 1.0);

            ylist_r[elem3 * idxu_max + jju] += betaj * ztmp_r;
            ylist_i[elem3 * idxu_max + jju] += betaj * ztmp_i;
          }
        } // end loop over jjz
      }

  }

}


