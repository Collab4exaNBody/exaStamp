#pragma once

#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>
#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  /* ----------------------------------------------------------------------
     compute dEidRj
  ------------------------------------------------------------------------- */
  ONIKA_HOST_DEVICE_FUNC
  static inline void snap_compute_deidrj( // READ ONLY
                                          int elem_duarray, int twojmax, int idxu_max
                                        , int const * __restrict__ idxu_block
                                        , double const * __restrict__ dulist_r
                                        , double const * __restrict__ dulist_i
                                        , double const * __restrict__ ylist_r
                                        , double const * __restrict__ ylist_i
                                          // ORIGINAL PARAMETERS
                                        , double * __restrict__ dedr)
  {

    for (int k = 0; k < 3; k++)
      dedr[k] = 0.0;

    int jelem = elem_duarray;
    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j];

      for (int mb = 0; 2*mb < j; mb++)
        for (int ma = 0; ma <= j; ma++) {

          //double const * dudr_r = dulist_r[jju];
          //double const * dudr_i = dulist_i[jju];
          double jjjmambyarray_r = YLIST_R(jelem*idxu_max+jju);
          double jjjmambyarray_i = YLIST_I(jelem*idxu_max+jju);

          for (int k = 0; k < 3; k++)
            dedr[k] +=
              DULIST_R(jju,k) * jjjmambyarray_r +
              DULIST_I(jju,k) * jjjmambyarray_i;
          jju++;
        } //end loop over ma mb

      // For j even, handle middle column

      if (j%2 == 0) {

        int mb = j/2;
        for (int ma = 0; ma < mb; ma++) {
          //double const * dudr_r = dulist_r[jju];
          //double const * dudr_i = dulist_i[jju];
          double jjjmambyarray_r = YLIST_R(jelem*idxu_max+jju);
          double jjjmambyarray_i = YLIST_I(jelem*idxu_max+jju);

          for (int k = 0; k < 3; k++)
            dedr[k] +=
              DULIST_R(jju,k) * jjjmambyarray_r +
              DULIST_I(jju,k) * jjjmambyarray_i;
          jju++;
        }

        //double const * dudr_r = dulist_r[jju];
        //double const * dudr_i = dulist_i[jju];
        double jjjmambyarray_r = YLIST_R(jelem*idxu_max+jju);
        double jjjmambyarray_i = YLIST_I(jelem*idxu_max+jju);

        for (int k = 0; k < 3; k++)
          dedr[k] +=
            (DULIST_R(jju,k) * jjjmambyarray_r +
             DULIST_I(jju,k) * jjjmambyarray_i)*0.5;
        // jju++;

      } // end if jeven

    } // end loop over j

    for (int k = 0; k < 3; k++)
      dedr[k] *= 2.0;

  }

}


