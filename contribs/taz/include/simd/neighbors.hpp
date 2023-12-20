/// @file
/// @brief Vectorization of neighbors computation

#ifndef __NEIGHBORS_HPP_INCLUDED
#define __NEIGHBORS_HPP_INCLUDED


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


  	/// @brief Vectorized computation of the neighbors of a particle
  	/// @tparam T Type of the variables
  	/// @param [in] ri_ Position of the particle
  	/// @param [in] rx_ Y components of the positions of the potential neighbors
  	/// @param [in] ry_ Y components of the positions of the potential neighbors
  	/// @param [in] rz_ Z components of the positions of the potential neighbors
  	/// @param [in] rcut_ Cutoff radius for each neighbor
  	/// @param [out] out_ Boolean results (is or isn't a neighbor)
  	/// @param [in] n Number of potential neighbors
  	/// @warning The boolean results are stored into a array of T (double)
    template <class T>
    inline void neighbors (const double *ri_, const T *rx_, const T *ry_, const T *rz_, const T *rcut_, T *out_, const uint n) {
      vector_t<T> rxi(ri_[0]);
      vector_t<T> ryi(ri_[1]);
      vector_t<T> rzi(ri_[2]);

      vector_t<T> rxj, ryj, rzj, rcj, out;

#if defined __evi_mic
      typename vector_t<T>::mask_t msk;
#endif

      for (uint j=0; j<n; j+=vector_t<T>::chunk_size) {

      	rxj.load(rx_  + j);
      	ryj.load(ry_  + j);
      	rzj.load(rz_  + j);
      	rcj.load(rcut_+ j);

      	rxj = rxi - rxj;
      	ryj = ryi - ryj;
      	rzj = rzi - rzj;

#if defined __evi_mic
      	msk = (rxj*rxj + ryj*ryj + rzj*rzj) < rcj;
      	out = blend(msk, vector_t<T>::zero(), vector_t<T>::one());
#else
        out = (rxj*rxj + ryj*ryj + rzj*rzj) < rcj;
#endif

      	out.store(out_+j);

      }

    }


  }


}

#endif // __NEIGHBORS_HPP_INCLUDED
