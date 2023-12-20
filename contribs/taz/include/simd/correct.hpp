/// @file
/// @brief SIMD implementation of the distance corrector

#ifndef __CORRECT_HPP_INCLUDED
#define __CORRECT_HPP_INCLUDED


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


  /// @brief Correct distances with vectorisation
  /// @tparam T Type of the distance to correct
  /// @param [in,out] a_ Distances to correct
  /// @param [in] D_ Extension of the system
  /// @param [in] n Number of distances to correct
  template <class T>
  inline void correct (T* a_, const T& D_, const uint& n) {
      
      vector_t<T> a, out, tmp1, tmp2;
      
      vector_t<T> D (     D_);
      vector_t<T> pd( 0.5*D_);
      vector_t<T> md(-0.5*D_);
      
      for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {
	
      	a.load(a_+i);

#ifdef __evi_mic
      	tmp1 = blend(a > pd, vector_t<T>::zero(), D);
      	tmp2 = blend(a < md, vector_t<T>::zero(), D);
#else
      	tmp1 = blendv(vector_t<T>::zero(), D, a > pd);
      	tmp2 = blendv(vector_t<T>::zero(), D, a < md);
#endif
	
      	out = tmp2 - tmp1 + a;
	
      	out.store(a_+i);
	
      }
      
    }

  /// @brief Correct positions with vectorisation
  /// @tparam T Type of the distance to correct
  /// @param [in] ref_ Positions to compare to
  /// @param [in,out] p_ Positions to correct
  /// @param [in] D_ Extension of the system
  /// @param [in] n Number of distances to correct
  template <class T>
  inline void correctP (T* ref_, T* p_, const T& D_, const uint& n) {

      vector_t<T> ref, p, out,tmpd, tmp1, tmp2;

      vector_t<T> D (     D_);
      vector_t<T> pd( 0.5*D_);
      vector_t<T> md(-0.5*D_);

      for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

      	ref.load(ref_+i);
      	p.load(p_+i);
      	tmpd=p-ref;

#ifdef __evi_mic
      	tmp1 = blend(tmpd > pd, vector_t<T>::zero(), D);
      	tmp2 = blend(tmpd < md, vector_t<T>::zero(), D);
#else
      	tmp1 = blendv(vector_t<T>::zero(), D, tmpd > pd);
      	tmp2 = blendv(vector_t<T>::zero(), D, tmpd < md);
#endif

      	out = tmp2 - tmp1 + p;

      	out.store(p_+i);

      }

    }

  }


}
  
#endif // __CORRECT_HPP_INCLUDED
