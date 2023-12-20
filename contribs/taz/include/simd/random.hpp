/// @file
/// @brief Vectorized generation of random numbers

#ifndef __RANDOM_HPP_INCLUDED
#define __RANDOM_HPP_INCLUDED


#include "libevi/simd.hpp"
#include "saru/saru.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Tool to handle vectorized random number generation
	/// @tparam T Type of the variables
  template <class T>
  class random {

  private:

    static vector_t<T> rnd_max_pls_two; ///< Vector of 2*RAND_MAX
    static vector_t<T> inv_2_rnd_max_pls_three; ///< Vector of inverse of 2*RAND_MAX+3

    static vector_t<T> minus_two; ///< Vector of -2
    static vector_t<T> two_pi; ///< Vector of two pi

  public:

    ///@brief Convert random integers in [-RAND_MAX, RAND_MAX[ with an uniform distribution
    /// into random reals of type T in ]0, 1] with an uniform distribution
    /// @tparam T Type of the variables
    /// @param [in] rnd Random integers
    /// @return Random reals of type T
    static inline vector_t<T> make_uniform(const double* rnd) {
      vector_t<T> v;
      v.load(rnd);
      return (v + rnd_max_pls_two) * inv_2_rnd_max_pls_three;
    }

    ///@brief Convert two independent sets of random integers in [-RAND_MAX, RAND_MAX[
    /// with an uniform distribution into random reals of type T in ]0, 1] with an normal distribution
    /// @tparam T Type of the variables
    /// @param [in] rnd1 Random integers
    /// @param [in] rnd2 Random integers
    /// @return Random reals of type T
    static inline vector_t<T> make_normal(const double* rnd1, const double* rnd2) {
      vector_t<T> u, v;
      u.load(rnd1);
      v.load(rnd2);
      u = (u + rnd_max_pls_two) * inv_2_rnd_max_pls_three;
      v = (v + rnd_max_pls_two) * inv_2_rnd_max_pls_three;
      return sqrt( minus_two * log(u) ) * cos( two_pi * v );
    }

  };


  template <class T> vector_t<T> random<T>::rnd_max_pls_two         = vector_t<T>( Saru::RAND_MAX_PLS_TWO );
  template <class T> vector_t<T> random<T>::inv_2_rnd_max_pls_three = vector_t<T>( Saru::INV_2_RAND_MAX_PLS_THREE );
  template <class T> vector_t<T> random<T>::minus_two               = vector_t<T>( -2.00 );
  template <class T> vector_t<T> random<T>::two_pi                  = vector_t<T>( Saru::TWO_PI );
  

}

#endif // __RANDOM_HPP_INCLUDED
