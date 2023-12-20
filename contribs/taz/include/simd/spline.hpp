/// @file
/// @brief Definition of a vectorized smoothstep function

#ifndef __SPLINE_HPP_INCLUDED
#define __SPLINE_HPP_INCLUDED


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


  	/// @brief Class gathering vectorized smoothstep functions
  	///
		/// The smoothstep function n is a 2n+1 order polynomial function
  	/// that takes value 0 at 0 and 1 at 1 and whose n first derivatives
  	/// take value 0 at 0 and 1
  	/// @tparam T Type of the variables
  	template <class T> class Spline {
  	public :

  		/// @brief Smoothstep 0 : identity (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> S0(vector_t<T> x) {
  			return x;
  		}

  		/// @brief Smoothstep 1 (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> S1(vector_t<T> x) {
  			return x*x * ( v3 - v2*x );
  		}

  		/// @brief Smoothstep 2 : smootherstep (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> S2(vector_t<T> x) {
  			vector_t<T> x2 = x*x;
  			return x2*x * ( v6*x2 - v15*x + v10 );
  		}

  		/// @brief Smoothstep 3 : smootheststep
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> S3(vector_t<T> x) {
  			vector_t<T> x2 = x*x;
  			return x2*x2 * ( v70*x2 - v84*x + v35 - v20*x2*x );
  		}

  		/// @brief Smoothstep 4 (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> S4(vector_t<T> x) {
  			vector_t<T> x2 = x *x;
  			vector_t<T> x3 = x2*x;
  			return x3*x2 * ( v70*x2*x2 - v315*x3 + v540*x2 - v420*x + v126 );
  		}

  		/// @brief Derivative of the smoothstep 0 (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> dS0(vector_t<T> x) {
  			return v1;
  		}

  		/// @brief Derivative of the smoothstep 1 (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> dS1(vector_t<T> x) {
  			return v6*x * ( v1 - x );
  		}

  		/// @brief Derivative of the smoothstep 2 (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> dS2(vector_t<T> x) {
  			vector_t<T> x2 = x*x;
  			return v30*x2 * ( x2 - v2*x + v1 );
  		}

  		/// @brief Derivative of the smoothstep 3
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> dS3(vector_t<T> x) {
  			vector_t<T> x2 = x *x;
  			vector_t<T> x3 = x2*x;
  			return v140*x3 * ( v1 - x3 + v3*x2 - v3*x );
  		}

  		/// @brief Derivative of the smoothstep 4 (not used)
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @return Output vector
  		static inline vector_t<T> dS4(vector_t<T> x) {
  			vector_t<T> x2 = x *x;
  			vector_t<T> x3 = x2*x;
  			vector_t<T> x4 = x2*x2;
  			return v630*x4 * ( x4 - v4*x3 + v6*x2 - v4*x + v1 );
  		}

  	private:

	  static const vector_t<T> v1; ///< Vector of 1
	  static const vector_t<T> v2; ///< Vector of 2
	  static const vector_t<T> v3; ///< Vector of 3
	  static const vector_t<T> v4; ///< Vector of 4
	  static const vector_t<T> v6; ///< Vector of 6
	  static const vector_t<T> v10; ///< Vector of 10
	  static const vector_t<T> v15; ///< Vector of 15
	  static const vector_t<T> v20; ///< Vector of 20
	  static const vector_t<T> v30; ///< Vector of 30
	  static const vector_t<T> v35; ///< Vector of 35
	  static const vector_t<T> v70; ///< Vector of 70
	  static const vector_t<T> v84; ///< Vector of 84
	  static const vector_t<T> v126; ///< Vector of 126
	  static const vector_t<T> v140; ///< Vector of 140
	  static const vector_t<T> v315; ///< Vector of 315
	  static const vector_t<T> v420; ///< Vector of 420
	  static const vector_t<T> v540; ///< Vector of 540
	  static const vector_t<T> v630; ///< Vector of 630

  	};


  	template <class T> const vector_t<T> Spline<T>::v1  (  1);
  	template <class T> const vector_t<T> Spline<T>::v2  (  2);
  	template <class T> const vector_t<T> Spline<T>::v3  (  3);
  	template <class T> const vector_t<T> Spline<T>::v4  (  4);
  	template <class T> const vector_t<T> Spline<T>::v6  (  6);
  	template <class T> const vector_t<T> Spline<T>::v10 ( 10);
  	template <class T> const vector_t<T> Spline<T>::v15 ( 15);
  	template <class T> const vector_t<T> Spline<T>::v20 ( 20);
  	template <class T> const vector_t<T> Spline<T>::v30 ( 30);
  	template <class T> const vector_t<T> Spline<T>::v35 ( 35);
  	template <class T> const vector_t<T> Spline<T>::v70 ( 70);
  	template <class T> const vector_t<T> Spline<T>::v84 ( 84);
  	template <class T> const vector_t<T> Spline<T>::v126(126);
  	template <class T> const vector_t<T> Spline<T>::v140(140);
  	template <class T> const vector_t<T> Spline<T>::v315(315);
  	template <class T> const vector_t<T> Spline<T>::v420(420);
  	template <class T> const vector_t<T> Spline<T>::v540(540);
  	template <class T> const vector_t<T> Spline<T>::v630(630);


  }


}

#endif // __SPLINE_HPP_INCLUDED
