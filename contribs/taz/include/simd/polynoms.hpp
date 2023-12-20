/*
 * polynoms.hpp
 *
 *  Created on: Mar 13, 2017
 *      Author: giarda
 */

/// @file
/// @brief Definition of a vectorized smoothstep function

#ifndef POLYNOMS_HPP_
#define POLYNOMS_HPP_


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


  	/// @brief Class gathering tools to compute various polynomial functions
  	///
  	/// @tparam T Type of the variables
  	template <class T> class Polynoms {
  	public :

  		/// @brief Emu : 1 + x + coeff * x^3
  		/// @param [in] x Input vector (only values between 0 and 1 are relevant)
  		/// @param [in] coeff Coefficients vector
  		/// @return Output vector
  		static inline vector_t<T> Emu(vector_t<T> x,vector_t<T> coeff) {
  			vector_t<T> x3 = x*x*x;
  			return v1+x+x3*coeff;
  		}


  		/// @brief C0+C1*x+C2*x^2+C3*x^3
  		/// @param [in] x Variable
  		/// @param [in] C0 Zero order coefficient
  		/// @param [in] C1 First order coefficient
  		/// @param [in] C2 Second order coefficient
  		/// @param [in] C3 Third order coefficient
  		/// @return Result
  		static inline vector_t<T> UpTo3(const vector_t<T>& x,
  				const vector_t<T>& C0, const vector_t<T>& C1, const vector_t<T>& C2, const vector_t<T>& C3) {
  			vector_t<T> x2 = x*x;
  			return C0+C1*x+C2*x2+C3*x*x2;
  		}


  		/// @brief Derivative of : C1*x+C2*x^2+C3*x^3
  		/// @param [in] x Variable
  		/// @param [in] C1 First order coefficient
  		/// @param [in] C2 Second order coefficient
  		/// @param [in] C3 Third order coefficient
  		/// @return Result
  		static inline vector_t<T> Derivative3(const vector_t<T>& x,
  				const vector_t<T>& C1, const vector_t<T>& C2, const vector_t<T>& C3) {
  			return C1+v2*C2*x+v3*C3*x*x;
  		}

  		/// @brief C0+C1*x+C2*x^2+C3*x^3+C4*x^4+C5*x^5+C6*x^6
  		/// @param [in] x Variable
  		/// @param [in] C0 Zero order coefficient
  		/// @param [in] C1 First order coefficient
  		/// @param [in] C2 Second order coefficient
  		/// @param [in] C3 Third order coefficient
  		/// @param [in] C4 Fourth order coefficient
  		/// @param [in] C5 Fifth order coefficient
  		/// @param [in] C6 Sixth order coefficient
  		/// @return Result
  		static inline vector_t<T> UpTo6(const vector_t<T>& x,
  				const vector_t<T>& C0, const vector_t<T>& C1, const vector_t<T>& C2, const vector_t<T>& C3, const vector_t<T>& C4, const vector_t<T>& C5, const vector_t<T>& C6) {
  			vector_t<T> x2 = x*x;
  			vector_t<T> x3 = x*x2;
  			return C0+C1*x+C2*x2+C3*x3+C4*x2*x2+C5*x3*x2+C6*x3*x3;
  		}


  		/// @brief Derivative of : C1*x+C2*x^2+C3*x^3+C4*x^4+C5*x^5+C6*x^6
  		/// @param [in] x Variable
  		/// @param [in] C1 First order coefficient
  		/// @param [in] C2 Second order coefficient
  		/// @param [in] C3 Third order coefficient
  		/// @param [in] C4 Fourth order coefficient
  		/// @param [in] C5 Fifth order coefficient
  		/// @param [in] C6 Sixth order coefficient
  		/// @return Result
  		static inline vector_t<T> Derivative6(const vector_t<T>& x,
  				const vector_t<T>& C1, const vector_t<T>& C2, const vector_t<T>& C3, const vector_t<T>& C4, const vector_t<T>& C5, const vector_t<T>& C6) {
  			vector_t<T> x2 = x*x;
  			vector_t<T> x3 = x*x2;
  			return C1+v2*C2*x+v3*C3*x2+v4*C4*x3+v5*C5*x2*x2+v6*C6*x2*x3;
  		}


  	private:

  		static const vector_t<T> v1; ///< Vector of 1
  		static const vector_t<T> v2; ///< Vector of 2
    	static const vector_t<T> v3; ///< Vector of 3
    	static const vector_t<T> v4; ///< Vector of 4
    	static const vector_t<T> v5; ///< Vector of 4
    	static const vector_t<T> v6; ///< Vector of 6
  	};


  	template <class T> const vector_t<T> Polynoms<T>::v1  (  1);
  	template <class T> const vector_t<T> Polynoms<T>::v2  (  2);
  	template <class T> const vector_t<T> Polynoms<T>::v3  (  3);
  	template <class T> const vector_t<T> Polynoms<T>::v4  (  4);
  	template <class T> const vector_t<T> Polynoms<T>::v5  (  5);
  	template <class T> const vector_t<T> Polynoms<T>::v6  (  6);

  }

}


#endif /* POLYNOMS_HPP_ */
