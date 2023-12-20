/*
 * bondDerivation.hpp
 *
 *  Created on: Feb 13, 2017
 *      Author: giarda
 */

/// @file
/// @brief Simd tools for the bond calculation and derivation

#ifndef BONDDERIVATION_HPP_
#define BONDDERIVATION_HPP_


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {


		/// @brief Compute the 2-norm of a 3D vector
		/// @param [in] x X component of the vector
		/// @param [in] y Y component of the vector
		/// @param [in] z Z component of the vector
		/// @return Norm
		template <class T> inline vector_t<T> norm2(const vector_t<T>& x, const vector_t<T>& y, const vector_t<T>& z) {
			return sqrt(x*x+y*y+z*z);
		}


		/// @brief Compute the partial derivatives of the 2-norm of a 3D vector
		/// @param [in,out] x X component of the vector then partial derivative with respect to x
		/// @param [in,out] y Y component of the vector then partial derivative with respect to y
		/// @param [in,out] z Z component of the vector then partial derivative with respect to z
		/// @param [in] r 2-norm of the vector
		/// @param [out] r Distance
		template <class T> inline void norm2Derivatives(vector_t<T>& x, vector_t<T>& y, vector_t<T>& z, vector_t<T>& r) {
			x=x/r;
			y=y/r;
			z=z/r;
		}


	}  // namespace kernels


}  // namespace simd

#endif /* BONDDERIVATION_HPP_ */
