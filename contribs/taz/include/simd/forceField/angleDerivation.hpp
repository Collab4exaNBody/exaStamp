/*
 * angleDerivation.hpp
 *
 *  Created on: Feb 23, 2017
 *      Author: giarda
 */

/// @file
/// @brief Simd tools for the angle calculation and derivation

#ifndef ANGLEDERIVATION_HPP_
#define ANGLEDERIVATION_HPP_


#include "libevi/simd.hpp"

#include "simd/forceField/bondDerivation.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {


		/// @brief Compute the scalar product of two distance vectors
		/// @param [in] rx1 X component of the first vector
		/// @param [in] ry1 Y component of the first vector
		/// @param [in] rz1 Z component of the first vector
		/// @param [in] rx2 X component of the second vector
		/// @param [in] ry2 Y component of the second vector
		/// @param [in] rz2 Z component of the second vector
		/// @return Scalar product
		template <class T> inline vector_t<T> scalar(const vector_t<T>& rx1, const vector_t<T>& ry1, const vector_t<T>& rz1,
				const vector_t<T>& rx2, const vector_t<T>& ry2, const vector_t<T>& rz2) {
			return rx1*rx2+ry1*ry2+rz1*rz2;
		}


		/// @brief Compute the cosine of the angle between two vectors
		/// @param [in] rx1 X component of the first vector
		/// @param [in] ry1 Y component of the first vector
		/// @param [in] rz1 Z component of the first vector
		/// @param [in] rx2 X component of the second vector
		/// @param [in] ry2 Y component of the second vector
		/// @param [in] rz2 Z component of the second vector
		/// @return Cosine of the angle
		template <class T> inline vector_t<T> cosine(const vector_t<T>& rx1, const vector_t<T>& ry1, const vector_t<T>& rz1,
				const vector_t<T>& rx2, const vector_t<T>& ry2, const vector_t<T>& rz2) {
			return scalar(rx1,ry1,rz1,rx2,ry2,rz2)/(norm2(rx1,ry1,rz1)*norm2(rx2,ry2,rz2));
		}


		/// @brief Compute the angle between two vectors
		/// @param [in] rx1 X component of the first vector
		/// @param [in] ry1 Y component of the first vector
		/// @param [in] rz1 Z component of the first vector
		/// @param [in] rx2 X component of the second vector
		/// @param [in] ry2 Y component of the second vector
		/// @param [in] rz2 Z component of the second vector
		/// @return Angle
		template <class T> inline vector_t<T> angle(const vector_t<T>& rx1, const vector_t<T>& ry1, const vector_t<T>& rz1,
				const vector_t<T>& rx2, const vector_t<T>& ry2, const vector_t<T>& rz2) {
			return acos(cosine(rx1, ry1, rz1, rx2, ry2, rz2));
		}


		/// @brief Compute the derivative of the cosine of the angle 123 relative to a coordinate of 1
		/// @param [in] r21a Component of the first vector relatively to which the derivation is done
		/// @param [in] r21b The second component of the first vector
		/// @param [in] r21c The third component of the first vector
		/// @param [in] r23a The first component of the second vector
		/// @param [in] r23b The second component of the second vector
		/// @param [in] r23c The third component of the second vector
		/// @return Angle cosine derivative to a coordinate of 1
		template <class T> inline vector_t<T> cosineDerivativeSide(const vector_t<T>& r21a, const vector_t<T>& r21b, const vector_t<T>& r21c,
				const vector_t<T>& r23a, const vector_t<T>& r23b, const vector_t<T>& r23c) {
			vector_t<T> cos, d12, d23;
			cos=cosine(r21a,r21b,r21c,r23a,r23b,r23c);
			d12=norm2(r21a,r21b,r21c);
			d23=norm2(r23a,r23b,r23c);
			return r21a/(d12*d12)*cos-r23a/(d12*d23);
		}


		/// @brief Compute the derivative of the cosine of the angle 123 relative to a coordinate of 2
		/// @param [in] r21a Component of the first vector relatively to which the derivation is done
		/// @param [in] r21b The second component of the first vector
		/// @param [in] r21c The third component of the first vector
		/// @param [in] r23a The first component of the second vector
		/// @param [in] r23b The second component of the second vector
		/// @param [in] r23c The third component of the second vector
		/// @return Angle cosine derivative to a coordinate of 2
		template <class T> inline vector_t<T> cosineDerivativeMiddle(const vector_t<T>& r21a, const vector_t<T>& r21b, const vector_t<T>& r21c,
				const vector_t<T>& r23a, const vector_t<T>& r23b, const vector_t<T>& r23c) {
			vector_t<T> cos, d12, d23;
			cos=cosine(r21a,r21b,r21c,r23a,r23b,r23c);
			d12=norm2(r21a,r21b,r21c);
			d23=norm2(r23a,r23b,r23c);
			return (r21a/(d12*d12)+r23a/(d23*d23))*cos-(r21a+r23a)/(d12*d23);
		}


		/// @brief Compute the opposite of the derivative of the angle 123 relative to a coordinate of 1
		/// @param [in] r21a Component of the first vector relatively to which the derivation is done
		/// @param [in] r21b The second component of the first vector
		/// @param [in] r21c The third component of the first vector
		/// @param [in] r23a The first component of the second vector
		/// @param [in] r23b The second component of the second vector
		/// @param [in] r23c The third component of the second vector
		/// @return Opposite angle derivative to a coordinate of 1
		template <class T> inline vector_t<T> mAngleDerivativeSide(const vector_t<T>& r21a, const vector_t<T>& r21b, const vector_t<T>& r21c,
				const vector_t<T>& r23a, const vector_t<T>& r23b, const vector_t<T>& r23c) {
			return cosineDerivativeSide(r21a,r21b,r21c,r23a,r23b,r23c)/sin(angle(r21a,r21b,r21c,r23a,r23b,r23c));
		}


		/// @brief Compute the opposite of the derivative of the angle 123 relative to a coordinate of 2
		/// @param [in] r21a Component of the first vector relatively to which the derivation is done
		/// @param [in] r21b The second component of the first vector
		/// @param [in] r21c The third component of the first vector
		/// @param [in] r23a The first component of the second vector
		/// @param [in] r23b The second component of the second vector
		/// @param [in] r23c The third component of the second vector
		/// @return Opposite angle derivative to a coordinate of 2
		template <class T> inline vector_t<T> mAngleDerivativeMiddle(const vector_t<T>& r21a, const vector_t<T>& r21b, const vector_t<T>& r21c,
				const vector_t<T>& r23a, const vector_t<T>& r23b, const vector_t<T>& r23c) {
			return cosineDerivativeMiddle(r21a,r21b,r21c,r23a,r23b,r23c)/sin(angle(r21a,r21b,r21c,r23a,r23b,r23c));
		}


	}  // namespace kernels


}  // namespace simd

#endif /* ANGLEDERIVATION_HPP_ */
