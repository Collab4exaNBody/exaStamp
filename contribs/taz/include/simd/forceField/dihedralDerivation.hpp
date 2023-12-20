/*
 * dihedralDerivation.hpp
 *
 *  Created on: Mar 29, 2017
 *      Author: giarda
 */

/// @file
/// @brief Simd tools for the dihedral angle calculation and derivation

#ifndef DIHEDRALDERIVATION_HPP_
#define DIHEDRALDERIVATION_HPP_


#include "libevi/simd.hpp"

#include "simd/forceField/bondDerivation.hpp"
#include "simd/forceField/angleDerivation.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {


		/// @brief Cross product of 2 vectors
		/// @param [in] xa X component of the first vector
		/// @param [in] ya Y component of the first vector
		/// @param [in] za Z component of the first vector
		/// @param [in] xb X component of the second vector
		/// @param [in] yb Y component of the second vector
		/// @param [in] zb Z component of the second vector
		/// @param [out] xcross X component of the cross product
		/// @param [out] ycross Y component of the cross product
		/// @param [out] zcross Z component of the cross product
		template <class T> inline void cross(const vector_t<T>& xa, const vector_t<T>& ya, const vector_t<T>& za,
				const vector_t<T>& xb, const vector_t<T>& yb, const vector_t<T>& zb,
				vector_t<T>& xcross, vector_t<T>& ycross, vector_t<T>& zcross) {
			xcross=ya*zb-za*yb;
			ycross=za*xb-xa*zb;
			zcross=xa*yb-ya*xb;
		}


		/// @brief Compute the derivative of the cosine of the angle between A and B
		/// Case where B is a constant
		/// @param [in] xA X component of the vector A
		/// @param [in] yA Y component of the vector A
		/// @param [in] zA Z component of the vector A
		/// @param [in] xB X component of the vector B
		/// @param [in] yB Y component of the vector B
		/// @param [in] zB Z component of the vector B
		/// @param [in] dxA X component of the derivative of vector A
		/// @param [in] dyA Y component of the derivative of vector A
		/// @param [in] dzA Z component of the derivative of vector A
		/// @param [in] cos Cosine of the angle between A and B
		/// @return Cosine derivative
		template <class T> inline vector_t<T> derivativeCosine1Variable(const vector_t<T>& xA, const vector_t<T>& yA, const vector_t<T>& zA,
				const vector_t<T>& xB, const vector_t<T>& yB, const vector_t<T>& zB,
				const vector_t<T>& dxA, const vector_t<T>& dyA, const vector_t<T>& dzA,
				const vector_t<T>& cos) {
			vector_t<T> normA;
			normA=norm2(xA,yA,zA);
			return scalar(dxA,dyA,dzA,xB,yB,zB)/(normA*norm2(xB,yB,zB))-cos*scalar(dxA,dyA,dzA,xA,yA,zA)/(normA*normA);
		}


		/// @brief Compute the derivative of the cosine of the angle between A and B
		/// Case where B is not a constant
		/// @param [in] xA X component of the vector A
		/// @param [in] yA Y component of the vector A
		/// @param [in] zA Z component of the vector A
		/// @param [in] xB X component of the vector B
		/// @param [in] yB Y component of the vector B
		/// @param [in] zB Z component of the vector B
		/// @param [in] dxA X component of the derivative of vector A
		/// @param [in] dyA Y component of the derivative of vector A
		/// @param [in] dzA Z component of the derivative of vector A
		/// @param [in] dxB X component of the derivative of vector B
		/// @param [in] dyB Y component of the derivative of vector B
		/// @param [in] dzB Z component of the derivative of vector B
		/// @param [in] cos Cosine of the angle between A and B
		/// @return Cosine derivative
		template <class T> inline vector_t<T> derivativeCosine2Variables(const vector_t<T>& xA, const vector_t<T>& yA, const vector_t<T>& zA,
				const vector_t<T>& xB, const vector_t<T>& yB, const vector_t<T>& zB,
				const vector_t<T>& dxA, const vector_t<T>& dyA, const vector_t<T>& dzA,
				const vector_t<T>& dxB, const vector_t<T>& dyB, const vector_t<T>& dzB,
				const vector_t<T>& cos) {
			vector_t<T> normA, normB;
			normA=norm2(xA,yA,zA);
			normB=norm2(xB,yB,zB);
			return (scalar(dxA,dyA,dzA,xB,yB,zB)+scalar(xA,yA,zA,dxB,dyB,dzB))/(normA*normB)-cos*scalar(dxA,dyA,dzA,xA,yA,zA)/(normA*normA)-cos*scalar(dxB,dyB,dzB,xB,yB,zB)/(normB*normB);
		}


	}  // namespace kernels


}  // namespace simd

#endif /* DIHEDRALDERIVATION_HPP_ */
