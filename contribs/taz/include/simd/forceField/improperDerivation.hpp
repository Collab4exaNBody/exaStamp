/*
 * improperDerivation.hpp
 *
 *  Created on: May 2, 2017
 *      Author: giarda
 */

/// @file
/// @brief Simd tools for the improper torsion calculation and derivation

#ifndef IMPROPERDERIVATION_HPP_
#define IMPROPERDERIVATION_HPP_


#include "libevi/simd.hpp"

#include "simd/forceField/bondDerivation.hpp"
#include "simd/forceField/angleDerivation.hpp"
#include "simd/forceField/dihedralDerivation.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {


		/// @brief Normalize a vector to get a unit vector
		/// @param [in,out] x X component of the vector
		/// @param [in,out] y Y component of the vector
		/// @param [in,out] z Z component of the vector
		template <class T> inline void unit(vector_t<T>& x, vector_t<T>& y, vector_t<T>& z) {

			vector_t<T> norm;

			norm=norm2(x,y,z);

			x=x/norm;
			y=y/norm;
			z=z/norm;
		}


		/// @brief Compute the derivative of the normalized version of a vector
		/// @param [in] x X component of the vector
		/// @param [in] y Y component of the vector
		/// @param [in] z Z component of the vector
		/// @param [in,out] dx X component of the derivative of the vector, then of the derivative of the unit vector
		/// @param [in,out] dy Y component of the derivative of the vector, then of the derivative of the unit vector
		/// @param [in,out] dz Z component of the derivative of the vector, then of the derivative of the unit vector
		template <class T> inline void derivativeUnit(const vector_t<T>& x, const vector_t<T>& y, const vector_t<T>& z,
				vector_t<T>& dx, vector_t<T>& dy, vector_t<T>& dz) {

			vector_t<T> s, norm;

			s=scalar(dx,dy,dz,x,y,z);
			norm=norm2(x,y,z);

			dx=dx/norm-s*x/(norm*norm*norm);
			dy=dy/norm-s*y/(norm*norm*norm);
			dz=dz/norm-s*z/(norm*norm*norm);
		}


		/// @brief Compute the projection of a vector on an unit vector
		/// @param [in] x X component of the vector
		/// @param [in] y Y component of the vector
		/// @param [in] z Z component of the vector
		/// @param [in] xu X component of the unit vector
		/// @param [in] yu Y component of the unit vector
		/// @param [in] zu Z component of the unit vector
		/// @param [out] xproj X component of the projection
		/// @param [out] yproj Y component of the projection
		/// @param [out] zproj Z component of the projection
		template <class T> inline void proj(const vector_t<T>& x, const vector_t<T>& y, const vector_t<T>& z,
				const vector_t<T>& xu, const vector_t<T>& yu, const vector_t<T>& zu,
				vector_t<T>& xproj, vector_t<T>& yproj, vector_t<T>& zproj) {

			vector_t<T> s;

			s=scalar(x,y,z,xu,yu,zu);

			xproj=xu*s;
			yproj=yu*s;
			zproj=zu*s;
		}


	}  // namespace kernels


}  // namespace simd


#endif /* IMPROPERDERIVATION_HPP_ */
