/*
 * acosto3.hpp
 *
 *  Created on: Oct 9, 2017
 *      Author: giarda
 */

/// @file
/// @brief Angle force computer where E=C0+C1*costheta+C2*cos2theta+C3*cos3theta

#ifndef ACOSTO3_HPP_
#define ACOSTO3_HPP_

#include "simd/forceField/angleDerivation.hpp"
#include "simd/polynoms.hpp"

/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {

		/// @brief Angle force computer where E=C0+C1*costheta+C2*cos2theta+C3*cos3theta
		template <class T> class acosto3{

		public :

			void operator() (T* fx1_, T* fy1_, T* fz1_, T* fx2_, T* fy2_, T* fz2_, T* e_,
					const T* rx21_, const T* ry21_, const T* rz21_, const T* rx23_, const T* ry23_, const T* rz23_,
					const T* c0_, const T* c1_, const T* c2_, const T* n_, const uint size) const;


		private :

  		static const vector_t<T> v2; ///< Vector of 2
    	static const vector_t<T> v3; ///< Vector of 3
    	static const vector_t<T> v4; ///< Vector of 4

		};

		template <class T> const vector_t<T> acosto3<T>::v2  	(  2.);
		template <class T> const vector_t<T> acosto3<T>::v3  	(  3.);
		template <class T> const vector_t<T> acosto3<T>::v4  	(  4.);


		/// @brief Function call operator : compute the angle forces and energies from the distances
		///
		/// Interaction between three particles 1, 2 and 3, where 2 is the central one
		/// @param [out] fx1_ X components of the forces between 1 and 2
		/// @param [out] fy1_ Y components of the forces between 1 and 2
		/// @param [out] fz1_ Z components of the forces between 1 and 2
		/// @param [out] fx2_ X components of the forces between 3 and 2
		/// @param [out] fy2_ Y components of the forces between 3 and 2
		/// @param [out] fz2_ Z components of the forces between 3 and 2
		/// @param [out] e_ Energies
		/// @param [in] rx21_ X components of the distance vector between 1 and 2
		/// @param [in] ry21_ Y components of the distance vector between 1 and 2
		/// @param [in] rz21_ Z components of the distance vector between 1 and 2
		/// @param [in] rx23_ X components of the distance vector between 3 and 2
		/// @param [in] ry23_ Y components of the distance vector between 3 and 2
		/// @param [in] rz23_ Z components of the distance vector between 3 and 2
		/// @param [in] c0_ Constants 0
		/// @param [in] c1_ Constants 1
		/// @param [in] c2_ Constants 2
		/// @param [in] c3_ Constants 3
		/// @param [in] size Number of elements
		template <class T> void acosto3<T>::operator() (T* fx1_, T* fy1_, T* fz1_, T* fx2_, T* fy2_, T* fz2_, T* e_,
				const T* rx21_, const T* ry21_, const T* rz21_, const T* rx23_, const T* ry23_, const T* rz23_,
				const T* c0_, const T* c1_, const T* c2_, const T* c3_, const uint size) const {
    	vector_t<T> rx21, ry21, rz21, rx23, ry23, rz23, e;
    	vector_t<T> costheta, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5;
    	for (uint i=0; i<size; i+=vector_t<T>::chunk_size) {
    		// Load data
    		rx21.load(rx21_+i);
    		ry21.load(ry21_+i);
    		rz21.load(rz21_+i);
    		rx23.load(rx23_+i);
    		ry23.load(ry23_+i);
    		rz23.load(rz23_+i);
    		tmp0.load(c0_+i);
    		tmp1.load(c1_+i);
    		tmp2.load(c2_+i);
    		tmp3.load(c3_+i);
    		// E=c0+c1*costheta+c2*cos2theta+c3*cos3theta
    		// with cos2theta and cos3theta are polynoms of costheta
    		costheta=cosine(rx21,ry21,rz21,rx23,ry23,rz23);
    		// Coefficients stored in tmp0, tmp1, tmp2, tmp3 and tmp4
    		tmp0=tmp0+tmp2;
    		tmp1=tmp1+v3*tmp3;
    		tmp2=v2*tmp2;
    		tmp3=v4*tmp3;
    		// Energy stored in n
    		e=Polynoms<T>::UpTo3(costheta, tmp0, tmp1, tmp2, tmp3);
    		// fx=dE/dcosx*dcosx/dx
    		// Energy derivative stored in costheta
    		costheta=Polynoms<T>::Derivative3(costheta, tmp1, tmp2, tmp3);
    		// Forces stored in tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, rz23;
    		tmp0=cosineDerivativeSide(rx21,ry21,rz21,rx23,ry23,rz23)*costheta;
    		tmp1=cosineDerivativeSide(ry21,rx21,rz21,ry23,rx23,rz23)*costheta;
    		tmp2=cosineDerivativeSide(rz21,ry21,rx21,rz23,ry23,rx23)*costheta;
    		tmp3=cosineDerivativeSide(rx23,ry23,rz23,rx21,ry21,rz21)*costheta;
    		tmp4=cosineDerivativeSide(ry23,rx23,rz23,ry21,rx21,rz21)*costheta;
    		tmp5=cosineDerivativeSide(rz23,ry23,rx23,rz21,ry21,rx21)*costheta;
    		// Store result
    		tmp0.store(fx1_+i);
    		tmp1.store(fy1_+i);
    		tmp2.store(fz1_+i);
    		tmp3.store(fx2_+i);
    		tmp4.store(fy2_+i);
    		tmp5.store(fz2_+i);
    		e.store(e_+i);
    	}
		}

	}  // namespace kernels


}  // namespace simd


#endif /* ACOSTO3_HPP_ */
