/*
 * bharm.hpp
 *
 *  Created on: Oct 9, 2017
 *      Author: giarda
 */

/// @file
/// @brief Bond force computer with harmonic fonctionnal force

#ifndef BHARM_HPP_
#define BHARM_HPP_

#include "simd/forceField/bondDerivation.hpp"

/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {

		/// @brief Bond force computer with harmonic fonctionnal force
		template <class T> class bharm {

		public :

			void operator () (T* fx_, T* fy_, T* fz_, T* e_,
							const T* rx_, const T* ry_, const T* rz_,
							const T* k_, const T* rij_, const uint n) const;


		private :

    	static const vector_t<T> vm2; ///< Vector of -2

		};

		template <class T> const vector_t<T> bharm<T>::vm2 	( -2.);

		/// @brief Function call operator : compute bond forces and energies from the distances
		/// @param [out] fx_ X components of the forces
		/// @param [out] fy_ Y components of the forces
		/// @param [out] fz_ Z components of the forces
		/// @param [out] e_ Energies
		/// @param [in] rx_ X components of the distance vector
		/// @param [in] ry_ Y components of the distance vector
		/// @param [in] rz_ Z components of the distance vector
		/// @param [in] k_ Force constants
		/// @param [in] rij_ Equilibrium distances
		/// @param [in] n Number of elements
		template <class T> void bharm<T>::operator ()(T* fx_, T* fy_, T* fz_, T* e_,
				const T* rx_, const T* ry_, const T* rz_,
				const T* k_, const T* rij_, const uint n) const {
    	vector_t<T> x, y, z, rij, k, tmp0, tmp1;
    	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {
    		// Load data
    		x.load(rx_+i);
    		y.load(ry_+i);
    		z.load(rz_+i);
    		rij.load(rij_+i);
    		k.load(k_+i);
    		// Compute distance (tmp0) and distance derivative (in x2, y2 and z2)
    		tmp0=norm2(x,y,z);
    		norm2Derivatives(x,y,z,tmp0);
    		// Compute energy = k(r-rij)^2 (tmp1)
    		tmp0=tmp0-rij; 		// tmp0=r-rij
    		tmp1=k*tmp0*tmp0;
//    		// Compute forces
//    		// fx = -2*k (r-rij) dr/dx (x1)
//    		// fy = -2*k (r-rij) dr/dy (y1)
//    		// fz = -2*k (r-rij) dr/dz (z1)
    		tmp0=vm2*k*tmp0;	// tmp0=-2k(r-rij)
    		x=tmp0*x;
    		y=tmp0*y;
    		z=tmp0*z;
    		x.store(fx_+i);
    		y.store(fy_+i);
    		z.store(fz_+i);
    		tmp1.store(e_+i);
    	}
		}

	}  // namespace kernels


}  // namespace simd

#endif /* BHARM_HPP_ */
