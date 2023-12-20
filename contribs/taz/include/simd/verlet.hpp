/// @file
/// @brief Vectorization of Verlet list condition
#ifndef __VERLET_HPP_INCLUDED
#define __VERLET_HPP_INCLUDED

#include "libevi/simd.hpp"

/// @brief Namespace for SIMD implementations 
namespace simd {


	/// @brief Namespace for vectorized computations
	namespace kernels {


		/// @brief Vectorized test of verlet list condition (An atom must not move more than 1/2 of the verlet radius)
		/// @tparam T Type of the variables
		/// @param [in] rx_ X components of the old positions of the particles
		/// @param [in] ry_ Y components of the old positions of the particles
		/// @param [in] rz_ Z components of the old positions of the particles
		/// @param [in] ox_ X components of the new positions of the particles
		/// @param [in] oy_ Y components of the new positions of the particles
		/// @param [in] oz_ Z components of the new positions of the particles
		/// @param [out] out_ Boolean results
		/// @param [in] n Number of potential neighbors
	  /// @param [in] check_ value of (1/2 of the verlet radius)²
		/// @warning The boolean results are stored into a array of T (double)
		template <class T>
			inline bool Verlet ( const T *rx_, const T *ry_, const T *rz_, const T *ox_, const T *oy_, const T *oz_, T * out_, const uint n, T check_) {

                if(n==0) return false;

				
				
        #pragma omp simd
        ////#pragma vector aligned
				for (size_t j=0; j<n; ++j) {
				  const double rx = rx_[j] - ox_[j];
					const double ry = ry_[j] - oy_[j];
					const double rz = rz_[j] - oz_[j];
					out_[j] = (rx*rx + ry*ry + rz*rz) >= check_;
        }
                    
				return *std::max_element(out_,out_+n); // rewrite later
			}




		/// @brief Vectorized test of verlet list condition (An atom must not move more than 1/2 of the verlet radius)
		/// @tparam T Type of the variables
		/// @param [in] rx_ X components of the positions of the particles
		/// @param [in] ry_ Y components of the positions of the particles
		/// @param [in] rz_ Z components of the positions of the particles
		/// @param [in] rcut Value of the raduis cut-off
		/// @param [out] out_ Boolean results
		/// @param [in] n Number of potential neighbors
		/// @warning The boolean results are stored into a array of T (double)
		template <class T>
			inline void Verlet_corrector ( const T *rx_, const T *ry_, const T *rz_, const T rcut_,  T *out_, const uint n) {
			
        #pragma omp simd
        ////#pragma vector aligned
				for (size_t j=0; j<n; ++j) 
					out_[j] = ( (rx_[j]*rx_[j] + ry_[j]*ry_[j] + rz_[j]*rz_[j]) <= rcut_ );

			}

	}
}


#endif
