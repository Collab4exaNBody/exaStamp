/// @file
/// @brief Vectorization of Lennard-Jones potential force computation

#ifndef __LJ_HPP_INCLUDED
#define __LJ_HPP_INCLUDED


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


    /// @brief Simd version for a Lennard-Jones potential
    /// @tparam T Type of the variables
    template <class T> class LennardJones {

    private:

      vector_t<T> _sigma2; ///< Squared finite distance at which the inter-particle potential is zero
      vector_t<T> _2epsilon; ///< Two times the depth of the potential well
      vector_t<T> _24epsilon; ///< Twenty four times the depth of the potential well
      vector_t<T> _minus_half_ecut; ///< Half of the energy at cutoff radius

    public:

      /// @brief Set the parameters of the potential
      /// @tparam T Type of the variables
      /// @param [in] sigma Finite distance at which the inter-particle potential is zero
      /// @param [in] epsilon Depth of the potential well
      /// @param [in] ecut Energy at cutoff radius
      inline void setParameters(const T& sigma, const T& epsilon, const T& ecut) {
	
      	_sigma2          = vector_t<T>(sigma*sigma);
      	_2epsilon        = vector_t<T>(2*epsilon);
      	_24epsilon       = vector_t<T>(24*epsilon);
      	_minus_half_ecut = vector_t<T>(-0.5*ecut);
	
      }

      /// @brief Function call operator : calculate the force and energy for the potential
      /// @tparam T Type of the variables
      /// @param [out] fx_ X components of the forces
      /// @param [out] fy_ Y components of the forces
      /// @param [out] fz_ Z components of the forces
      /// @param [out] en_ Potential energies
      /// @param [in] drx_ X components of the distances
      /// @param [in] dry_ Y components of the distances
      /// @param [in] drz_ Z components of the distances
      /// @param [in] n Number of elements
      void operator () (T* fx_, T* fy_, T* fz_, T* en_,
			const T* drx_, const T* dry_, const T* drz_, const uint n) {

      	vector_t<T> drx, dry, drz, tmp0, tmp1, tmp2;

      	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

      		// Load data
      		drx.load(drx_+i);
      		dry.load(dry_+i);
      		drz.load(drz_+i);

      		// We want to compute (lennard-jones) :
      		//  e =     4*epsilon * (   (sigma/r)^12 - (sigma/r)^6 ) - ecut;
      		// de =  - 24*epsilon * ( 2*(sigma/r)^12 - (sigma/r)^6 )/r;

      		// We will compute (more usefull for force computation) :
      		//  e' =   e/2
      		// de' = -de/r

      		// tmp0 = 1/r^2
      		tmp0 = vector_t<T>::one() / (drx*drx + dry*dry + drz*drz);

      		// tmp1 = (sigma/r)^2
      		// tmp2 = (sigma/r)^6
      		// tmp1 = (sigma/r)^12
      		tmp1 = _sigma2 * tmp0;
      		tmp2 = tmp1  * tmp1 * tmp1;
      		tmp1 = tmp2  * tmp2;

      		// let R2 = 1/r^2, R6 = (sigma/r)^6, R12 = (sigma/r)^12

      		// tmp2 =   R12 - R6
      		// tmp1 = 2*R12 - R6
      		tmp2 = tmp1 - tmp2;
      		tmp1 = tmp2 + tmp1;

      		// tmp2 =  2*epsilon*(R12 - R6) - ecut
      		// tmp1 = 24*epsilon*(2*R12 - R6)/R2
      		tmp2 = fmadd(_2epsilon, tmp2, _minus_half_ecut);
      		tmp1 = _24epsilon * tmp1 * tmp0;

      		// So we have tmp[3]=e' and tmp4=de'
      		// We want df = de' * dr
      		drx = drx * tmp1;
      		dry = dry * tmp1;
      		drz = drz * tmp1;

      		// Store back the results
      		drx .store(fx_+i);
      		dry .store(fy_+i);
      		drz .store(fz_+i);
      		tmp2.store(en_+i);

      	}

      }

    };



  }


}

#endif // __LJ_HPP_INCLUDED
