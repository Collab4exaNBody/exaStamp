/// @file
/// @brief Vectorization of Gaussian potential force computation

#ifndef __GPOT_HPP_INCLUDED
#define __GPOT_HPP_INCLUDED


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


    /// @brief Simd version for a Gaussian potential
  	/// @tparam T Type of the variables
    template <class T> class Gaussian {

    private:

      vector_t<T> _r2Att; ///< Squared attractive radius
      vector_t<T> _minus_r2Rep; ///< Squared repulsive radius
      vector_t<T> _minus_half_ir2Att; ///< Inverse Squared attractive radius
      vector_t<T> _minus_half_ir2Rep; ///< Inverse Squared repulsive radius
      vector_t<T> _ratio; ///< Attractive/repulsive ratio
      vector_t<T> _epsilon; ///< force amplitude
      vector_t<T> _minus_half_ecut; ///< Half of the energy at cutoff radius

    public:

      /// @brief Set the parameters of the potential
      /// @tparam T Type of the variables
      /// @param [in] ra2 Squared attractive radius
      /// @param [in] rb2 Squared repulsive radius
      /// @param [in] A Attractive/repulsive ratio
      /// @param [in] e Epsilon
      /// @param [in] ecut Energy at cutoff radius
      inline void setParameters(const T& ra2, const T& rb2, const T& A, const T& e, const T& ecut) {
	
      	_r2Att              = vector_t<T>(ra2);
	_minus_r2Rep        = vector_t<T>(-rb2);
	_minus_half_ir2Att  = vector_t<T>(-0.5/ra2);
	_minus_half_ir2Rep  = vector_t<T>(-0.5/rb2);
	_ratio              = vector_t<T>(A);
      	_epsilon            = vector_t<T>(e);
      	_minus_half_ecut    = vector_t<T>(-0.5*ecut);
	
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

      		// We want to compute (gaussian) :
      		//  e =     epsilon * ( A*ra^2exp(-0.5*r^2/ra^2) - rb^2*exp(-0.5*r^2/rb^2) ) - ecut;
      		// de =     epsilon * ( A*exp(-0.5*r^2/ra^2) - exp(-0.5*r^2/rb^2) );

      		// We will compute (more usefull for force computation) :
      		//  e' =   e/2
      		// de' = -de/r

      		// tmp0 = r^2
      		tmp0 = (drx*drx + dry*dry + drz*drz);

      		// tmp1 = A*exp(-0.5*r^2/ra^2)
      		// tmp2 = exp(-0.5*r^2/rb^2)
      		tmp1 = _ratio*exp(tmp0*_minus_half_ir2Att);
      		tmp2 = exp(tmp0*_minus_half_ir2Rep);

      		// tmp0 =  e'
      		// tmp1 = de'
      		tmp0 = _epsilon*fmadd(_r2Att,tmp1,_minus_r2Rep*tmp2);
		tmp1 = _epsilon*(tmp1-tmp2);;

      		// So we have tmp[3]=e' and tmp4=de'
      		// We want df = de' * dr
      		drx = drx * tmp1;
      		dry = dry * tmp1;
      		drz = drz * tmp1;

      		// Store back the results
      		drx .store(fx_+i);
      		dry .store(fy_+i);
      		drz .store(fz_+i);
      		tmp0.store(en_+i);

      	}

      }

    };



  }


}

#endif // __GPOT_HPP_INCLUDED
