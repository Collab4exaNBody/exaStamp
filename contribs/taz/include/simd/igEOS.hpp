/// @file
/// @brief Vectorization of ideal gas equation of state

#ifndef __IG_EOS_HPP_INCLUDED
#define __IG_EOS_INCLUDED


#include "libevi/simd.hpp"

#include "utils/stampUnits.hpp"


namespace simd {

  namespace kernels {



    /// @brief Compute thermodynamical quantities for ideal gas equation of state
    /// @tparam T Numeric type
    template <class T> class IdealGas {

    private:

      vector_t<T> _cv; ///< Heat capacity
      vector_t<T> _cvk; ///< Heat capacity (in \f[ k_B \f])
      vector_t<T> _icv; ///< Inverse heat capacity
      vector_t<T> _icvk; ///< Inverse heat capacity (in \f[ \frac1{k_B} \f])
      vector_t<T> _cd; ///< Coefficient of the density contribution to entropy
      vector_t<T> _cp; ///< Constant for pressure



    public:

      /// @brief Set parameters
      /// @param [in] size Mesoparticle size
      /// @param [in] mass Mesoparticle mass
      inline void setParameters(const uint& size, const T& mass) {
	
	_cv          = vector_t<T>(1.5*(size-1)*Stamp_Constant::boltzmann);
	_cvk         = vector_t<T>(1.5*(size-1));
	_icv         = vector_t<T>(1./(1.5*(size-1)*Stamp_Constant::boltzmann));
	_icvk        = vector_t<T>(1./(1.5*(size-1)));
	_cd          = vector_t<T>((size-1)*Stamp_Constant::boltzmann);
	_cp          = vector_t<T>(2./3./mass);
	
      }



      /// @brief Compute entropy, temperature and pressure from density and energy
      /// @param [out] s_ Entropy
      /// @param [out] b_ Inverse temperature
      /// @param [out] p_ Pressure
      /// @param [out] ic_ Inverse heat capacity
      /// @param [in] d_ Density
      /// @param [in] e_ Energy
      /// @param [in] n Number of particles
      void SBP (T* s_, T* b_, T* p_, T* ic_, const T* d_, const T* e_, const uint n) {

	vector_t<T> ss, bb, pp, dd, ee;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ee.load(e_+i);

	  // We want to compute (ideal gas) : 
	  //  s =  cv*log(e) - cd*log(d);
	  //  t =  e*icv;
	  //  p =  d*e*cp

	  ss = _cv*log(ee) - _cd*log(dd);
	  bb = _cvk/ee;
	  pp = dd*ee*_cp;

	  // store back the results
	  ss .store(s_+i);
	  bb .store(b_+i);
	  pp .store(p_+i);
	  _icvk .store(ic_+i);

	}

      }




      /// @brief Compute energy, inverse temperature and pressure from density and entropy
      /// @param [out] e_ Energy
      /// @param [out] b_ Inverse temperature
      /// @param [out] p_ Pressure
      /// @param [out] ic_ Inverse heat capacity
      /// @param [in] d_ Density
      /// @param [in] s_ Entropy
      /// @param [in] n Number of particles
      void EBP (T* e_, T* b_, T* p_, T* ic_, const T* d_, const T* s_, const uint n) {

	vector_t<T> ee, bb, pp, dd, ss;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ss.load(s_+i);

	  // We want to compute (ideal gas) : 
	  //  e =  exp(s/cv)*pow(d,2/3);
	  //  t =  e*icv;
	  //  p =  d*e*cp

	  ee = exp(ss*_icv)*cbrt(dd*dd);
	  bb = _cvk/ee;
	  pp = dd*ee*_cp;

	  // store back the results
	  ee .store(e_+i);
	  bb .store(b_+i);
	  pp .store(p_+i);
	  _icvk .store(ic_+i);

	}

      }



      /// @brief Compute only energy from density and entropy
      /// @param [out] e_ Energy
      /// @param[in] d_ Density
      /// @param[in] s_ entropy
      /// @param[in] n Number of particles
      void energy (T* e_, const T* d_, const T* s_, const uint n) {

	vector_t<T> ee, dd, ss;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ss.load(s_+i);

	  // We want to compute (ideal gas) : 
	  //  e =  exp(s/cv)*pow(d,2/3);

	  ee = exp(ss*_icv)*cbrt(dd*dd);

	  // store back the results
	  ee .store(e_+i);

	}

      }

    };



  }

}



#endif // __IG_EOS_HPP_INCLUDED
