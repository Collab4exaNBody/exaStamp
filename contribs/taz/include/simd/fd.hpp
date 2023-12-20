/// @file
/// @brief Vectorization of fluctuation-dissipation

#ifndef __FD_HPP_INCLUDED
#define __FD_HPP_INCLUDED


#include "libevi/simd.hpp"


namespace simd {

  namespace kernels {



    /// @brief Handles the vectorized computation of fluctuation in dissipative dynamics
    /// @tparam T Numeric type
    template <class T>
    class DPDFluctuation {

    public:
      
      /// @brief Set parameters
      /// @param [in] ircut Inverse cut-off radius
      /// @param [in] ad1 Friction coefficient in the parallel direction
      /// @param [in] ad2 Friction coefficient in the orthoganal directions
      /// @param [in] time Timestep
      /// @param [in] beta Inverse temperature
      inline void setParameters(const T& ircut, const T& ad1, const T& ad2, const T& time, const T& beta = 1) {

	_ircut = vector_t<T>(ircut);
	_ad2   = vector_t<T>(-1.0*ad2);
	_af2   = vector_t<T>(auxSqrt(4.0*ad2/time/beta));
	_ad1m2 = vector_t<T>((ad2-ad1));
	_af1m2 = vector_t<T>(auxSqrt(4.0*ad1/time/beta))-_af2;

      }

      /// @brief Compute FD force (DPD)
      /// @param [out] fx_ Force in the x-direction
      /// @param [out] fy_ Force in the y-direction
      /// @param [out] fz_ Force in the z-direction
      /// @param [in] drx_ Distance in the x-direction
      /// @param [in] dry_ Distance in the y-direction
      /// @param [in] drz_ Distance in the z-direction
      /// @param [in] dvx_ Velocity difference in the x-direction
      /// @param [in] dvy_ Velocity difference in the y-direction
      /// @param [in] dvz_ Velocity difference in the z-direction
      /// @param [in] sgn Sign prefactor to ensure correct parallelisation
      /// @param [in] rndx0_ Random number
      /// @param [in] rndx1_ Random number
      /// @param [in] rndy0_ Random number
      /// @param [in] rndy1_ Random number
      /// @param [in] rndz1_ Random number
      /// @param [in] n Number of particles
      void operator () (T* fx_, T* fy_, T* fz_, 
			const T* drx_, const T* dry_, const T* drz_,
			const T* dvx_, const T* dvy_, const T* dvz_,
			const T* sgn, 
			const T* rndx0_, const T* rndx1_, const T* rndy0_, const T* rndy1_, const T* rndz1_, 
			const uint n) {

	vector_t<T> ex, ey, ez, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  ex.load(drx_+i);
	  ey.load(dry_+i);
	  ez.load(drz_+i);

	  // tmp0 = chi
	  tmp0 = inv_sqrt( ex*ex + ey*ey + ez*ez );  
	  ex = ex * tmp0;
	  ey = ey * tmp0;
	  ez = ez * tmp0;

	  tmp0 = vector_t<T>::one() - _ircut/tmp0;

	  // tmp4 = chi*sign
	  tmp4.load(sgn+i);
	  tmp4 = tmp4*tmp0; 

	  tmp1 = tmp4*random<T>::make_normal(rndx0_+i, rndx1_+i);
	  tmp2 = tmp4*random<T>::make_normal(rndy0_+i, rndy1_+i);
	  tmp3 = tmp4*random<T>::make_normal(rndx0_+i, rndz1_+i);


	  // tmp0 = ad1m2*chi*chi
	  // tmp8 = ad2*chi*chi
	  tmp8 = _ad2*tmp0*tmp0; 
	  tmp0 = _ad1m2*tmp0*tmp0; 

	  // tmp(5-7) = -g*v+s*G (para-ortho) 
	  // tmp(1-3) = -g*v+s*G (ortho)
	  tmp4.load(dvx_+i);
	  tmp5 = fmadd(tmp0,tmp4,_af1m2*tmp1);
	  tmp1 = fmadd(tmp8,tmp4,_af2*tmp1);

	  tmp4.load(dvy_+i);
	  tmp6 = fmadd(tmp0,tmp4,_af1m2*tmp2);
	  tmp2 = fmadd(tmp8,tmp4,_af2*tmp2);

	  tmp4.load(dvz_+i);
	  tmp7 = fmadd(tmp0,tmp4,_af1m2*tmp3);
	  tmp3 = fmadd(tmp8,tmp4,_af2*tmp3);

	  tmp0 = fmadd(ex,tmp5,fmadd(ey,tmp6,ez*tmp7));
	  tmp1 = fmadd(ex,tmp0,tmp1);
	  tmp2 = fmadd(ey,tmp0,tmp2);
	  tmp3 = fmadd(ez,tmp0,tmp3);

	  // store back the results
	  tmp1.store(fx_+i);
	  tmp2.store(fy_+i);
	  tmp3.store(fz_+i);

	}

      }

      /// @brief Compute FD force (DPDE)
      /// @param [out] fx_ Force in the x-direction
      /// @param [out] fy_ Force in the y-direction
      /// @param [out] fz_ Force in the z-direction
      /// @param [in] drx_ Distance in the x-direction
      /// @param [in] dry_ Distance in the y-direction
      /// @param [in] drz_ Distance in the z-direction
      /// @param [in] dvx_ Velocity difference in the x-direction
      /// @param [in] dvy_ Velocity difference in the y-direction
      /// @param [in] dvz_ Velocity difference in the z-direction
      /// @param [in] sgn Sign prefactor to ensure correct parallelisation
      /// @param [in] bij_ Mean inverse temperature \f[ \beta_{ij} \f]
      /// @param [in] rndx0_ Random number
      /// @param [in] rndx1_ Random number
      /// @param [in] rndy0_ Random number
      /// @param [in] rndy1_ Random number
      /// @param [in] rndz1_ Random number
      /// @param [in] n Number of particles
      void operator () (T* fx_, T* fy_, T* fz_, 
			const T* drx_, const T* dry_, const T* drz_,
			const T* dvx_, const T* dvy_, const T* dvz_,
			const T* sgn, const T* bij_,
			const T* rndx0_, const T* rndx1_, const T* rndy0_, const T* rndy1_, const T* rndz1_, 
			const uint n) {

	vector_t<T> ex, ey, ez, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  ex.load(drx_+i);
	  ey.load(dry_+i);
	  ez.load(drz_+i);

	  // tmp0 = chi
	  tmp0 = inv_sqrt( ex*ex + ey*ey + ez*ez );  
	  ex = ex * tmp0;
	  ey = ey * tmp0;
	  ez = ez * tmp0;

	  tmp0 = vector_t<T>::one() - _ircut/tmp0;

	  // tmp4 = chi*sign
	  tmp4.load(sgn+i);
	  tmp4 = tmp4*tmp0; 

	  tmp1 = tmp4*random<T>::make_normal(rndx0_+i, rndx1_+i);
	  tmp2 = tmp4*random<T>::make_normal(rndy0_+i, rndy1_+i);
	  tmp3 = tmp4*random<T>::make_normal(rndx0_+i, rndz1_+i);


	  // tmp8 = chi*chi*bij
	  tmp8.load(bij_+i);
	  tmp8 = tmp0*tmp0*tmp8; 

	  // tmp(5-7) = -g*v+s*G (para-ortho) 
	  // tmp(1-3) = -g*v+s*G (ortho)
	  tmp4.load(dvx_+i);
	  tmp4 = tmp8*tmp4;
	  tmp5 = fmadd(_ad1m2,tmp4,_af1m2*tmp1);
	  tmp1 = fmadd(_ad2,tmp4,_af2*tmp1);

	  tmp4.load(dvy_+i);
	  tmp4 = tmp8*tmp4;
	  tmp6 = fmadd(_ad1m2,tmp4,_af1m2*tmp2);
	  tmp2 = fmadd(_ad2,tmp4,_af2*tmp2);

	  tmp4.load(dvz_+i);
	  tmp4 = tmp8*tmp4;
	  tmp7 = fmadd(_ad1m2,tmp4,_af1m2*tmp3);
	  tmp3 = fmadd(_ad2,tmp4,_af2*tmp3);

	  tmp0 = ex*tmp5 + ey*tmp6 + ez*tmp7;
	  tmp1 = fmadd(ex,tmp0,tmp1);
	  tmp2 = fmadd(ey,tmp0,tmp2);
	  tmp3 = fmadd(ez,tmp0,tmp3);

	  // store back the results
	  tmp1.store(fx_+i);
	  tmp2.store(fy_+i);
	  tmp3.store(fz_+i);



	}

      }




      /// @brief Compute FD force (SDPD)
      /// @param [out] fx_ Force in the x-direction
      /// @param [out] fy_ Force in the y-direction
      /// @param [out] fz_ Force in the z-direction
      /// @param [in] drx_ Distance in the x-direction
      /// @param [in] dry_ Distance in the y-direction
      /// @param [in] drz_ Distance in the z-direction
      /// @param [in] dvx_ Velocity difference in the x-direction
      /// @param [in] dvy_ Velocity difference in the y-direction
      /// @param [in] dvz_ Velocity difference in the z-direction
      /// @param [in] sgn Sign prefactor to ensure correct parallelisation
      /// @param [in] bij_ Mean inverse temperature \f[ \beta_{ij} \f]
      /// @param [in] sij_ Fluctuation amplitude
      /// @param [in] chi2_ Squared cut-off function
      /// @param [in] rndx0_ Random number
      /// @param [in] rndx1_ Random number
      /// @param [in] rndy0_ Random number
      /// @param [in] rndy1_ Random number
      /// @param [in] rndz1_ Random number
      /// @param [in] n Number of particles
      void operator () (T* fx_, T* fy_, T* fz_, 
			const T* drx_, const T* dry_, const T* drz_,
			const T* dvx_, const T* dvy_, const T* dvz_,
			const T* sgn, const T* bij_, const T* sij_, const T* chi2_, 
			const T* rndx0_, const T* rndx1_, const T* rndy0_, const T* rndy1_, const T* rndz1_, 
			const uint n) {

	vector_t<T> ex, ey, ez, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  ex.load(drx_+i);
	  ey.load(dry_+i);
	  ez.load(drz_+i);

	  tmp0 = inv_sqrt( ex*ex + ey*ey + ez*ez );  
	  ex = ex * tmp0;
	  ey = ey * tmp0;
	  ez = ez * tmp0;

	  // tmp0 = chi
	  tmp0.load(chi2_+i);

	  // tmp4 = chi*sign*sij
	  tmp4.load(sgn+i);
	  tmp5.load(sij_+i);
	  tmp4 = tmp4*sqrt(tmp0*tmp5); 

	  tmp1 = tmp4*random<T>::make_normal(rndx0_+i, rndx1_+i);
	  tmp2 = tmp4*random<T>::make_normal(rndy0_+i, rndy1_+i);
	  tmp3 = tmp4*random<T>::make_normal(rndx0_+i, rndz1_+i);

	  // tmp8 = chi*chi*(1-bij)
	  tmp8.load(bij_+i);
	  tmp8 = tmp0*(vector_t<T>::one()-tmp8); 

	  // tmp(5-7) = -g*v+s*G (para-ortho) 
	  // tmp(1-3) = -g*v+s*G (ortho)
	  tmp4.load(dvx_+i);
	  tmp4 = tmp8*tmp4;
	  tmp5 = fmadd(_ad1m2,tmp4,_af1m2*tmp1);
	  tmp1 = fmadd(_ad2,tmp4,_af2*tmp1);

	  tmp4.load(dvy_+i);
	  tmp4 = tmp8*tmp4;
	  tmp6 = fmadd(_ad1m2,tmp4,_af1m2*tmp2);
	  tmp2 = fmadd(_ad2,tmp4,_af2*tmp2);

	  tmp4.load(dvz_+i);
	  tmp4 = tmp8*tmp4;
	  tmp7 = fmadd(_ad1m2,tmp4,_af1m2*tmp3);
	  tmp3 = fmadd(_ad2,tmp4,_af2*tmp3);

	  tmp0 = ex*tmp5 + ey*tmp6 + ez*tmp7;
	  tmp1 = fmadd(ex,tmp0,tmp1);
	  tmp2 = fmadd(ey,tmp0,tmp2);
	  tmp3 = fmadd(ez,tmp0,tmp3);

	  // store back the results
	  tmp1.store(fx_+i);
	  tmp2.store(fy_+i);
	  tmp3.store(fz_+i);

	}

      }

    private:

      vector_t<T> _ad1m2; ///< Difference of the dissipation amplitudes \f[ A_d^1 - A_d^2 \f]
      vector_t<T> _af1m2; ///< Difference of the fluctuation amplitudes \f[ A_f^1 - A_f^2 \f]
      vector_t<T> _ad2; ///< Amplitude of dissipation in the orthogonal directions 
      vector_t<T> _af2; ///< Amplitude of fluctuation in the orthogonal directions 
      vector_t<T> _ircut; ///< Inverse cut-off radius

    };



  }

}



#endif // __FD_HPP_INCLUDED
