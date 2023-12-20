/// @file
/// @brief Vectorization of EAM VNIITF potential force computation

#ifndef __EAM_VNIITF_HPP_INCLUDED
#define __EAM_VNIITF_HPP_INCLUDED


#include "libevi/simd.hpp"

#include "simd/spline.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


    /// @brief Simd version for a EAM VNIITF potential
  	/// @tparam T Type of the variables
    template <class T> class EamVniitf {

    private:

      vector_t<T> half; ///< Vector of 1/2
      vector_t<T> minus_one; ///< Vector of -1
      vector_t<T> minus_beta; ///< Opposite of the attenuation rate
      vector_t<T> rt0; ///< Characteristic radius
      vector_t<T> irt0; ///< Inverse of the characteristic radius
      vector_t<T> iZ; ///< Inverse of the number of neighbors in the reference structure
      vector_t<T> rmax; ///< Maximal distance for density contribution
      vector_t<T> irdiff; ///< Inverse of the difference between the max and min distances for density contribution
      vector_t<T> A; ///< -2*Ecoh/Z (not empirical parameter A)
      vector_t<T> alpha; ///< Parameter alpha
      vector_t<T> minus_alpha; ///< Opposite of the parameter alpha
      vector_t<T> eta; ///< Parameter eta
      vector_t<T> mu; ///< Parameter mu
      vector_t<T> B; ///< alpha*alpha*alpha*D*rt0
      vector_t<T> E0; ///< Energy parameter
      vector_t<T> two_eta; ///< Two times parameter eta
      vector_t<T> three_mu; ///< Three times parameter mu
      vector_t<T> three_ialpha; ///< Three times the inverse of parameter alpha
      vector_t<T> rhoCut; ///< Density contribution at cutoff radius
      vector_t<T> phiCut; ///< \f$ \phi(r) \f$ term at cutoff radius

    public:

      /// @brief Set the parameters of the potential
      /// @tparam T Type of the variables
      /// @param [in] rmax_ Maximal distance for density contribution
      /// @param [in] rmin_ Minimal distance for density contribution
      /// @param [in] rt0_ Characteristic radius
      /// @param [in] Ecoh_ Cohesive energy
      /// @param [in] E0_ Energy parameter
      /// @param [in] beta_ Attenuation rate
      /// @param [in] Z_ Number of neighbors in the reference structure
      /// @param [in] alpha_ Parameter alpha
      /// @param [in] D_ Parameter D
      /// @param [in] eta_ Parameter eta
      /// @param [in] mu_ Parameter mu
      /// @param [in] rhoCut_ Density contribution at cutoff radius
      /// @param [in] phiCut_ \f$ \phi(r) \f$ term at cutoff radius
      inline void setParameters(const T& rmax_, const T& rmin_, const T& rt0_, const T& Ecoh_, const T& E0_, 
				const T& beta_, const T& Z_, const T& alpha_, const T& D_, const T& eta_, 
				const T& mu_, const T& rhoCut_, const T& phiCut_) {

      	half         = vector_t<T>(0.5);
      	minus_one    = vector_t<T>(-1.);

      	rt0          = vector_t<T>(rt0_);
      	irt0         = vector_t<T>(1./rt0_);
      	iZ           = vector_t<T>(1/Z_);
      	rmax         = vector_t<T>(rmax_);
      	irdiff       = vector_t<T>(1./(rmax_-rmin_));
      	A            = vector_t<T>(-2*Ecoh_/Z_);
      	alpha        = vector_t<T>(alpha_);
      	eta          = vector_t<T>(eta_);
      	mu           = vector_t<T>(mu_);
      	B            = vector_t<T>(alpha_*alpha_*alpha_*D_*rt0_);
      	E0           = vector_t<T>(E0_);
      	minus_alpha  = vector_t<T>(-1*alpha_);
      	minus_beta   = vector_t<T>(-1*beta_);
      	two_eta      = vector_t<T>(2*eta_);
      	three_mu     = vector_t<T>(3*mu_);
      	three_ialpha = vector_t<T>(3./alpha_);
      	rhoCut       = vector_t<T>(rhoCut_);
      	phiCut       = vector_t<T>(phiCut_);

      }

      /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$)
      /// @tparam T Type of the variables
      /// @param [out] out_ Density contribution
      /// @param [in] drx_ X components of the distances
      /// @param [in] dry_ Y components of the distances
      /// @param [in] drz_ Z components of the distances
      /// @param [in] n Number of elements
      void rho (T* out_, const T* drx_, const T* dry_, const T* drz_, const uint n) {

      	vector_t<T> tmp0, tmp1, tmp2, drx, dry, drz;

      	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

      		// Load data
      		drx.load(drx_+i);
      		dry.load(dry_+i);
      		drz.load(drz_+i);

      		// tmp0 = r
      		// tmp1 = - beta * (r/rt0 - 1)
      		// tmp2 = exp(tmp1)/Z
      		tmp0 = sqrt(drx*drx + dry*dry + drz*drz);
      		tmp1 = minus_beta * fmadd(tmp0, irt0, minus_one);
      		tmp2 = exp(tmp1)*iZ;

      		// tmp1 = (rmax-r) / (rmax-rmin)
      		// tmp0 = S<3>(tmp1)
      		tmp1 = (rmax-tmp0) * irdiff;
      		tmp0 = Spline<T>::S3(tmp1);

      		// tmp1 = rho
      		tmp1 = fmsub(tmp2, tmp0, rhoCut);

      		// Now store back the results
      		tmp1.store(out_+i);

      	}

      }

      /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and force
      /// @tparam T Type of the variables
      /// @param [out] fx_ X components to the forces
      /// @param [out] fy_ Y components to the forces
      /// @param [out] fz_ Z components to the forces
      /// @param [out] en_ Potential energies
      /// @param [in] drx_ X components of the distances
      /// @param [in] dry_ Y components of the distances
      /// @param [in] drz_ Z components of the distances
      /// @param [in] emb_ Embedding terms for the neighbors
      /// @param [in] e Embedding term for the particle
      /// @param [in] n Number of elements
      void phi (T* fx_, T* fy_, T* fz_, T* en_,
      		const T* drx_, const T* dry_, const T* drz_,
      		const T* emb_, const T& e, const uint n) {

      	vector_t<T> tmp0, tmp1, tmp2, tmp3, tmp4, drx, dry, drz;
      	vector_t<T> f1, df1, S, dS;

      	vector_t<T> emb(e);

      	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

      		// Load data
      		drx.load(drx_+i);
      		dry.load(dry_+i);
      		drz.load(drz_+i);

      		// tmp0 = r, r = [drx, dry, drz]
      		tmp0 = sqrt( drx*drx + dry*dry + drz*drz );

      		// tmp1 = 1/r
      		// tmp2 = dr  = (r-rt0)  / rt0
      		// tmp3 = drS = (rmax-r) / (rmax-rmin)
      		// tmp0 = dr^2
      		tmp1 = vector_t<T>::one()/tmp0;
      		tmp2 = fmadd(tmp0, irt0, minus_one);
      		tmp3 = (rmax-tmp0) * irdiff;
      		tmp0 = tmp2*tmp2;

      		// Pre-computing
      		f1  = A*( vector_t<T>::one() + alpha*tmp2 + eta*tmp0 + fmadd(B, tmp1, mu)*tmp0*tmp2 );
      		df1 = A*irt0*( alpha + two_eta*tmp2 + three_mu*tmp0 + B*(three_ialpha-rt0*tmp2*tmp1)*tmp0*tmp1 );
      		S   = Spline<T>:: S3(tmp3);
      		dS  = Spline<T>::dS3(tmp3) * minus_one * irdiff;

      		// tmp4 = f2
      		tmp4 = exp(minus_alpha*tmp2);
      		tmp0 = fmadd(f1, tmp4, E0);

      		// phi  = tmp0 * S;
      		// tmp3 = 0.5 * (phi - phiCut)
      		tmp3 = half * fmsub(S, tmp0, phiCut);
      		tmp3 . store(en_+i);

      		// tmp3 = dphi = tmp0 * dS + (f1*df2 + f2*df1) * S;
      		tmp3 = tmp0 * dS + (df1 - alpha*irt0*f1) * tmp4*S;

      		// tmp2 = drho
      		tmp0 = iZ*exp(minus_beta*tmp2);
      		tmp2 = tmp0 * (minus_beta*irt0*S + dS);

      		// tmp4 = de = drho(r) * (emb+emb[i]) + dphi(r)
      		// tmp3 = df = -de/r
      		tmp0 . load(emb_+i);
      		tmp4 = tmp2 * (tmp0 + emb) + tmp3;
      		tmp3 = minus_one * tmp4 * tmp1;

      		//
      		tmp0 = tmp3 * drx;
      		tmp1 = tmp3 * dry;
      		tmp2 = tmp3 * drz;

      		tmp0.store(fx_+i);
      		tmp1.store(fy_+i);
      		tmp2.store(fz_+i);

      	}

      }

    };


  }


}

#endif // __EAM_VNIITF_HPP_INCLUDED
