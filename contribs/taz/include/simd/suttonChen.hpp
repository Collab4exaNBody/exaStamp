/// @file
/// @brief Vectorization of Sutton-Chen potential force computation

#ifndef __SUTTON_CHEN_HPP_INCLUDED
#define __SUTTON_CHEN_HPP_INCLUDED


#include "libevi/simd.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


    /// @brief Simd version for a Sutton-Chen potential
  	/// @tparam T Type of the variables
    template <class T> class SuttonChen {
      
    private:

      vector_t<T> half; ///< Vector of 1/2
      vector_t<T> minus_one; ///< Vector of -1
      vector_t<T> epsilon; ///< Energy scale
      vector_t<T> a0; ///< Length scale
      vector_t<T> m; ///< Shape parameter m
      vector_t<T> n; ///< Shape parameter n
      vector_t<T> rhoCut; ///< Density contribution at cutoff radius
      vector_t<T> embCut; ///< \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term at cutoff radius
      vector_t<T> phiCut; ///< \f$ \phi(r) \f$ term at cutoff radius
      
    public:

      /// @brief Set the parameters of the potential
      /// @tparam T Type of the variables
      /// @param [in] epsilon_ Energy scale
      /// @param [in] a0_ Length scale
      /// @param [in] m_ Shape parameter m
      /// @param [in] n_ Shape parameter n
      /// @param [in] rhoCut_ Density contribution at cutoff radius
      /// @param [in] embCut_ \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term at cutoff radius
      /// @param [in] phiCut_ \f$ \phi(r) \f$ term at cutoff radius
      inline void setParameters(const T& epsilon_, const T& a0_, const T& m_, const T& n_,
				const T& rhoCut_, const T& embCut_, const T& phiCut_) {
	
      	half      = vector_t<T>(0.5);
      	minus_one = vector_t<T>(-1.0);

      	epsilon   = vector_t<T>(epsilon_);
      	a0        = vector_t<T>(a0_);
      	m         = vector_t<T>(m_);
      	n         = vector_t<T>(n_);
	
      	rhoCut    = vector_t<T>(rhoCut_);
      	embCut    = vector_t<T>(embCut_);
      	phiCut    = vector_t<T>(phiCut_);
	
      }

      /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$)
      /// @tparam T Type of the variables
      /// @param [out] out_ Density contribution
      /// @param [in] drx_ X components of the distances
      /// @param [in] dry_ Y components of the distances
      /// @param [in] drz_ Z components of the distances
      /// @param [in] N Number of elements
      void rho (T* out_, const T* drx_, const T* dry_, const T* drz_, const uint N) {

      	vector_t<T> half_m = half * m;
      	vector_t<T> a0_2   = a0 * a0 ;
    
      	vector_t<T> tmp0, tmp1, drx, dry, drz;

      	for (uint i=0; i<N; i+=vector_t<T>::chunk_size) {

      		// Load data
      		drx.load(drx_+i);
      		dry.load(dry_+i);
      		drz.load(drz_+i);

      		// Compute a0/norm(r), r = [drx, dry, drz]
      		// tmp0 = r^2
      		// tmp1 = (a0/r)^m
      		tmp0 = drx*drx + dry*dry + drz*drz;
      		tmp1 = pow(a0_2/tmp0, half_m) - rhoCut;

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
      /// @param [in] N Number of elements
      void phi (T* fx_, T* fy_, T* fz_, T* en_,
      		const T* drx_, const T* dry_, const T* drz_,
      		const T* emb_, const T& e, const uint N) {

      	vector_t<T> emb(e);
      	vector_t<T> minus_m = m * minus_one;
      	vector_t<T> minus_n = n * minus_one;

      	vector_t<T> tmp0, tmp1, tmp2, tmp3, drx, dry, drz;

      	for (uint i=0; i<N; i+=vector_t<T>::chunk_size) {

      		// Load data
      		drx.load(drx_+i);
      		dry.load(dry_+i);
      		drz.load(drz_+i);

      		// tmp0 = 1/r, r = [drx, dry, drz]
      		// tmp1 = a0/r
      		tmp0 = inv_sqrt( drx*drx + dry*dry + drz*drz );
      		tmp1 = a0 * tmp0;

      		// tmp2 = phi(r) = epsilon * (a0/r)^n
      		tmp2 = epsilon * pow(tmp1, n);

      		// tmp3 = 0.5 * (phi(r) - phiCut)
      		tmp3 = half * ( tmp2 - phiCut );
      		tmp3 . store(en_+i);

      		// tmp2 = -n * phi(r) / r = dphi(r)
      		tmp3 = minus_n * tmp2 * tmp0;

      		// tmp1 = drho(r) = -m * (a0/r)^m / r
      		tmp1 = minus_m * pow(tmp1, m) * tmp0;

      		// tmp3 = de = drho(r) * (emb+emb[i]) + dphi(r)
      		// tmp2 = df = -de/r
      		tmp2 . load(emb_+i);
      		tmp1 = fmadd(tmp1, tmp2 + emb, tmp3);
      		tmp2 = minus_one * tmp1 * tmp0;

      		//
      		tmp0 = tmp2 * drx;
      		tmp1 = tmp2 * dry;
      		tmp3 = tmp2 * drz;

      		tmp0.store(fx_+i);
      		tmp1.store(fy_+i);
      		tmp3.store(fz_+i);

      	}

      }

    };


  }


}

#endif // __SUTTON_CHEN_HPP_INCLUDED
