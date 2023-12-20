/// @file
/// @brief Vectorization of sdpd kernels

#ifndef __SDPD_KERNEL_HPP_INCLUDED
#define __SDPD_KERNEL_HPP_INCLUDED


#include "libevi/simd.hpp"
#include "utils/kernelFunction/lucy.hpp"


namespace simd {

  namespace kernels {



    /// @brief Compute kernel relative quantities
    /// @tparam T Numeric type
    /// @tparam KernelFunctionA Kernel type for particle type A
    /// @tparam KernelFunctionB Kernel type for particle type B
    template <class T, class KernelFunctionA, class KernelFunctionB> class SdpdKernel {
    public:
      /// @brief Set the kernel parameters
      inline void setParameters(KernelFunctionA* A, KernelFunctionB* B, T ma_, T mb_, T param_=0) {}
      /// @brief Evaluate the kernel
      void operator () (T* rhoA_, T* rhoB_, const T* drx_, const T* dry_, const T* drz_, const uint n) {}
      /// @brief Evaluate the kernel derivative
      void operator () (T* fx_, T* fy_, T* fz_, const T* drx_, const T* dry_, const T* drz_, const uint n) {}
      /// @brief Compute the cut-off function
      void computeChi (T* chi_,  const T* drx_, const T* dry_, const T* drz_, const T* rhoj_, const uint n) {}
      };
      



    /// @brief Compute kernel related quantities for the Lucy kernel
    /// @tparam T Numeric type
    template <class T> class SdpdKernel<T,Lucy,Lucy> {

    private:

      vector_t<T> ma; ///< Mass of type A particles
      vector_t<T> mb; ///< Mass of type B particles
      vector_t<T> iha; ///< Inverse smoothing length for type A particle
      vector_t<T> ihb; ///< Inverse smoothing length for type B particle
      vector_t<T> awa; ///< Kernel normalisation factor for type A particles
      vector_t<T> awb; ///< Kernel normalisation factor for type B particles
      vector_t<T> afa; ///< Kernel derivative normalisation factor for type B particles
      vector_t<T> afb; ///< Kernel derivative normalisation factor for type B particles
      vector_t<T> three; ///< 3
      vector_t<T> pir2; ///< Pressure over squared density for the considered particle


    public:

      SdpdKernel<T,Lucy,Lucy>() : three(3.) {}

      /// @brief Set parameters
      inline void setParameters(Lucy* A, Lucy* B, T ma_, T mb_, T pir2_=0) {
	
	ma = vector_t<T>(ma_);
	mb = vector_t<T>(mb_);
	iha = vector_t<T>(A->getInvSmoothingLength());
	ihb = vector_t<T>(B->getInvSmoothingLength());
	awa = vector_t<T>(A->getAlphaW());
	awb = vector_t<T>(B->getAlphaW());
	afa = vector_t<T>(A->getAlphaF());
	afb = vector_t<T>(B->getAlphaF());
	pir2 = vector_t<T>(pir2_);

      }



      /// @brief compute mass*kernel
      void operator () (T* rhoA_, T* rhoB_, const T* drx_, const T* dry_, const T* drz_, const uint n) {

	vector_t<T> drx, dry, drz, tmp0, tmp1, tmp2;

 	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  drx.load(drx_+i);
	  dry.load(dry_+i);
	  drz.load(drz_+i);

	  // We want to compute mb*awa*wab and ma*awb*wba
	  // where wij = (1+3*r/hi)*(1-r/hi)^3
	  tmp0 = sqrt(drx*drx + dry*dry + drz*drz);

	  tmp1 = tmp0 * iha;
	  tmp2 = vector_t<T>::one() - tmp1;
	  tmp2 = mb*awa*fmadd(three, tmp1, vector_t<T>::one())*tmp2*tmp2*tmp2;

	  tmp0 = tmp0 * ihb;
	  tmp1 = vector_t<T>::one() - tmp0;
	  tmp1 = ma*awb*fmadd(three, tmp0, vector_t<T>::one())*tmp1*tmp1*tmp1;
	  
	  tmp2.store(rhoA_+i);
	  tmp1.store(rhoB_+i);
	}

      }

      /// @brief Compute kernel derivative
      void operator () (T* fx_, T* fy_, T* fz_, const T* drx_, const T* dry_, const T* drz_, const T* prj_, const T* rhoj_, const uint n) {

	vector_t<T> drx, dry, drz, tmp0, tmp1, tmp2, tmp3;

 	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  drx.load(drx_+i);
	  dry.load(dry_+i);
	  drz.load(drz_+i);

	  // We want to compute ma*mb*(pri/rhoi^2*fij + prj/rho_j^2*fji)*rij
	  // where fij = afi*(1-r/hi)^2
	  tmp0 = sqrt(drx*drx + dry*dry + drz*drz);

	  tmp1 = vector_t<T>::one() - tmp0*iha;
	  tmp1 = pir2*afa*tmp1*tmp1;
	    
	  tmp2.load(prj_+i);
	  tmp3.load(rhoj_+i);
	  tmp2 = tmp2/(tmp3*tmp3);

	  tmp0 = vector_t<T>::one() - tmp0*ihb;
	  tmp1 = fmadd(tmp2,afb*tmp0*tmp0,tmp1);
	  tmp1 = ma*mb*tmp1;

	  drx = tmp1*drx;
	  dry = tmp1*dry;
	  drz = tmp1*drz;

	  drx.store(fx_+i);
	  dry.store(fy_+i);
	  drz.store(fz_+i);

	}
	
      }

      
      /// @brief Compute force with wall
      template <bool iIsWall>
      void wall (T* fx_, T* fy_, T* fz_, const T* drx_, const T* dry_, const T* drz_, const T* prj_, const T* rhoj_, const uint n) {

	vector_t<T> drx, dry, drz, tmp0, tmp1, tmp2, tmp3;

 	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  drx.load(drx_+i);
	  dry.load(dry_+i);
	  drz.load(drz_+i);

	  // We want to compute ma*mb*(pri/rhoi^2*fij + prj/rho_j^2*fji)*rij
	  // where fij = afi*(1-r/hi)^2
	  tmp0 = sqrt(drx*drx + dry*dry + drz*drz);

	  tmp2.load(prj_+i);
	  tmp3.load(rhoj_+i);
	  tmp2 = tmp2/(tmp3*tmp3);

	  if(iIsWall)
	    pir2 = tmp2;
	  else
	    tmp2 = pir2;
	  
	  tmp1 = vector_t<T>::one() - tmp0*iha; 
	  tmp1 = pir2*afa*tmp1*tmp1;

	  tmp0 = vector_t<T>::one() - tmp0*ihb; 
	  tmp1 = fmadd(tmp2,afb*tmp0*tmp0,tmp1);
	  tmp1 = ma*mb*tmp1;

	  drx = tmp1*drx;
	  dry = tmp1*dry;
	  drz = tmp1*drz;

	  drx.store(fx_+i);
	  dry.store(fy_+i);
	  drz.store(fz_+i);

	}
	
      }

      
      /// @brief Compute cut-off function chi for fluctuation-dissipation
      void computeChi (T* chi_,  const T* drx_, const T* dry_, const T* drz_, const T* rhoj_, const uint n) {

	vector_t<T> drx, dry, drz, tmp0, tmp1, tmp2, tmp3;

 	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  drx.load(drx_+i);
	  dry.load(dry_+i);
	  drz.load(drz_+i);

	  // We want to compute ma*mb*(fij+fji)/2/rhoi^2
	  // where fij = afi*(1-r/hi)^2
	  tmp0 = sqrt(drx*drx + dry*dry + drz*drz);

	  tmp1 = vector_t<T>::one() - tmp0*iha;
	  tmp1 = afa*tmp1*tmp1;
	    
	  // tmp0 = vector_t<T>::one() - tmp0*ihb;
	  // tmp1 = afb*tmp0*tmp0+tmp1;//fmadd(afb,tmp0*tmp0,tmp1);

	  tmp2.load(rhoj_+i);
	  tmp2 = vector_t<T>(0.5)*ma*mb*tmp1*pir2/tmp2;

	  tmp2.store(chi_+i);

	}
	
      }

    };



  }

}



#endif // __SDPD_KERNELS_HPP_INCLUDED
