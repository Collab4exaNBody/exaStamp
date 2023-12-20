
/// @file
/// @brief Vectorization of the Mie-Gruneisen equation of state

#ifndef __MG_EOS_HPP_INCLUDED
#define __MG_EOS_INCLUDED


#include "libevi/simd.hpp"

#include "utils/stampUnits.hpp"


namespace simd {

  namespace kernels {



    /// @brief Compute thermodynamical quantities for the Mie-Gruneisen equation of state
    /// @tparam T Numeric type
    template <class T> class MieGruneisen {

    private:

      vector_t<T> _rho0; ///< Parameter \f[ \rho_0 \f] of hte Mie-Gruneisen EoS
      vector_t<T> _q; ///< Parameter \f[ q \f] of the Mie-Gruneisen EoS
      vector_t<T> _miq; ///< Parameter \f[ -\frac{1}{q} \f] of the Mie-Gruneisen EoS
      vector_t<T> _ginf; ///< Parameter \f[ \Gamma_{\infty} \f] of the Mie-Gruneisen EoS
      vector_t<T> _g0i; ///< Parameter \f[ \Gamma_0 - \Gamma_{\infty} \f] of the Mie-Gruneisen EoS
      vector_t<T> _itheta0; ///< Parameter \f[ \frac{1}{\theta_0} \f] of the Mie-Gruneisen EoS
      vector_t<T> _mrhos; ///< Parameter \f[ -\rho_s \f] of the Mie-Gruneisen EoS
      vector_t<T> _irhos; ///< Parameter \f[ \frac1{\rho_s} \f] of the Mie-Gruneisen EoS
      vector_t<T> _npu; ///< Parameter \f[ N_s+1 \f] of the Mie-Gruneisen EoS
      vector_t<T> _inpu; ///< Parameter \f[ \frac1{N_s+1} \f] of the Mie-Gruneisen EoS
      vector_t<T> _kin; ///< Parameter \f[ \frac{K_s}{N_s+1} \f] of the Mie-Gruneisen EoS
      vector_t<T> _cvr; ///< Parameter \f[ Cv_r \f] of the Mie-Gruneisen EoS
      vector_t<T> _icvr; ///< Parameter \f[ \frac{1}{Cv_r} \f] of the Mie-Gruneisen EoS
      vector_t<T> _icv; ///< Parameter \f[ \frac{1}{m \times Cv_r} \f] of the Mie-Gruneisen EoS
      vector_t<T> _ur; ///< Parameter \f[ u_r \f] of the Mie-Gruneisen EoS
      vector_t<T> _iur; ///< Parameter \f[ \frac{1}{u_r} \f] of the Mie-Gruneisen EoS
      vector_t<T> _one; ///< 1
      vector_t<T> _kb; ///< Boltzmann constant
      vector_t<T> _mass; ///< Mass
      vector_t<T> _imass; ///< Inverse mass


    public:

      /// @brief Set parameters
      inline void setParameters(const T& mass, const T& g0, const T& ginf, const T& t0, const T& q,
				const T& d0, const T& ks, const T& npu, const T& ds, const T& ur, const T& cvr) {

	_rho0        = vector_t<T>(d0);
	_q           = vector_t<T>(q);
	_miq         = vector_t<T>(-1./q);
	_ginf        = vector_t<T>(ginf);
	_g0i         = vector_t<T>(g0-ginf);
	_itheta0     = vector_t<T>(1./t0);
	
	_mrhos       = vector_t<T>(-ds);
	_irhos       = vector_t<T>(1./ds);
	_npu         = vector_t<T>(npu);
	_inpu        = vector_t<T>(1./npu);
	_kin         = vector_t<T>(ks/npu);

	_cvr         = vector_t<T>(cvr);
	_icvr        = vector_t<T>(1./cvr);
	_icv         = vector_t<T>(1./cvr/mass);
	_ur          = vector_t<T>(ur);
	_iur         = vector_t<T>(1./ur);
	
	_one         = vector_t<T>(1.);
	_kb          = vector_t<T>(Stamp_Constant::boltzmann);
	_mass        = vector_t<T>(mass);
	_imass       = vector_t<T>(1./mass);
	
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

	vector_t<T> ss, bb, pp, dd, ee, ic;
	vector_t<T> tmp0, tmp1, tmp2, gamma, itheta, pote, potp, redu;
	
	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ee.load(e_+i);

	  tmp0 = _one/dd;

	  // log(rho0/rho)
	  tmp1 = log(_rho0*tmp0);

	  // gamma(rho)
	  // gammaInf + (gamma0-gammaInf)*exp(q*log(rho0/rho))
	  // theta(rho)
	  // theta0*exp(-gammaInf*log(rho0/rho) + (gamma0-gammaInf)/q*(1-exp(q*log(rho0/rho))))
	  tmp2 = exp(_q*tmp1);
	  gamma = fmadd(_g0i,tmp2,_ginf);
	  tmp2 = _one-tmp2;
	  tmp1 = fmadd(_ginf,tmp1,_g0i*_miq*tmp2);
	  itheta = _itheta0*exp(tmp1);

 
	  // pote(rho)
	  // ks/rhos*(exp(npu*(1-rhos/rho))/npu - (1-rhos/rho))/npu
	  // potp(rho)
	  // ks/npu*(exp(npu*(1-rhos/rho))-1);
	  tmp0 = fmadd(_mrhos,tmp0,_one);
	  tmp1 = exp(_npu*tmp0);
	  
	  pote = _kin*_irhos*(tmp1*_inpu-tmp0);
	  potp = _kin*(tmp1-_one);

	  // reduced energy
	  // (e/m - pote(rho))/cvr/theta(rho) + ur
	  tmp0 = ee*_imass;
	  tmp0 = tmp0 - pote;
	  redu = fmadd(tmp0 , _icvr*itheta, _ur);

	  // We want to compute (Mie-Gruneisen)
	  // s = cvr*log(redu/ur)
	  // t = redu*theta
	  // p = potp + rho*gamma*(e/m - pote)
	  ss = _cvr*_mass*log(redu*_iur);
	  bb = itheta/(redu*_kb);
	  pp = potp + dd*gamma*tmp0;
	  ic = _icv*_kb;
	  
	  // store back the results
	  ss .store(s_+i);
	  bb .store(b_+i);
	  pp .store(p_+i);
	  ic .store(ic_+i);

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

	vector_t<T> ee, bb, pp, dd, ss, ic;
	vector_t<T> tmp0, tmp1, tmp2, gamma, itheta, pote, potp;
	
	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ss.load(s_+i);

	  tmp0 = _one/dd;
	  
	  // log(rho0/rho)
	  tmp1 = log(_rho0*tmp0);

	  // gamma(rho)
	  // gammaInf + (gamma0-gammaInf)*exp(q*log(rho0/rho))
	  // theta(rho)
	  // theta0*exp(-gammaInf*log(rho0/rho) + (gamma0-gammaInf)/q*(1-exp(q*log(rho0/rho))))
	  tmp2 = exp(_q*tmp1);
	  gamma = fmadd(_g0i,tmp2,_ginf);
	  tmp2 = _one-tmp2;
	  tmp1 = fmadd(_ginf,tmp1,_g0i*_miq*tmp2);
	  itheta = _itheta0*exp(tmp1);

	  // pote(rho)
	  // ks/rhos*(exp(npu*(1-rhos/rho))/npu - (1-rhos/rho))/npu
	  // potp(rho)
	  // ks/npu*(exp(npu*(1-rhos/rho))-1);
	  tmp0 = fmadd(_mrhos,tmp0,_one);
	  tmp1 = exp(_npu*tmp0);
	  
	  pote = _kin*_irhos*(tmp1*_inpu-tmp0);
	  potp = _kin*(tmp1-_one);

	  // We want to compute (Mie-Gruneisen)
	  // e/m = cvr*ur*theta*(exp(S/cvr)-1) + pote
	  // t = ur*theta*exp(S/cvr)
	  // p = potp + rho*gamma*cvr*ur*theta*(exp(S/cvr)-1)
	  tmp0 = _cvr*_ur/itheta;
	  tmp1 = exp(ss*_icv);
	  tmp0 = tmp0*(tmp1-_one);

	  ee = (tmp0 + pote)*_mass;
	  bb = _iur*itheta/(tmp1*_kb);
	  pp = dd*gamma*tmp0 + potp;
	  ic = _icv*_kb;
	  
	  // store back the results
	  ee .store(e_+i);
	  bb .store(b_+i);
	  pp .store(p_+i);
	  ic .store(ic_+i);

	}

      }



      /// @brief Compute only energy from density and entropy
      /// @param [out] e_ Energy
      /// @param[in] d_ Density
      /// @param[in] s_ entropy
      /// @param[in] n Number of particles
      void energy (T* e_, const T* d_, const T* s_, const uint n) {

	vector_t<T> ee, dd, ss;
	vector_t<T> tmp0, tmp1, tmp2, gamma, itheta, pote, potp;
	
	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ss.load(s_+i);

	  tmp0 = _one/dd;
	  
	  // log(rho0/rho)
	  tmp1 = log(_rho0*tmp0);

	  // theta(rho)
	  // theta0*exp(-gammaInf*log(rho0/rho) + (gamma0-gammaInf)/q*(1-exp(q*log(rho0/rho))))
	  tmp2 = _one-exp(_q*tmp1);
	  tmp1 = fmadd(_ginf,tmp1,_g0i*_miq*tmp2);
	  itheta = _itheta0*exp(tmp1);
	  
	  // pote(rho)
	  // ks/rhos*(exp(npu*(1-rhos/rho))/npu - (1-rhos/rho))/npu
	  // potp(rho)
	  // ks/npu*(exp(npu*(1-rhos/rho))-1);
	  tmp0 = fmadd(_mrhos,tmp0,_one);
	  tmp1 = exp(_npu*tmp0);
	  
	  pote = _kin*_irhos*(tmp1*_inpu-tmp0);

	  // We want to compute (Mie-Gruneisen)
	  // e/m = cvr*ur*theta*(exp(S/cvr)-1) + pote
	  // t = ur*theta*exp(S/cvr)
	  // p = potp + rho*gamma*cvr*ur*theta*(exp(S/cvr)-1)
	  tmp0 = _cvr*_ur/itheta;
	  tmp1 = tmp0*(exp(ss*_icv)-_one);

	  ee = (tmp1 + pote)*_mass;

	  // store back the results
	  ee .store(e_+i);

	}

      }

    };



  }

}



#endif // __MG_EOS_HPP_INCLUDED
