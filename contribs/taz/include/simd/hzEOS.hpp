/// @file
/// @brief Vectorization of the HZ equation of state

#ifndef __HZ_EOS_HPP_INCLUDED
#define __HZ_EOS_INCLUDED

#include "libevi/simd.hpp"

#include "utils/stampUnits.hpp"


namespace simd {

  namespace kernels {



    /// @brief Compute thermodynamical quantities for the HZ equation of state
    /// @tparam T Numeric type
    template <class T> class HZ {

    private:

      vector_t<T> _rho0; ///< Reference density for the HZ equation of state
      vector_t<T> _g0; ///< Gruneisen parameter \f[ \Gamma_0 \f] of the HZ EoS
      vector_t<T> _rg0; ///< Gruneisen parameter \f[ \rho_0\Gamma_0 \f] of the HZ EoS
      vector_t<T> _c02; ///< Squared parameter \f[ c_0^2 \f] of the HZ EoS
      vector_t<T> _cv; ///< Heat capacity \f[ c_v \f] of the HZ EoS
      vector_t<T> _icv; ///< Inverse heat capacity of the HZ EoS
      vector_t<T> _icvm; ///< Parameter \f[ \frac{1}{m \times c_v} \f] of the HZ EoS
      vector_t<T> _s; ///< Parameter \f[ s \f] of the HZ EoS
      vector_t<T> _sgs; ///< Parameter \f[ s(\Gamma_0-s) \f] of the HZ EoS
      vector_t<T> _sg43ms; ///< Parameter \f[ (\frac43\Gamma_0-s) \f] of the HZ EoS
      vector_t<T> _delt0; ///< Temperature shift \f[ T_0 - T_{00} \f] of the HZ EoS
      vector_t<T> _third; ///< 1/3
      vector_t<T> _half; ///< 1/2
      vector_t<T> _one; ///< 1
      vector_t<T> _kb; ///< Boltzmann constant
      vector_t<T> _mass; ///< Mass
      vector_t<T> _imass; ///< Inverse mass


    public:

      /// @brief Set parameters
      inline void setParameters(const T& mass, const T& g0, const T& rho0, const T& c02,
				const T& cv, const T& s, const T& delt0) {

	_rho0        = vector_t<T>(rho0);
	_g0          = vector_t<T>(g0);
	_rg0         = _rho0*_g0;
	_c02         = vector_t<T>(c02);
	_cv          = vector_t<T>(cv);
	_icv         = vector_t<T>(1./cv);
	_icvm        = vector_t<T>(1./cv/mass);
	_s           = vector_t<T>(s);
	_sgs         = vector_t<T>(s*(g0-s));
	_sg43ms      = vector_t<T>(s*(4./3*g0-s));
	_delt0       = vector_t<T>(delt0);
	
	_third       = vector_t<T>(1./3);
	_half        = vector_t<T>(0.5);
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
	
	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ee.load(e_+i);

	  __SBP(dd,ee,ss,bb,pp,ic);
	  
	  ss.store(s_+i);
	  bb.store(b_+i);
	  pp.store(p_+i);
	  ic.store(ic_+i);
	  
	}

      }


      /// @brief Compute entropy, temperature and pressure from density and energy
      ///
      /// Vectorized computation
      /// @param [in] dd Density
      /// @param [in] ee Energy
      /// @param [out] ss Entropy
      /// @param [out] bb Inverse temperature
      /// @param [out] pp Pressure
      /// @param [out] ic Inverse heat capacity
      void __SBP(const vector_t<T>& dd, const vector_t<T>& ee, vector_t<T>& ss, vector_t<T>& bb, vector_t<T>& pp, vector_t<T>& ic) {

	vector_t<T> yy, tmp0, tmp1, er, pr;

	yy = _rho0/dd;

	epref(yy,er,pr);
	
	// deltaT0*exp(gamma0*(1-rho0/rho)
	tmp0 = _delt0*exp(_g0*(_one-yy));
	// (e/mass-eref(rho))/cv
	tmp1 = fmsub(ee,_imass,er)*_icv;
	
	ss = _cv*_mass*fmadd(_g0,yy,log(tmp0+tmp1));
	bb = _one/(_kb*(tmp0+tmp1));
	pp = fmadd(_rg0,tmp1*_cv,pr);
	ic = _icvm*_kb;
	
      }


      /// @brief Compute temperature and pressure and their derivatives from density and energy
      /// @param [in] dd Density
      /// @param [in] ee Energy
      /// @param [out] tt Temperature
      /// @param [out] drt Derivative of temperature with respect to density
      /// @param [out] det Derivative of temperature with respect to energy
      /// @param [out] pp Pressure
      /// @param [out] drp Derivative of pressure with respect to density
      /// @param [out] dep Derivative of pressure with respect to energy
      void __E2TPAndDerivatives(const vector_t<T>& dd, const vector_t<T>& ee, vector_t<T>& tt, vector_t<T>& drt, vector_t<T>& det, vector_t<T>& pp, vector_t<T>& drp, vector_t<T>& dep) {

	vector_t<T> yy, tmp0, tmp1, tmp2, er, pr, dpr;

	yy = _rho0/dd;

	epdref(yy,er,pr,dpr);
		
	// deltaT0*exp(gamma0*(1-rho0/rho)
	tmp0 = _delt0*exp(_g0*(_one-yy));
	// (e/mass-eref(rho))/cv
	tmp1 = fmsub(ee,_imass,er)*_icv;
	// 1/rho/rho
	tmp2 = _one/(dd*dd);
	
	tt  = tmp0+tmp1;
	drt = fmsub(_rg0,tmp0,pr*_icv)*tmp2;
	det = _icvm;
	pp  = fmadd(_rg0,tmp1*_cv,pr);
	drp = fnmadd(_rg0,pr*tmp2,dpr);
	dep = _rg0*_imass;
	
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
	
	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ss.load(s_+i);

	  __EBP(dd,ss,ee,bb,pp,ic);
	  
	  ee.store(e_+i);
	  bb.store(b_+i);
	  pp.store(p_+i);
	  ic.store(ic_+i);

	}

      }


      /// @brief Compute energy, temperature and pressure from density and entropy
      ///
      /// Vectorized computation
      /// @param [in] dd Density
      /// @param [in] ss Entropy
      /// @param [out] ee Energy
      /// @param [out] bb Inverse temperature
      /// @param [out] pp Pressure
      /// @param [out] ic Inverse heat capacity
      void __EBP(const vector_t<T>& dd, const vector_t<T>& ss, vector_t<T>& ee, vector_t<T>& bb, vector_t<T>& pp, vector_t<T>& ic) {

	vector_t<T> yy, tmp0, tmp1, er, pr;
	
	yy = _rho0/dd;

	epref(yy,er,pr);
	
	// deltaT0*exp(gamma0*(1-rho0/rho)
	tmp0 = _delt0*exp(_g0*(_one-yy));

	// e-eref(rho) = cv*exp(S/mass/cv-gamma0*rho0/rho) -cv*deltaT0*egx
	tmp1 = fmsub(ss,_icvm,_g0*yy);
	tmp1 = _cv*(exp(tmp1)-tmp0);

	ee = _mass*(er+tmp1);
	bb = _one/(_kb*fmadd(tmp1,_icv,tmp0));
	pp = fmadd(_rg0,tmp1,pr);
	ic = _icvm*_kb;

      }

      
      /// @brief Compute temperature and pressure and their derivatives from density and entropy
      /// @param [in] dd Density
      /// @param [in] ss Entropy
      /// @param [out] tt Temperature
      /// @param [out] drt Derivative of temperature with respect to density
      /// @param [out] dst Derivative of temperature with respect to entropy
      /// @param [out] pp Pressure
      /// @param [out] drp Derivative of pressure with respect to density
      /// @param [out] dsp Derivative of pressure with respect to entropy
      void __S2TPAndDerivatives(const vector_t<T>& dd, const vector_t<T>& ss, vector_t<T>& tt, vector_t<T>& drt, vector_t<T>& dst, vector_t<T>& pp, vector_t<T>& drp, vector_t<T>& dsp) {

	vector_t<T> yy, tmp0, tmp1, tmp2, tmp3, er, pr, dpr;

	yy = _rho0/dd;

	epdref(yy,er,pr,dpr);
	
	// deltaT0*exp(gamma0*(1-rho0/rho)
	tmp0 = _delt0*exp(_g0*(_one-yy));

	// e-eref(rho) = cv*exp(S/mass/cv-gamma0*rho0/rho) -cv*deltaT0*egx
	tmp1 = fmsub(ss,_icvm,_g0*yy);
	tmp2 = exp(tmp1)-tmp0;

	//1./rho/rho
	tmp3 = _one/dd/dd;
	  
	tt  = tmp1;
	drt = _rg0*tmp1*tmp3;
	dst = _icvm*tmp0;
	pp  = fmadd(_rg0,_cv*tmp2,pr);
	drp = fmadd(_rg0*_rg0,_cv*tmp2*tmp3,dpr);
	dsp = _rg0*_imass*tmp1;

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

	  __energy(dd,ss,ee);
	  
	  ee.store(e_+i);
	  
	}

      }


      /// @brief Compute only energy from density and entropy
      /// @param[in] dd Density
      /// @param[in] ss entropy
      /// @param [out] ee Energy
      void __energy (const vector_t<T>& dd, const vector_t<T>& ss, vector_t<T>& ee) {

	vector_t<T> yy, tmp0, tmp1, er;
	
	yy = _rho0/dd;

	eref(yy,er);
 
	// deltaT0*exp(gamma0*(1-rho0/rho)
	tmp0 = _delt0*exp(_g0*tmp0);
	
	// e-eref(rho) = cv*exp(S/mass/cv-gamma0*rho0/rho) -cv*deltaT0*egx
	tmp1 = fmsub(ss,_icvm,_g0*dd);
	tmp1 = _cv*(exp(tmp1)-tmp0);

	ee = _mass*(er+tmp1);

      }

      
    private:
      /// @brief Compute reference energy for HZ EOS
      /// @param [in] yy Reduced specific volume \f[ \frac{\rho0}{\rho} \f]
      /// @param [out] er Reference energy
      void eref(const vector_t<T>& yy, vector_t<T>& er) {

	vector_t<T> tmp0, tmp1, tmp2, tmp3;
		
	tmp0 = _one-yy;
	tmp1 = tmp0*tmp0;
	tmp2 = _one/fnmadd(_s,tmp0,_one);

	// eref(rho)
	tmp3 = _half*_c02*tmp1*tmp2;
	er   = fmadd(_s,tmp0*_third,_one);
	er   = fnmadd(_sgs,tmp1*_half*_third,er);
	
	er = blendv(tmp3,tmp3*er,yy<_one);
	
      }

      
      /// @brief Compute reference energy and pressure for HZ EOS
      /// @param [in] yy Reduced specific volume \f[ \frac{\rho0}{\rho} \f]
      /// @param [out] er Reference energy
      /// @param [out] pr Reference pressure
      void epref(const vector_t<T>& yy, vector_t<T>& er, vector_t<T>& pr) {

	vector_t<T> tmp0, tmp1, tmp2, tmp3;
	
	tmp0 = _one-yy;
	tmp1 = tmp0*tmp0;
	tmp2 = _one/fnmadd(_s,tmp0,_one);

	// eref(rho)
	tmp3 = _half*_c02*tmp1*tmp2;
	er   = fmadd(_s,tmp0*_third,_one);
	er   = fnmadd(_sgs,tmp1*_half*_third,er);
	
	er = blendv(tmp3,tmp3*er,yy<_one);
	
	// pref(rho)
	tmp3 = _rho0*_c02*tmp0*tmp2*tmp2;
	pr   = fnmadd(_s*_g0,tmp1*_third,_one);
	pr   = fmadd(_s*_sgs,tmp1*tmp0*_half*_half,pr);
	
	pr = tmp3*blendv(fnmadd(_s,tmp0*_half,_one),pr,yy<_one);

      }

	
      /// @brief Compute reference energy, pressure and pressure derivative for HZ EOS
      /// @param [in] yy Reduced specific volume \f[ \frac{\rho0}{\rho} \f]
      /// @param [out] er Reference energy
      /// @param [out] pr Reference pressure
      /// @param [out] dpr Derivative of the reference pressure
      void epdref(const vector_t<T>& yy, vector_t<T>& er, vector_t<T>& pr, vector_t<T>& dpr) {

	vector_t<T> tmp0, tmp1, tmp2, tmp3;
	
	tmp0 = _one-yy;
	tmp1 = tmp0*tmp0;
	tmp2 = _one/fnmadd(_s,tmp0,_one);

	// eref(rho)
	tmp3 = _half*_c02*tmp1*tmp2;
	er   = fmadd(_s,tmp0*_third,_one);
	er   = fnmadd(_sgs,tmp1*_half*_third,er);
	
	er = blendv(tmp3,tmp3*er,yy<_one);
	
	// pref(rho)
	tmp3 = _rho0*_c02*tmp0*tmp2*tmp2;
	pr   = fnmadd(_s*_g0,tmp1*_third,_one);
	pr   = fmadd(_s*_sgs,tmp1*tmp0*_half*_half,pr);
	
	pr = tmp3*blendv(fnmadd(_s,tmp0*_half,_one),pr,yy<_one);
	
	// dpref(rho)
	tmp3 = yy*yy*_c02*tmp2*tmp2*tmp2;
	tmp2 = _s*tmp1;
	dpr  = fmadd(_s,tmp0,_one);
	dpr  = fnmadd(tmp2,_g0,dpr);
	dpr  = fmadd(tmp2,_sg43ms*tmp0,dpr);
	dpr  = fnmadd(_half*tmp2,_sgs*tmp2,dpr);

	dpr = blendv(tmp3,tmp3*dpr,yy<_one);
	
      }
      
    };



  }

}



#endif // __HZ_EOS_HPP_INCLUDED
