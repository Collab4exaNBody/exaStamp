/// @file
/// @brief Vectorization of the JWL equation of state

#ifndef __JWL_EOS_HPP_INCLUDED
#define __JWL_EOS_INCLUDED

#include "libevi/simd.hpp"

#include "utils/stampUnits.hpp"


namespace simd {

  namespace kernels {



    /// @brief Compute thermodynamical quantities for the JWL equation of state
    /// @tparam T Numeric type
    template <class T> class JWL {

    private:

      vector_t<T> _rho0; ///< Reference density for the JWL equation of state
      vector_t<T> _g0; ///< Gruneisen parameter \f[ \Gamma_0 \f] of the JWL EoS

      vector_t<T> _cv; ///< Heat capacity \f[ c_v \f] of the JWL EoS
      vector_t<T> _icv; ///< Inverse heat capacity of the JWL EoS
      vector_t<T> _icvm; ///< Parameter \f[ \frac{1}{m \times c_v} \f] of the JWL EoS

      vector_t<T> _a; ///< Parameter \f[ a \f] of the JWL EoS
      vector_t<T> _b; ///< Parameter \f[ b \f] of the JWL EoS
      vector_t<T> _r1; ///< Parameter \f[ R_1 \f] of the JWL EoS
      vector_t<T> _r2; ///< Parameter \f[ R_2 \f] of the JWL EoS
      
      vector_t<T> _cek; ///< Parameter \f[ C_{ek} \f] of the JWL EoS
      vector_t<T> _k; ///< Temperature shift \f[ K \f] of the JWL EoS

      vector_t<T> _minus; ///< -1
      vector_t<T> _one; ///< 1
      vector_t<T> _kb; ///< Boltzmann constant
      vector_t<T> _mass; ///< Mass
      vector_t<T> _imass; ///< Inverse mass


    public:

      /// @brief Set parameters
      inline void setParameters(const T& mass, const T& g0, const T& rho0, const T& cv,
				const T& a, const T& b, const T& R1, const T& R2,
				const T& cek, const T& k) {

	_rho0        = vector_t<T>(rho0);
	_g0          = vector_t<T>(g0);

	_cv          = vector_t<T>(cv);
	_icv         = vector_t<T>(1./cv);
	_icvm        = vector_t<T>(1./cv/mass);

	_a           = vector_t<T>(a);
	_b           = vector_t<T>(b);
	_r1          = vector_t<T>(R1);
	_r2          = vector_t<T>(R2);

	_cek         = vector_t<T>(cek);
	_k           = vector_t<T>(k);

	_minus       = vector_t<T>(-1);
	_one         = vector_t<T>(1);
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

	vector_t<T> tmp0, ek, pk;

	tmp0 = _rho0/dd;
	
	epref(tmp0,ek,pk);
	
	// (e/mass-ek(rho))/cv
	tmp0 = fmsub(ee,_imass,ek)*_icv;
	
	ss = _cv*_mass*fnmadd(_g0,log(dd),log(tmp0));
	bb = _one/(tmp0*_kb);
	pp = fmadd(dd*_g0,tmp0*_cv,pk);
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

	vector_t<T> tmp0, ek, pk, dpk;

	tmp0 = _rho0/dd;
	
	epdref(tmp0,ek,pk,dpk);

	tmp0 = fmsub(ee,_imass,ek);
	
	tt  = tmp0*_icv;
	drt = _minus*pk*_icv/(dd*dd);
	det = _icvm;
	pp  = fmadd(dd*_g0,tmp0,pk);
	drp = fmadd(_g0,tmp0,dpk);
	drp = fnmadd(_g0/dd,pk,drp);
	dep = dd*_g0*_imass;

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

	vector_t<T> tmp0, ek, pk;
	
	tmp0 = _rho0/dd;
	
	epref(tmp0,ek,pk);
	
	// e/mass-ek(rho)
	tmp0 = fmadd(ss,_imass*_icv,_g0*log(dd));
	tmp0 = _cv*exp(tmp0);
	
	ee = _mass*(ek+tmp0);
	bb = _one/(tmp0*_icv*_kb);
	pp = fmadd(dd*_g0,tmp0,pk);
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

	vector_t<T> tmp0, ek, pk,dpk;
	
	tmp0 = _rho0/dd;
	
	epdref(tmp0,ek,pk,dpk);
	
	// e/mass-ek(rho)
	tmp0 = fmadd(ss,_imass*_icv,_g0*log(dd));
	tmp0 = _cv*exp(tmp0);
	
	tt = tmp0*_icv;
	drt = _g0/dd*tt;
	dst = _icvm*tt;

	tmp0 = _g0*tmp0;
	
	pp = fmadd(dd*_g0,tmp0,pk);
	drp = fmadd(_g0,tmp0,dpk+tmp0);
	dsp = dd*_g0*_imass*tt;
	
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

	vector_t<T> tmp0, ek;

	tmp0 = _rho0/dd;
	
	eref(tmp0,ek);
	
	// e/mass-ek(rho)
	tmp0 = fmadd(ss,_imass,_g0*log(dd));
	tmp0 = _cv*exp(tmp0);
	
	ee = _mass*(ek+tmp0);

      }
      

    private:
      /// @brief Compute reference energy for JWL EOS
      /// @param [in] yy Reduced specific volume \f[ \frac{\rho0}{\rho} \f]
      /// @param [out] ek Reference energy
      void eref(const vector_t<T>& yy, vector_t<T>& ek) {

	vector_t<T> tmp0, tmp1, tmp2;
	
	// ek(rho)
	tmp0 = exp(_minus*_r1*yy);
	tmp1 = exp(_minus*_r2*yy);
	tmp2 = pow(yy,fmsub(_minus,_g0,_one));

	ek = fmadd(_a/_r1,tmp0,_b/_r2*tmp1);
	ek = fmadd(_k/_g0,tmp2*yy,ek)/_rho0 + _cek;

      }

      
      /// @brief Compute reference energy, pressure for JWL EOS
      /// @param [in] yy Reduced specific volume \f[ \frac{\rho0}{\rho} \f]
      /// @param [out] ek Reference energy
      /// @param [out] pk Reference pressure
      void epref(const vector_t<T>& yy, vector_t<T>& ek, vector_t<T>& pk) {

	vector_t<T> tmp0, tmp1, tmp2;
	
	// ek(rho)
	tmp0 = exp(_minus*_r1*yy);
	tmp1 = exp(_minus*_r2*yy);
	tmp2 = pow(yy,fmsub(_minus,_g0,_one));

	ek = fmadd(_a/_r1,tmp0,_b/_r2*tmp1);
	ek = fmadd(_k/_g0,tmp2*yy,ek)/_rho0 + _cek;

	// pk(rho)
	pk = fmadd(_a,tmp0,_b*tmp1);
	pk = fmadd(_k,tmp2,pk);

      }
      
      
      /// @brief Compute reference energy, pressure and pressure derivative for JWL EOS
      /// @param [in] yy Reduced specific volume \f[ \frac{\rho0}{\rho} \f]
      /// @param [out] ek Reference energy
      /// @param [out] pk Reference pressure
      /// @param [out] dpk Derivative of the reference pressure
      void epdref(const vector_t<T>& yy, vector_t<T>& ek, vector_t<T>& pk, vector_t<T>& dpk) {

	vector_t<T> tmp0, tmp1, tmp2;
		
	// ek(rho)
	tmp0 = exp(_minus*_r1*yy);
	tmp1 = exp(_minus*_r2*yy);
	tmp2 = pow(yy,fmsub(_minus,_g0,_one));

	ek = fmadd(_a/_r1,tmp0,_b/_r2*tmp1);
	ek = fmadd(_k/_g0,tmp2*yy,ek)/_rho0 + _cek;

	// pk(rho)
	pk = fmadd(_a,tmp0,_b*tmp1);
	pk = fmadd(_k,tmp2,pk);

	// dpk(rho)
	dpk = fmadd(_a*_r1,tmp0,_b*_r2*tmp1)*yy*yy;
	dpk = fmadd(_k*(_g0+_one),tmp2*yy,dpk)/_rho0;
	
      }
      
    };



  }

}



#endif // __JWL_EOS_HPP_INCLUDED
