/// @file
/// @brief Vectorization of the reactive equation of state

#ifndef __SIMD_REACTIVE_EOS_HPP_INCLUDED
#define __SIMD_REACTIVE_EOS_HPP_INCLUDED

#include "libevi/simd.hpp"

#include "utils/stampUnits.hpp"


namespace simd {

  namespace kernels {

    /// @brief Compute thermodynamical quantities for the reactive equation of state
    /// @tparam T Numeric type
    /// @tparam EOS0 First equation of state
    /// @tparam EOS1 Second equation of state
    template <class T, class EOS0, class EOS1> class Reactive {

    private:
      int imax; ///< Maximum number of iteration
      vector_t<T> tol; ///< Tolerance for the convergence (not used for now)
      
      typename EOS0::simd_opt_t* eos0; ///< First equation of state (simd version)
      typename EOS1::simd_opt_t* eos1; ///< Second equation of state (simd version)

      typedef void (EOS0::simd_opt_t::*FThermo0)(const vector_t<T>&, const vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&);
      typedef void (EOS1::simd_opt_t::*FThermo1)(const vector_t<T>&, const vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&, vector_t<T>&);

      /// @brief Find reactant state corresponding to mixed state
      /// @param [in] progress Progress variable
      /// @param [in] rho Density
      /// @param [in] var Internal variable (energy or entropy)
      /// @param [out] rho0 Density for the reactant
      /// @param [out] var0 Internel variable (energy or entropy) for the reactant
      /// @param [in] fct0 Function to minimize for the first equation of state
      /// @param [in] fct1 Function to minimize for the second equation of state
      void newton(const vector_t<T>& progress, const vector_t<T>& rho, const vector_t<T>& var, vector_t<T>& rho0, vector_t<T>& var0, FThermo0 fct0, FThermo1 fct1) {
	
	// Initialize with progress = 0
	rho0 = rho;
	var0 = var;

	vector_t<T> ratio = (vector_t<T>::one()-progress)/progress;
	vector_t<T> rl = rho/progress;
	vector_t<T> vl = var/progress;
	
	vector_t<T> t0, p0, drt0, dvt0, drp0, dvp0;
	vector_t<T> t1, p1, drt1, dvt1, drp1, dvp1;

	(eos0->*fct0)(rho0,var0,t0,drt0,dvt0,p0,drp0,dvp0);
	(eos1->*fct1)(fnmadd(ratio,rho0,rl),fnmadd(ratio,var0,vl),t1,drt1,dvt1,p1,drp1,dvp1);
	
	t0 = (t0-t1)*(t0-t1);
	p0 = (p0-p1)*(p0-p1);

	int iter = 0;
	//	bool crit = bool(movemask(blendv(t0+p0,tol,progress==vector_t<T>::zero()) > tol)) && iter<imax;
	
	//	while (crit) {
	while (iter<imax) {
	
	  // Get temperature, pressure and their derivatives
	  (eos0->*fct0)(rho0,var0,t0,drt0,dvt0,p0,drp0,dvp0);
	  (eos1->*fct1)(fnmadd(ratio,rho0,rl),fnmadd(ratio,var0,vl),t1,drt1,dvt1,p1,drp1,dvp1);

	  // Error function
	  t0 = t0-t1;
	  p0 = p0-p1;

	  // Jacobian matrix
	  drt0 = fmadd(ratio,drt1,drt0);
	  dvt0 = fmadd(ratio,dvt1,dvt0);
	  drp0 = fmadd(ratio,drp1,drp0);
	  dvp0 = fmadd(ratio,dvp1,dvp0);

	  t1 = vector_t<T>::one()/fmsub(dvt0,drp0,drt0*dvp0);

	  // Update the variables
	  rho0 = fmadd(fmsub(dvp0,t0,dvt0*p0),t1,rho0);
	  var0 = fmadd(fnmadd(drp0,t0,drt0*p0),t1,var0);

	  iter++;
	  //	  crit = bool(movemask(blendv(t0+p0,tol,progress==vector_t<T>::zero()) > tol)) && iter<imax;
	  
	}

	rho0 = blendv(rho0,rho,progress==vector_t<T>::zero());
	var0 = blendv(var0,var,progress==vector_t<T>::zero());
	
      }
      
    public :

      /// @brief Constructor
      Reactive<T,EOS0,EOS1>() {
	eos0 = new typename EOS0::simd_opt_t();
	eos1 = new typename EOS1::simd_opt_t();
      }

      /// @brief Destructor
      ~Reactive<T,EOS0,EOS1>() {
	delete eos0;
	delete eos1;
      }
      

      /// @brief Set the parameters of the related equation of states
      /// @param [in] eq0 First equation of state
      /// @param [in] eq1 Second equation of state
      void setParameters(EOS0* eq0, EOS1* eq1) {
	eq0->setParameters(*eos0);
	eq1->setParameters(*eos1);

	imax = 25;
	tol = vector_t<T>(1e-15);
      }
	
      
      /// @brief Compute entropy, temperature and pressure from density and energy
      /// @param [out] s_ Entropy
      /// @param [out] b_ Inverse temperature
      /// @param [out] p_ Pressure
      /// @param [in,out] ic_ Progress variable (in) / Inverse heat capacity (out)
      /// @param [in] d_ Density
      /// @param [in] e_ Energy
      /// @param [in] n Number of particles
      void SBP (T* s_, T* b_, T* p_, T* ic_, const T* d_, const T* e_, const uint n) {

	vector_t<T> dd, ee, dd0, ee0;
	vector_t<T> ss0, bb0, pp0, ic0, ss1, ic1;
	vector_t<T> progress, ratio;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ee.load(e_+i);

	  progress.load(ic_+i);
	  // T *tmp=(T*)(&progress.data());
	  // for (uint k = i+vector_t<T>::chunk_size; i < n; i++)
	  //   tmp[k] = 0.;

	  // Invert to find the mixed state
	  newton(progress,dd,ee,dd0,ee0,&EOS0::simd_opt_t::__E2TPAndDerivatives,&EOS1::simd_opt_t::__E2TPAndDerivatives);
	  
	  ratio = (vector_t<T>::one()-progress)/progress;

	  // Compute the component state
	  // Temperature and pressure should be the same !
	  eos1->__SBP(dd/progress-ratio*dd0,ee/progress-ratio*ee0,ss1,bb0,pp0,ic1);
	  eos0->__SBP(dd0,ee0,ss0,bb0,pp0,ic0);
  
	  ss0 = blendv( fmadd(vector_t<T>::one()-progress,ss0,progress*ss1), ss0, progress==vector_t<T>::zero() );
	  ic0 = blendv( vector_t<T>::one()/((vector_t<T>::one()-progress)/ic0+progress/ic1), ic0, progress==vector_t<T>::zero() );

	  ss0.store(s_+i);
	  bb0.store(b_+i);
	  pp0.store(p_+i);
	  ic0.store(ic_+i);
	  
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

	vector_t<T> dd, ss, dd0, ss0;
	vector_t<T> ee0, bb0, pp0, ic0, ee1, ic1;
	vector_t<T> progress, ratio;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ss.load(s_+i);

	  progress.load(ic_+i);
	  // T *tmp=(T*)(&progress.data());
	  // for (uint k = i+vector_t<T>::chunk_size; i < n; i++)
	  //   tmp[k] = 0.;

	  // Invert to find the mixed state
	  newton(progress,dd,ss,dd0,ss0,&EOS0::simd_opt_t::__S2TPAndDerivatives,&EOS1::simd_opt_t::__S2TPAndDerivatives);
	  
	  ratio = (vector_t<T>::one()-progress)/progress;

	  // Compute the component state
	  // Temperature and pressure should be the same !
	  eos1->__EBP(dd/progress-ratio*dd0,ss/progress-ratio*ss0,ee1,bb0,pp0,ic1);
	  eos0->__EBP(dd0,ss0,ee0,bb0,pp0,ic0);
	  
	  ee0 = blendv( fmadd(vector_t<T>::one()-progress,ee0,progress*ee1), ee0, progress==vector_t<T>::zero() );
	  ic0 = blendv( vector_t<T>::one()/((vector_t<T>::one()-progress)/ic0+progress/ic1), ic0, progress==vector_t<T>::zero() );

	  ee0.store(e_+i);
	  bb0.store(b_+i);
	  pp0.store(p_+i);
	  ic0.store(ic_+i);
	  
	}
	
      }

      /// @brief Compute only energy from density and entropy
      /// @param [in,out] e_ Progress (in) / Energy (out)
      /// @param [in] d_ Density
      /// @param [in] s_ entropy
      /// @param [in] n Number of particles
      void energy (T* e_, const T* d_, const T* s_, const uint n) {

	vector_t<T> dd, ss, dd0, ss0;
	vector_t<T> ee0, ee1;
	vector_t<T> progress, ratio;

	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  // load data
	  dd.load(d_+i);
	  ss.load(s_+i);

	  progress.load(e_+i);
	  for (uint k = i+vector_t<T>::chunk_size; i < n; i++)
	    progress[i] = 0.;

	  // Invert to find the mixed state
	  newton(progress,dd,ss,dd0,ss0,&EOS0::simd_opt_t::__S2TPAndDerivatives,&EOS1::simd_opt_t::__S2TPAndDerivatives);
	  
	  ratio = (vector_t<T>::one()-progress)/progress;

	  // Compute the component state
	  // Temperature and pressure should be the same !
	  eos1->__energy(dd-ratio*dd0,ss-ratio*ss0,ee1);
	  eos0->__energy(dd0,ss0,ee0);
	  
	  ee0 = blendv( fmadd(vector_t<T>::one()-progress,ee0,progress*ee1), ee0, progress==vector_t<T>::zero() );

	  ee0.store(e_+i);
		  
	}

      }
      
    };
    
  }
  
}

#endif // __SIMD_REACTIVE_EOS_HPP_INCLUDED
