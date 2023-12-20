/// @file
/// @brief Vectorization of chemical reactions

#ifndef __SIMD_CHEMISTRY_INCLUDED
#define __SIMD_CHEMISTRY_INCLUDED

#include "libevi/simd.hpp"


namespace simd {

  namespace kernels {

    /// @brief Update the progress variable
    /// @tparam T Type of the variables
    /// @param [in] time Time step
    /// @param [in] fchem Chemical evolution
    /// @param [in,out] progress Progress of the reaction
    /// @param [in,out] ei Internal energy
    /// @param [in] eexo Exothermicity
    /// @param [in] n Number of particles
    template <class T>
    inline void pushReaction(const T& time, const T* fchem, T* progress, T* ei, const T* eexo, const uint& n) {

      const vector_t<T> _minus(-1.);
      vector_t<T> dt(time);
      vector_t<T> idt(1./time);
      vector_t<T> evo, prog, eint, exo;
      
      for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	evo.load(fchem+i);
	prog.load(progress+i);
	eint.load(ei+i);
	exo.load(eexo+i);

	// Ensure that progress remains between 0 and 1
	evo = min(evo, (vector_t<T>::one()-prog)*idt);
	evo = max(evo, _minus*prog*idt);

	// Update progress and internal energy
	prog = fmadd(dt,evo,prog);
	eint = fmadd(evo,dt*exo,eint);

	prog.store(progress+i);
	eint.store(ei+i);
		
      }

    }


    /// @class SecondOrderReaction
    /// @brief Simd version of second order reactions kinetics
    /// @tparam T Type of the variables
    template <class T>
    class SecondOrderReaction {

    private:

      vector_t<T> _zab; ///< Arrhenius prefactor (direct reaction: A->B)
      vector_t<T> _zba; ///< Arrhenius prefactor (reverse reaction: A<-B)
      vector_t<T> _meab; ///< Activation energy (direct reaction: A->B)
      vector_t<T> _meba; ///< Activation energy (reverse reaction: A<-B)
            
    public:

      /// @brief Set parameters
      /// @param [in] zab Arrhenius prefactor (direct reaction: A->B)
      /// @param [in] zba Arrhenius prefactor (reverse reaction: A<-B)
      /// @param [in] eab Activation energy (direct reaction: A->B)
      /// @param [in] eba Activation energy (reverse reaction: A<-B)
      inline void setParameters(const T& zab, const T& zba, const T& eab, const T& eba) {

	_zab = vector_t<T>(zab);
	_zba = vector_t<T>(zba);
	_meab = vector_t<T>(-0.5*eab);
	_meba = vector_t<T>(-0.5*eba);

      }

      /// @brief Compute kinetics evolution
      /// @param [in] bij_ Mean inverse temperature \f[ \frac{1}{2} \frac{1}{k_{\rm B}} \left( \frac{1}{T_i} + \frac{1}{T_j} \right) \f]
      /// @param [in] wij_ Weight of the interaction
      /// @param [in] prog0_ Progress of the reaction for the first particle
      /// @param [in] prog1_ Progress of the reaction for the second particle
      /// @param [out] evo_ Evolution of the reaction
      /// @param [in] n number of neighbours
      void operator()(const T* bij_, const T* wij_, const T prog0_, const T* prog1_, T* evo_, const uint n) {

	vector_t<T> prog0(prog0_);
	vector_t<T> bij, wij, prog1, evo;
	vector_t<T> r0, r1;
	
	for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

	  bij.load(bij_+i);
	  wij.load(wij_+i);
	  prog1.load(prog1_+i);

	  // Compute reaction rate
	  r0 = _zab*exp(_meab*bij);
	  r1 = _zba*exp(_meba*bij);
	  
	  // Kinetics
	  evo = fmsub(r0,(vector_t<T>::one()-prog0)*(vector_t<T>::one()-prog1),r1*prog0*prog1);
	  evo = wij*evo;

	  evo.store(evo_+i);
	  	  
	}
	
      }
      
    };
      
  }

}


#endif //__SIMD_CHEMISTRY_INCLUDED
