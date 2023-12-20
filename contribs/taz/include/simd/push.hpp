/// @file
/// @brief Vectorization of "push" functions

#ifndef __PUSH_HPP_INCLUDED
#define __PUSH_HPP_INCLUDED


#include "libevi/simd.hpp"
#include "simd/random.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
  namespace kernels {


    /// @brief Vectorized reset of an array
    /// @tparam T Type of the variables
    /// @param [in,out] t_ Array to reset
    /// @param [in] n Size of the array
    template <class T>
    inline void reset (T* t_, const uint n) {

      for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {
      	vector_t<T>::zero().store(t_+i);
      }

    }


    /// @brief Vectorized first order push (var += step * derivative)
    /// @tparam T Type of the variables
    /// @param [in] a_ Step
    /// @param [in] b_ Derivatives
    /// @param [in,out] c_ Variables
    /// @param [in] n Size of the arrays
    template <class T>
    inline void push1stOrder (const T& a_, const T *b_, T *c_, const uint n) {

      vector_t<T> a(a_);
      vector_t<T> b, c;

      for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

      	b.load(b_+i);
      	c.load(c_+i);
    
      	c = fmadd(a, b, c);
      
      	c.store(c_+i);
    
      }

    }


    /// @brief Vectorized first order push (var += step * 1st_derivative + 0.5*step*step * 2nd_derivative)
    /// @tparam T Type of the variables
    /// @param [in] a_ Step
    /// @param [in] b_ Second derivatives
    /// @param [in] c_ First derivatives
    /// @param [in,out] d_ Variables
    /// @param [in] n Size of the arrays
    template <class T>
    inline void push2ndOrder (const T& a_, const T *b_, const T *c_, T *d_, const uint n) {

      vector_t<T> a1(a_);
      vector_t<T> a2(0.5*a_*a_);
      vector_t<T> b, c, d;

      for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

      	b.load(b_+i);
      	c.load(c_+i);
      	d.load(d_+i);
    
      	d = fmadd(a1, c, d);
      	d = fmadd(a2, b, d);
      
      	d.store(d_+i);
    
      }

    }


    /// @brief Vectorized Langevin fluctuation (\f$ var = exp(-\gamma*dt)*var + sqrt( (1-exp(-2*\gamma*dt)) * mass/beta ) * gaussian(0,1) \f$ )
    /// @tparam T Type of the variables
    /// @param [in] a_ \f$ -\gamma*\Delta t \f$
    /// @param [in,out] v_ Velocities
    /// @param [in] b_ Inverse temperature \f$\beta\f$
    /// @param [in] invMass_ Inverse masses
    /// @param [in] rnd1_ First random array
    /// @param [in] rnd2_ Second random array
    /// @param [in] n Size of the arrays
    template <class T>
    inline void pushFluctuationLangevin (const T& a_, T* v_, const T& b_, const T* invMass_, const T* rnd1_ , const T* rnd2_, const uint n) {

      vector_t<T> a(a_);
      vector_t<T> b(1./b_);
      
      vector_t<T> v, alpha, tmp;

      for (uint i=0; i<n; i+=vector_t<T>::chunk_size) {

      	v  .load(v_+i);
      	tmp.load(invMass_+i);

	alpha = exp(a);
	tmp = (vector_t<T>::one()-alpha*alpha)*tmp*b;
	  
      	v = alpha*v + sqrt(tmp)*random<T>::make_normal(rnd1_+i, rnd2_+i);

      	v.store(v_+i);

      }

    }


  }


}

#endif // __PUSH_HPP_INCLUDED
