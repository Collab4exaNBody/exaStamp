/// @file
/// @brief Vectorization of "push" functions

#ifndef __SIMDAMR_PUSH_HPP_INCLUDED
#define __SIMDAMR_PUSH_HPP_INCLUDED


#include "saru/saru.hpp"
#include <cmath>
/// @brief Namespace for SIMD implementations
namespace simdAMR {


  /// @brief Namespace for vectorized computations
  namespace kernels {


    ///@brief Convert two independent sets of random integers in [-RAND_MAX, RAND_MAX[
    /// with an uniform distribution into random reals of type T in ]0, 1] with an normal distribution
    /// @tparam T Type of the variables
    /// @param [in] rnd1 Random integers
    /// @param [in] rnd2 Random integers
    /// @return Random reals of type T
    static inline double make_normal(const double rnd1, const double rnd2) {

      double r2 = Saru::RAND_MAX_PLS_TWO;
      double ir3 = Saru::INV_2_RAND_MAX_PLS_THREE;

      double u = (rnd1 + r2) * ir3;
      double v = (rnd2 + r2) * ir3;
      return std::sqrt( -2.0 * log(u) ) * std::cos( Saru::TWO_PI * v );
    }


    /// @brief Vectorized reset of an array
    /// @tparam T Type of the variables
    /// @param [in,out] t_ Array to reset
    /// @param [in] n Size of the array
    template <class T>
    inline void reset (T* t_, const uint n) {
      #pragma omp simd
      ////#pragma vector aligned
      for (uint i=0; i<n;++i) 
      	t_[i]=0.;

    }


    /// @brief Vectorized first order push (var += step * derivative)
    /// @tparam T Type of the variables
    /// @param [in] a_ Step
    /// @param [in] b_ Derivatives
    /// @param [in,out] c_ Variables
    /// @param [in] n Size of the arrays
    template <class T>
    inline void push1stOrder (const T& a_, const T *b_, T *c_, const uint n) {
      #pragma omp simd
      ////#pragma vector aligned
      for (uint i=0; i<n; ++i) 
        c_[i] += a_*b_[i];
    

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

      T a2(0.5*a_*a_);

      #pragma omp simd
      ////#pragma vector aligned
      for (uint i=0; i<n; ++i) {

        d_[i] += a_ * c_[i];
        d_[i] += a2 * b_[i];
    
      }

    }


    inline double make_normal_amr(const double rnd1, const double rnd2) {
      double u = rnd1;
      double v = rnd2; 
      u = (u + Saru::RAND_MAX_PLS_TWO) * Saru::INV_2_RAND_MAX_PLS_THREE ;
      v = (v + Saru::RAND_MAX_PLS_TWO) * Saru::INV_2_RAND_MAX_PLS_THREE ;
      return std::sqrt( -2.0 * std::log(u) ) * std::cos( Saru::TWO_PI * v );
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

      T b(1./b_);
      
      T alpha, tmp;

      for (uint i=0; i<n; ++i) {
	      alpha = std::exp(a_);
	      tmp   = (1-alpha*alpha)*invMass_[i]*b;
      	      v_[i] = alpha*v_[i] + std::sqrt(tmp)* make_normal(rnd1_[i], rnd2_[i]);
      }

    }


  }
}

#endif // __PUSH_HPP_INCLUDED
