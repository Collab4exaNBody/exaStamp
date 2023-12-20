/// @file
/// @brief Definition of a smoothstep functions

#ifndef SPLINE_HPP_INCLUDED
#define SPLINE_HPP_INCLUDED


/// @brief Namespace for non-vectorized smoothstep functions
namespace Spline {


	/// @brief Template gathering smoothstep functions
	///
	/// The smoothstep function n is a 2n+1 order polynomial function
	/// that takes value 0 at 0 and 1 at 1 and whose n first derivatives
	/// take value 0 at 0 and 1
	/// @tparam N Number of the smoothstep function
  template <uint8_t N> double  S(double x);


	/// @brief Template gathering the first derivative of smoothstep functions
	/// @tparam N Number of the smoothstep function
  template <uint8_t N> double dS(double x);


	/// @brief Smoothstep 0 : identity (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double S<0>(double x) {
    return x;
  }


	/// @brief Smoothstep 1 (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double S<1>(double x) {
    return x*x * ( -2*x + 3 );
  }


	/// @brief Smoothstep 2 : smootherstep (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double S<2>(double x) {
    double x2 = x*x;
    return x2*x * ( 6*x2 - 15*x + 10 );
  }


	/// @brief Smoothstep 3 : smootheststep
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double S<3>(double x) {
    double x2 = x*x;
    return x2*x2 * ( -20*x2*x + 70*x2 - 84*x + 35 );
  }


	/// @brief Smoothstep 4 (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double S<4>(double x) {
    double x2 = x *x;
    double x3 = x2*x;
    return x3*x2 * ( 70*x2*x2 - 315*x3 + 540*x2 -420*x + 126 );
  }


	/// @brief Derivative of the smoothstep 0 (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double dS<0>(double x) {
    return 1;
  }


	/// @brief Derivative of the smoothstep 1 (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double dS<1>(double x) {
    return 6*x * ( -1*x + 1 );
  }


	/// @brief Derivative of the smoothstep 2 (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double dS<2>(double x) {
    double x2 = x*x;
    return 30*x2 * ( x2 - 2*x + 1 );
  }


	/// @brief Derivative of the smoothstep 3
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double dS<3>(double x) {
    double x2 = x *x;
    double x3 = x2*x;
    return 140*x3 * ( -1*x3 + 3*x2 - 3*x + 1 );
  }


	/// @brief Derivative of the smoothstep 4 (not used)
	/// @param [in] x Input (only values between 0 and 1 are relevant)
	/// @return Output
  template <> 
  inline double dS<4>(double x) {
    double x2 = x *x;
    double x3 = x2*x;
    double x4 = x2*x2;
    return 630*x4 * ( x4 - 4*x3 + 6*x2 -4*x + 1 );
  }


}

#endif // SPLINE_HPP_INCLUDED
