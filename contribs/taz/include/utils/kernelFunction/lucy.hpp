/// @file
/// @brief Lucy function

#ifndef __LUCY_HPP_INCLUDED
#define __LUCY_HPP_INCLUDED


#include "kernelFunction.hpp"


/// @brief Lucy kernel function
class Lucy : public KernelFunction {
public:
  /// @brief Constructor
  /// @param [in] h_ Cut-off radius
  Lucy(double h_) : ih(1./h_), aW(105./(16*M_PI*h_*h_*h_)), aF(315./(4*M_PI*pow(h_,5))) {}

  /// @brief Destructor
  virtual ~Lucy() {}

  /// @brief Get the kernel type
  /// @return Kernel type
  virtual inline Type getType() { return LUCY; }

  /// @brief Get the inverse smoothing length
  /// @return Inverse smoothing length
  inline double getInvSmoothingLength() { return ih; }

  /// @brief Get the normalization factor of the kernel W
  /// @return Kernel normalisation factor
  inline double getAlphaW() { return aW; }

  /// @brief Get the normalization factor of the kernel derivative F
  /// @return Normalisation factor for F
  inline double getAlphaF() { return aF; }

  virtual double operator()(double r);

  virtual double F(double r);

private:
  double ih; ///< inverse smoothing length
  double aW; ///< normalization constant for the kernel
  double aF; ///<normalization constant for F = -W'/r
};

/// @brief Evaluate the kernel function
/// @param [in] r Distance
/// @return Kernel value at distance r
inline double Lucy::operator()(double r) {
  double rh = r*ih;
  double umrh = 1-rh;

  return aW*(1+3*rh)*umrh*umrh*umrh;
}

/// @brief Evaluate the kernel derivative
/// @param [in] r Distance
/// @return F defined by \f[ \nabla_r W (r) = -F(r)r \f]
inline double Lucy::F(double r) {
  double umrh = 1-r*ih;

  return aF*umrh*umrh;
}

#endif // __LUCY_HPP_INCLUDED
