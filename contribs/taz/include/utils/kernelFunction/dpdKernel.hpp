/// @file
/// @brief Lucy function

#ifndef __DPD_KERNEL_HPP_INCLUDED
#define __DPD_KERNEL_HPP_INCLUDED


#include "kernelFunction.hpp"


/// @brief Simple linear kernel for DPD and DPDE
class DPDKernel : public KernelFunction {
public:
  /// @brief Constructor
  /// @param [in] h_ Cut-off radius
  DPDKernel(double h_) : ih(1./h_) {}

  /// @brief Destructor
  virtual ~DPDKernel() {}

  /// @brief Get the kernel type
  /// @return Kernel type
  virtual inline Type getType() { return DPD_KERNEL; }

  /// @brief Get the inverse smoothing length
  /// @return Inverse smoothing length
  inline double getInvSmoothingLength() { return ih; }

  virtual inline double operator()(double r);

  virtual inline double F(double r);

private:
  double ih; ///< Inverse smoothing length
};

/// @brief Evaluate the kernel function
/// @param [in] r Distance
/// @return Kernel value at distance r
double DPDKernel::operator()(double r) {
  double umrh = 1.-r*ih;

  return umrh;
}

/// @brief Evaluate the kernel derivative
/// @param [in] r Distance
/// @return F defined by \f[ \nabla_r W (r) = -F(r)r \f]
double DPDKernel::F(double r) {
  return -ih;
}

#endif // __DPD_KERNEL_HPP_INCLUDED
