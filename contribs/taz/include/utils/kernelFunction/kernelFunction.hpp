/// @file
/// @brief Kernel function for SDPD

#ifndef __KERNEL_FUNCTION_HPP_INCLUDED
#define __KERNEL_FUNCTION_HPP_INCLUDED


/// @brief Kernel function
class KernelFunction {
public:
  /// @brief Kernel types
  enum Type {
    LUCY, ///< Lucy kernel
    DPD_KERNEL, ///< Simple linear kernel
    NONE, ///< No kernel
  };

  /// @brief Destructor
  virtual ~KernelFunction() {}

  /// @brief Get Kernel type
  virtual Type getType() = 0;
  /// @brief Evaluate the kernel at distance r
  /// @param [in] r Distance
  virtual double operator()(double r) = 0;
  /// @brief Evaluate the kernel derivative
  /// @param [in] r distance
  virtual double F(double r) = 0;

};

#endif // __KERNEL_FUNCTION_HPP_INCLUDED
