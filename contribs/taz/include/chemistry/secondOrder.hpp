///@file secondOrder.hpp
///@brief Second order kinetics
#ifndef __SECOND_ORDER_HPP
#define __SECOND_ORDER_HPP


#include "chemistry/kinetics.hpp"
#include "simd/chemistry.hpp"


/// @class SecondOrder
/// @brief Class for second order chemical kinetics
class SecondOrder : public IKinetics {
  
public:
  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that kinetics
  /// @brief Shorcut for the simd version of the kinetics
  typedef simd::kernels::SecondOrderReaction<double> simd_opt_t;

  
  /// @brief Constructor
  /// @param [in] zAB_ Arrhenius prefactor (direct reaction: A->B)
  /// @param [in] zBA_ Arrhenius prefactor (reverse reaction: A<-B)
  /// @param [in] actAB_ Activation energy (direct reaction: A->B)
  /// @param [in] actBA_ Activation energy (reverse reaction: A<-B)
  SecondOrder(double zAB_, double zBA_, double actAB_, double actBA_) : IKinetics(zAB_,zBA_,actAB_,actBA_) {}

  
  /// @brief Destructor
  virtual ~SecondOrder() {}

  /// @brief Get Kinetics type
  virtual Type getType() const {
    return Type::SECOND_ORDER;
  }

  /// @brief Get Kinetics name
  virtual std::string getName() const {
    return "Second order";
  }

  /// @brief Set the parameters in the simd version of the kinetics
  void setParameters(simd::kernels::SecondOrderReaction<double>& opt) {
    opt.setParameters(zAB,zBA,actAB,actBA);
  }
  
  /// @brief Compute kinetics
  /// @param [in] T Temperature
  /// @param [in] r Distance
  /// @param [in] progress0 Progress variable for first particle
  /// @param [in] progress1 Progress variable for second particle
  /// @param [in] K Kernel function
  virtual double computeKinetics(double T, double r, double progress0, double progress1, KernelFunction* K) const {
    return (*K)(r) * (directReactionRate(T)*(1-progress0)*(1-progress1) - reverseReactionRate(T)*progress0*progress1);
  }
};

#endif
