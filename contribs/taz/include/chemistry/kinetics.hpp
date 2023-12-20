/// @file kinetics.hpp
///@brief Interface for chemical kinetics object

#ifndef __KINETICS_HPP_INCLUDED
#define __KINETICS_HPP_INCLUDED

#include <string>

#include "utils/kernelFunction/kernelFunction.hpp"

#include "utils/stampUnits.hpp"

/// @class IKinetics
/// @brief Interface for a @c Kinetics
class IKinetics {
  
public:
  /// @brief Kinetics type
  enum Type {
    SECOND_ORDER, ///< Second order
  };

  /// @brief Constructor
  /// @param [in] zAB_ Arrhenius prefactor (direct reaction: A->B)
  /// @param [in] zBA_ Arrhenius prefactor (reverse reaction: A<-B)
  /// @param [in] actAB_ Activation energy (direct reaction: A->B)
  /// @param [in] actBA_ Activation energy (reverse reaction: A<-B)
  IKinetics(double zAB_, double zBA_, double actAB_, double actBA_) :
    zAB(zAB_), zBA(zBA_), actAB(actAB_), actBA(actBA_), eexo(actBA_-actAB_) {}
  
  /// @brief Destructor
  virtual ~IKinetics() {}

  /// @brief Get Kinetics type
  virtual Type getType() const = 0;

  /// @brief Get Kinetics name
  virtual std::string getName() const = 0;

  /// @brief Get reaction rate for direct reaction (A->B)
  /// @param [in] T Temperature
  inline double directReactionRate(double T) const {
    return zAB*exp(-actAB/Stamp_Constant::boltzmann/T);
  }

  /// @brief Get reaction rate for reverse reaction (A->B)
  /// @param [in] T Temperature
  inline double reverseReactionRate(double T) const {
    return zAB*exp(-actBA/Stamp_Constant::boltzmann/T);
  }

  /// @brief Get exothermicity
  inline double getExothermicity() {
    return eexo;
  }

  /// @brief Compute kinetics
  /// @param [in] T Temperature
  /// @param [in] r Distance
  /// @param [in] progress0 Progress variable for first particle
  /// @param [in] progress1 Progress variable for second particle
  /// @param [in] K Kernel function
  virtual double computeKinetics(double T, double r, double progress0, double progress1, KernelFunction* K) const = 0;
  
protected:

  double zAB; ///< Arrhenius prefactor (direct reaction: A->B)
  double zBA; ///< Arrhenius prefactor (reverse reaction: A<-B)
  double actAB; ///< Activation energy (direct reaction: A->B)
  double actBA; ///< Activation energy (reverse reaction: A<-B)
  double eexo; ///< Exothermicity

};

#endif // __KINETICS_HPP_INCLUDED
