/// @file
/// @brief Definition of ideal gas (as a potential)

#ifndef __IDEAL_GAS_HPP_INCLUDED
#define __IDEAL_GAS_HPP_INCLUDED


#include "potential/potential.hpp"


/// @brief Ideal gas interaction potential
///
/// Ideal gas interaction is no interaction so this potential just do nothing
class IdealGas : public IPotential {

public:

  /// @brief Default constructor
  IdealGas() {}

  /// @brief Destructor (nothing to do)
  virtual ~IdealGas() {}

  /// @brief Get the subclass of the potential
  /// @return Subclass
  virtual Traversal getTraversal() const { 
    return Traversal::IDEAL_GAS;
  }

  /// @brief Get the type of the potential
  /// @return Type
  virtual Type getType() const {
    return Type::IDEAL_GAS_; 
  }

  /// @brief Get the name of the potential
  /// @return Name
  virtual std::string getName() const {
    return "Ideal-Gas";
  }

  /// @brief Get the cost to apply the potential
  ///
	/// No interaction -> no cost
  /// @return Cost
  virtual double cost() const {
    return 0.;
  }

};

#endif // __IDEAL_GAS_HPP_INCLUDED
