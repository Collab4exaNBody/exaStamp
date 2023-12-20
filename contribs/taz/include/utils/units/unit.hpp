/// @file 
/// @brief Class Unit

#ifndef __UNIT_HPP_INCLUDED
#define __UNIT_HPP_INCLUDED


#include "utils/constants.hpp"
#include "utils/units/quantity.hpp"


/// @brief Class to handle units
///
/// A unit consists of a physical dimension and a value in SI units
class Unit {

public:

	/// @brief Destructor (nothing to do)
  ~Unit() {}

  /// @brief Constructor from a dimension and value
	/// @param [in] value_ Value in SI units
  /// @param [in] quantity_ Physical dimension
  Unit(double value_, const Quantity& quantity_) : ref(quantity_), value(value_) {} 
  /// @brief Copy constructor
  /// @param [in] unit Unit to copy
  Unit(const Unit& unit) : ref(unit.ref), value(unit.value) {}

  Unit& operator = (const Unit& unit);

  Unit operator * (double value_) const;
  Unit operator / (double value_) const;

  Unit operator * (const Unit& unit) const;
  Unit operator / (const Unit& unit) const;

  Unit inv() const;

  const Quantity& quantity() const;
  double getValue() const;

  void print() const;

private:

  Unit();

  Quantity ref; ///< Physical dimension
  double value; ///< Value in SI units

};


/// @brief Multiplication operator between a constant and a unit
/// @param [in] value Constant to multiply
/// @param [in] unit Unit to multiply
/// @return New unit of the same dimension
inline Unit operator * (double value, const Unit& unit) {
  return unit*value;
}


/// @brief Main units in International System (and other)
namespace SI_Units_base {


  using namespace Quantity_base;


  const Unit meter   (1., length); ///< Length SI unit : meter
  const Unit kilogram(1., mass); ///< Mass SI unit : kilogram
  const Unit second  (1., Quantity_base::time); ///< Time SI unit : second
  const Unit ampere  (1., electric_current); ///< Electric current SI unit : Ampere
  const Unit kelvin  (1., temperature); ///< Temperature SI unit : Kelvin
  const Unit mol     (1., amount_of_substance); ///< Amount of substance SI unit : mole
  const Unit candela (1., luminous_intensity); ///< Luminous intensity SI unit : candela

  const Unit coulomb = second * ampere; ///< Electric charge SI unit : Coulomb
  const Unit joule(1., energy); ///< Energy Si unit : Joule

  const Unit electronVolt(Constant::elementaryCharge, energy); ///< Electronvolt
  const Unit kcalPerMol(6.9477e-21, energy);	///< kcal/mol
  const Unit angstrom(1.0e-10,length);	///< Angstrom

  const Unit pascal(1., pressure); ///< Pressure Si unit : Pascal
  
}


/// @brief Assignment operator
/// @param [in] unit Unit to copy
inline Unit& Unit::operator = (const Unit& unit) {
  ref = unit.ref;
  value = unit.value;
  return *this;
}


/// @brief Multiplication operator to a constant
/// @param [in] value_ Constant to multiply to
/// @return New unit with the same dimension
inline Unit Unit::operator * (double value_) const {
  return Unit(value*value_, ref);
}


/// @brief Division operator by a constant
/// @param [in] value_ Constant to divide by
/// @return New unit with the same dimension
inline Unit Unit::operator / (double value_) const {
  return Unit(value/value_, ref);
}


/// @brief Multiplication operator to a unit
/// @param [in] unit Unit to multiply to
/// @return New unit with multiplied dimensions
inline Unit Unit::operator * (const Unit& unit) const {
  return Unit(value*unit.value, ref*unit.ref);
}


/// @brief Division operator to a unit
/// @param [in] unit Unit to divide by
/// @return New unit with divided dimensions
inline Unit Unit::operator / (const Unit& unit) const {
  return Unit(value/unit.value, ref/unit.ref);
}


/// @brief Inverse unit
/// @return New unit with inverse value and inverse dimension
inline Unit Unit::inv() const  {
  return Unit(1/value, Quantity_base::one/ref);
}


/// @brief Accessor to the dimension
inline const Quantity& Unit::quantity() const {
  return ref;
}


/// @brief Accessor to the value
inline double Unit::getValue() const {
  return value;
}


/// @brief Print for debug
inline void Unit::print() const {
  std::cout<< value << " | ";
  ref.print();
}

#endif // __UNIT_HPP_INCLUDED
