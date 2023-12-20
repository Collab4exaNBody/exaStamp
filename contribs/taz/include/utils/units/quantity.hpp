/// @file 
/// @brief Class quantity (a base for class Unit)

#ifndef __QUANTITY_HPP_INCLUDED
#define __QUANTITY_HPP_INCLUDED


#include <vector>
#include <iomanip>
#include <iostream>


/// @brief Number of basic physical dimensions
#define NUM_BASE_QUANTITIES 7


/// @brief Represents the dimension of a physical quantity
///
/// The dimension of a physical quantity can be expressed as a product
/// of the basic physical dimensions raised to a rational power
class Quantity {

public : 

	/// @brief Default constructor
	///
	///
  Quantity() : vec(NUM_BASE_QUANTITIES, 0) {}

  /// @brief Destructor (nothing to do)
  ~Quantity() {}

  Quantity(int L, int M, int T, int I, int Theta, int N, int J);
  Quantity(const Quantity& quantity);

  Quantity& operator = (const Quantity& quantity);

  int& operator [] (int index);
  const int& operator [] (int index) const;

  Quantity operator * (const Quantity& quantity) const;
  Quantity operator / (const Quantity& quantity) const;

  bool operator == (const Quantity& quantity) const;
  bool operator != (const Quantity& quantity) const;

  void print() const;

private :

  std::vector<int> vec; ///< Vector to store the powers to apply to the basic physical dimensions to get this dimension

};


/// @brief The seven basic physical dimensions and essential others
namespace Quantity_base {

  const Quantity one                 (0, 0, 0, 0, 0, 0, 0); ///< Dimensionless/constants

  const Quantity length              (1, 0, 0, 0, 0, 0, 0); ///< Length
  const Quantity mass                (0, 1, 0, 0, 0, 0, 0); ///< Mass
  const Quantity time                (0, 0, 1, 0, 0, 0, 0); ///< Time
  const Quantity electric_current    (0, 0, 0, 1, 0, 0, 0); ///< Electric current
  const Quantity temperature         (0, 0, 0, 0, 1, 0, 0); ///< Temperature
  const Quantity amount_of_substance (0, 0, 0, 0, 0, 1, 0); ///< Amount of substance
  const Quantity luminous_intensity  (0, 0, 0, 0, 0, 0, 1); ///< Luminous intensity

  const Quantity velocity = length/time; ///< Velocity
  const Quantity acceleration = velocity/time; ///< Acceleration

  const Quantity force = mass*acceleration; ///< Force
  const Quantity energy = length*force; ///< Energy
  const Quantity pressure = energy/length/length/length; ///< Pressure

}


/// @brief Constructor from the power on each basic dimension
/// @param [in] L Power on the length
/// @param [in] M Power on the mass
/// @param [in] T Power on time
/// @param [in] I Power on the electric current
/// @param [in] Theta Power on the temperature
/// @param [in] N Power on the amount of substance
/// @param [in] J Power on the luminous intensity
inline Quantity::Quantity(int L, int M, int T, int I, int Theta, int N, int J)
  : vec(NUM_BASE_QUANTITIES, 0) {
    vec[0] = L;
    vec[1] = M;
    vec[2] = T;
    vec[3] = I;
    vec[4] = Theta;
    vec[5] = N;
    vec[6] = J;
}


/// @brief Copy constructor
/// @param [in] quantity Another quantity
inline Quantity::Quantity(const Quantity& quantity) : vec(quantity.vec) {
}


/// @brief Assignment operator
/// @param [in] quantity Quantity to copy
inline Quantity& Quantity::operator = (const Quantity& quantity) {
  vec = quantity.vec;
  return *this;
}


/// @brief Access operator
///
/// Get the power on a basic dimension
/// @param [in] index Index of the dimension
inline int& Quantity::operator [] (int index) {
  return vec[index];
}


/// @brief Const access operator
///
/// Get the power on a basic dimension
/// @param [in] index Index of the basic dimension
inline const int& Quantity::operator [] (int index) const {
  return vec[index];
}


/// @brief Multiplication operator
///
/// Add the powers on each basic dimension
/// @param [in] quantity Dimension to multiply to
inline Quantity Quantity::operator * (const Quantity& quantity) const {
  Quantity temp(*this);
  for (unsigned int i=0; i<vec.size(); ++i) temp[i]+=quantity.vec[i];
  return temp;
}


/// @brief Division operator
///
/// Subtract the powers on each basic dimension
/// @param [in] quantity Dimension to divide by
inline Quantity Quantity::operator / (const Quantity& quantity) const {
  Quantity temp(*this);
  for (unsigned int i=0; i<vec.size(); ++i) temp[i]-=quantity.vec[i];
  return temp;
}


/// @brief Equality operator
/// @param [in] quantity Quantity to compare to
inline bool Quantity::operator == (const Quantity& quantity) const {
  return \
    vec[0] == quantity[0] &&
    vec[1] == quantity[1] &&
    vec[2] == quantity[2] &&
    vec[3] == quantity[3] &&
    vec[4] == quantity[4] &&
    vec[5] == quantity[5] &&
    vec[6] == quantity[6] ;
}


/// @brief Inequality operator
/// @param [in] quantity Quantity to compare to
inline bool Quantity::operator != (const Quantity& quantity) const {
  return !(*this==quantity);
}


/// @brief Print to debug
inline void Quantity::print() const {
  for (auto i : vec) std::cout<< std::setw(3) << i << " " ;
}

#endif // __QUANTITY_HPP_INCLUDED
