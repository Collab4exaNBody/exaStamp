/// @file 
/// @brief Basic physical and mathematical constants

#ifndef __CONSTANTS_HPP_INCLUDED
#define __CONSTANTS_HPP_INCLUDED


#include<math.h>


/// @brief Physical constants in SI units
namespace ConstantPhys {

  const double atomicMass       = 1.66053904020e-27;  ///< Dalton atomic mass unit in kg
  const double elementaryCharge = 1.6021892e-19;      ///< Elementary charge in Coulomb
  const double boltzmann        = 1.380662e-23;       ///< Boltzmann constant in joules per kelvin
};


/// @brief Mathematical constants
namespace ConstantMath { 

  const double e   = 2.718281828;    ///< Euler's number
  const double phi = 1.61803399;     ///< Golden ratio
  const double pi  = 4. * atan(1.);  ///< Pi

};


/// @brief Common namespace for all the constants
namespace Constant {

  using namespace ConstantMath;
  using namespace ConstantPhys;

};

#endif // __CONSTANTS_HPP_INCLUDED
