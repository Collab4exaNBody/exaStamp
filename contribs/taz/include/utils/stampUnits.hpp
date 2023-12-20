/// @file 
/// @brief Definition of the unit used in ExaStamp (and Stamp) and associated functions

#ifndef __STAMP_UNITS_HPP_INCLUDED
#define __STAMP_UNITS_HPP_INCLUDED


#include "utils/constants.hpp"
#include "utils/units/unit.hpp"


/// @brief Units used in the code in order to have values closer to 1
namespace Stamp_Units {


  using namespace SI_Units_base;


  const Unit length = 1.0e-09 * meter; ///< Stamp length : nanometer
  const Unit mass = Constant::atomicMass * kilogram; ///< Stamp mass : Dalton atomic mass
  const Unit time = 1.0e-12 * second; ///< Stamp time : picosecond
  const Unit temperature = 1.0e+00 * kelvin; ///< Stamp temperature : kelvin
  const Unit charge = Constant::elementaryCharge*coulomb; ///< Stamp charge : elementary charge


  //const Unit timeInverse = time.inv();


  const Unit speed = length / time; ///< Stamp speed
  const Unit force = mass*length/(time*time); ///< Stamp force
  const Unit energy = length*force; ///< Stamp energy
  const Unit pressure = energy / length / length /length; ///< Stamp pressure  

  // to be continued ...
  

}


/// @brief Convert a value with a unit in another unit
///
/// In case of wrong units (ie they do not have the same reference Quantity, the
/// return value is 0)
/// @param [in] val Value to convert
/// @param [in] valUnit Initial unit
/// @param [in] convertUnit Final unit
/// @return Value in the new unit
inline double convert(double val, const Unit& valUnit, const Unit& convertUnit) {

  if (valUnit.quantity()==convertUnit.quantity()) \
    return val*valUnit.getValue()/convertUnit.getValue();

  else return 0.;

}


/// @brief Physical constants in ExaStamp units
namespace Stamp_Constant {


  const double boltzmann = convert(Constant::boltzmann, SI_Units_base::joule/SI_Units_base::kelvin, Stamp_Units::energy/Stamp_Units::temperature); ///< Boltzmann constant in ExaStamp units


}

#endif // __STAMP_UNITS_HPP_INCLUDED
