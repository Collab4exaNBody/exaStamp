#pragma once

#include <onika/physics/units.h>

namespace exaStamp
{
  static inline constexpr onika::physics::UnitSystem UNIT_SYSTEM = { {
    onika::physics::angstrom,
    onika::physics::atomic_mass_unit,
    onika::physics::picosecond,
    onika::physics::elementary_charge,
    onika::physics::kelvin,
    onika::physics::particle,
    onika::physics::candela,
    onika::physics::radian
    } };
    
  inline constexpr double to_internal_units( const onika::physics::Quantity & q )
  {
    return q.convert( UNIT_SYSTEM );
  }
    
}

#define EXASTAMP_QUANTITY( E ) ONIKA_QUANTITY( E ).convert( ::exaStamp::UNIT_SYSTEM )
#define EXASTAMP_CONST_QUANTITY( E ) ONIKA_CONST_QUANTITY( E ).convert( ::exaStamp::UNIT_SYSTEM )

