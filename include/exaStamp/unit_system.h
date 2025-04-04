#pragma once

#include <onika/physics/units.h>
#include <onika/cuda/cuda.h>

namespace exaStamp
{

# define EXASTAMP_UNIT_SYSTEM_UNITS      \
    ::onika::physics::angstrom,          \
    ::onika::physics::atomic_mass_unit,  \
    ::onika::physics::picosecond,        \
    ::onika::physics::elementary_charge, \
    ::onika::physics::kelvin,            \
    ::onika::physics::particle,          \
    ::onika::physics::candela,           \
    ::onika::physics::radian
    
  static inline constexpr onika::physics::UnitSystem UNIT_SYSTEM = { { EXASTAMP_UNIT_SYSTEM_UNITS } };

  ONIKA_HOST_DEVICE_FUNC
  inline constexpr double to_internal_units( const onika::physics::Quantity & q )
  {
    constexpr onika::physics::UnitSystem target_units = { { EXASTAMP_UNIT_SYSTEM_UNITS } };
    return q.convert( target_units );
  }

}

#define EXASTAMP_QUANTITY( E ) ::exaStamp::to_internal_units( ONIKA_QUANTITY( E ) )
#define EXASTAMP_CONST_QUANTITY( E )  ::exaStamp::to_internal_units( ONIKA_CONST_QUANTITY( E ) )

