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
}

