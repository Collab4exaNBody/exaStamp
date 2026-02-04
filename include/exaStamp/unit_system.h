/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

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

