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

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/physics/constants.h>
#include <exaStamp/thermo_state/thermodynamic_state.h>

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  // These exist purely for human-readable progress printing (see cg_progress_screen in
  // config_conjugate_gradient.msp) - the CG algorithm itself (Armijo bound, direction, force_max
  // convergence test) works entirely in internal units and must NOT use these, converting only
  // energy or only force there would silently break the Armijo formula's dimensional consistency.
  //
  // Same conversion factor as print_thermodynamic_state.cpp (internal units -> eV): matching it
  // exactly here rather than re-deriving/hardcoding the constant in YAML.
  static inline double internal_to_ev_factor()
  {
    return 1.e4 * onika::physics::atomicMass / onika::physics::elementaryCharge;
  }

  // Bare scalar internal-units -> eV conversion, e.g. for force_max (already expressed "per
  // angstrom" since the internal length unit is angstrom, so no separate length conversion needed).
  class InternalToEVNode : public OperatorNode
  {
    ADD_SLOT( double , value    , INPUT , REQUIRED );
    ADD_SLOT( double , value_ev , OUTPUT );

  public:
    inline void execute () override final
    {
      *value_ev = (*value) * internal_to_ev_factor();
    }
  };

  // potential_energy() per atom, in eV - mirrors print_thermodynamic_state.cpp's POTENTIAL_E
  // column exactly, so it can be compared directly against it.
  class ThermoPotentialEnergyEVPerAtomNode : public OperatorNode
  {
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT , REQUIRED );
    ADD_SLOT( double             , value_ev            , OUTPUT );

  public:
    inline void execute () override final
    {
      *value_ev = thermodynamic_state->potential_energy() / thermodynamic_state->particle_count() * internal_to_ev_factor();
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(thermo_energy_display)
  {
    OperatorNodeFactory::instance()->register_factory( "internal_to_ev", make_compatible_operator< InternalToEVNode > );
    OperatorNodeFactory::instance()->register_factory( "thermo_potential_energy_ev_per_atom", make_compatible_operator< ThermoPotentialEnergyEVPerAtomNode > );
  }

}
