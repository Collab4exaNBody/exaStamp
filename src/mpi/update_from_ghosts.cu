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
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/mpi/update_from_ghosts.h>

namespace exaStamp
{
  using namespace exanb;
    
  // === register factory ===
  template<typename GridT> using UpdateForceEnergyFromGhosts = UpdateFromGhosts< GridT , FieldSet<field::_fx,field::_fy,field::_fz, field::_ep>, UpdateValueAdd >;
  template<typename GridT> using UpdateFlatForceEnergyFromGhosts = UpdateFromGhosts< GridT , FieldSet<field::_flat_fx,field::_flat_fy,field::_flat_fz, field::_flat_ep>, UpdateValueAdd >;
  template<typename GridT> using UpdateVirialForceEnergyFromGhosts = UpdateFromGhosts< GridT , FieldSet<field::_fx,field::_fy,field::_fz, field::_ep, field::_virial>, UpdateValueAdd >;
  template<typename GridT> using UpdateFromGhostsTestId = UpdateFromGhosts< GridT , FieldSet<field::_id>, UpdateValueAssertEqual >;

  ONIKA_AUTORUN_INIT(update_from_ghosts)
  {
    OperatorNodeFactory::instance()->register_factory( "update_force_energy_from_ghost", make_grid_variant_operator<UpdateForceEnergyFromGhosts> );
    OperatorNodeFactory::instance()->register_factory( "flat_force_energy_from_ghost", make_grid_variant_operator<UpdateFlatForceEnergyFromGhosts> );
    OperatorNodeFactory::instance()->register_factory( "update_virial_force_energy_from_ghost", make_grid_variant_operator<UpdateVirialForceEnergyFromGhosts> );
    OperatorNodeFactory::instance()->register_factory( "update_from_ghost_check_id", make_grid_variant_operator<UpdateFromGhostsTestId> );
  }

}

