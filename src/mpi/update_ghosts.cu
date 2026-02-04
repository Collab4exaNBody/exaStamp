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
#include <exanb/mpi/update_ghosts.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace UpdateGhostsUtils;

  // === register factory ===
  template<typename GridT> using UpdateGhostsRandVandVir = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz,  field::_vx, field::_vy, field::_vz, field::_virial > , false >;
  template<typename GridT> using UpdateGhostsRandV = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz,  field::_vx, field::_vy, field::_vz > , false >;

  template<typename GridT> using UpdateGhostsRandRf = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz, field::_rxf, field::_ryf, field::_rzf > , false >;  
  template<typename GridT> using UpdateGhostsRandRfandV = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz, field::_rxf, field::_ryf, field::_rzf, field::_vx, field::_vy, field::_vz > , false >;
  
  template<typename GridT> using UpdateGhostsRQ = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz , field::_orient > , false >;
  template<typename GridT> using UpdateGhostsIdMol = UpdateGhostsNode< GridT , FieldSet<field::_idmol> , false >;

  ONIKA_AUTORUN_INIT(update_ghosts)
  {
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_v_vir",   make_grid_variant_operator<UpdateGhostsRandVandVir> );
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_v",       make_grid_variant_operator<UpdateGhostsRandV> );
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_rf",      make_grid_variant_operator<UpdateGhostsRandRf> );    
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_rf_v",    make_grid_variant_operator<UpdateGhostsRandRfandV> );    
    OperatorNodeFactory::instance()->register_factory( "ghost_update_rq",        make_grid_variant_operator<UpdateGhostsRQ> );
    OperatorNodeFactory::instance()->register_factory( "ghost_update_idmol",     make_grid_variant_operator<UpdateGhostsIdMol> );
  }

}

