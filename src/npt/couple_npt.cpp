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
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_stream.h>
#include <exaStamp/npt/npt.h>
#include <exanb/core/domain.h>

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>

#include <exaStamp/compute/thermodynamic_state.h>

#include <iostream>
#include <string>

namespace exaStamp
{
  using namespace exanb;
  
  struct CoupleNPTNode : public OperatorNode
  {
    
    ADD_SLOT( NPTContext              , npt_ctx    , INPUT_OUTPUT );
    ADD_SLOT( ThermodynamicState      , thermodynamic_state , INPUT );

    inline void execute () override final
    {
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      static constexpr double conv_pressure = EXASTAMP_CONST_QUANTITY( 1.e4 * onika::physics::atomicMass * ( m^3 ) );
      static constexpr double pascal_to_bar = 1.e-5;
      
      Mat3d tensor = pascal_to_bar * sim_info.full_stress_tensor() * conv_pressure;
      double scalar = pascal_to_bar * sim_info.pressure_scal() * conv_pressure;
      double ave = 0.0;
      
      if (npt_ctx->pstyle == "ISO")
      	npt_ctx->p_current[0] = npt_ctx->p_current[1] = npt_ctx->p_current[2] = scalar;
      else if (npt_ctx->pcouple == "XYZ") {
      	ave = 1.0/3.0 * (tensor.m11 + tensor.m22 + tensor.m33);
      	npt_ctx->p_current[0] = npt_ctx->p_current[1] = npt_ctx->p_current[2] = ave;
      } else if (npt_ctx->pcouple == "XY") {
      	ave = 0.5 * (tensor.m11 + tensor.m22);
      	npt_ctx->p_current[0] = npt_ctx->p_current[1] = ave;
      	npt_ctx->p_current[2] = tensor.m33;
      } else if (npt_ctx->pcouple == "YZ") {
      	ave = 0.5 * (tensor.m22 + tensor.m33);
      	npt_ctx->p_current[1] = npt_ctx->p_current[2] = ave;
      	npt_ctx->p_current[0] = tensor.m11;
      } else if (npt_ctx->pcouple == "XZ") {
      	ave = 0.5 * (tensor.m11 + tensor.m33);
      	npt_ctx->p_current[0] = npt_ctx->p_current[2] = ave;
      	npt_ctx->p_current[1] = tensor.m22;
      } else {
      	npt_ctx->p_current[0] = tensor.m11;
      	npt_ctx->p_current[1] = tensor.m22;
      	npt_ctx->p_current[2] = tensor.m33;
      }

      if (npt_ctx->pstyle == "TRICLINIC") {
        npt_ctx->p_current[3] = tensor.m23;
        npt_ctx->p_current[4] = tensor.m13;
        npt_ctx->p_current[5] = tensor.m12;
      }
      
    }
  };
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(couple_npt)
  {
   OperatorNodeFactory::instance()->register_factory( "couple_npt", make_compatible_operator< CoupleNPTNode > );
  }

}
