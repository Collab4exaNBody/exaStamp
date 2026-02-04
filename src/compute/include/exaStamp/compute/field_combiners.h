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

#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types.h>
#include <exanb/compute/math_functors.h>
#include <onika/soatl/field_combiner.h>
#include <exaStamp/compute/physics_functors.h>

ONIKA_DECLARE_FIELD_COMBINER( exaStamp, VelocityNorm2Combiner , vnorm2 , exanb::Vec3Norm2Functor , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz )
ONIKA_DECLARE_FIELD_COMBINER( exaStamp, VelocityNormCombiner  , vnorm  , exanb::Vec3NormFunctor  , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz )

ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MonomatKineticEnergyCombiner  , mv2 , exaStamp::KineticEnergyFunctor , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz )
ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MultimatKineticEnergyCombiner , mv2 , exaStamp::KineticEnergyFunctor , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz , exanb::field::_type )

ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MonomatKineticEnergyTensorCombiner  , mv2tensor , exaStamp::KineticEnergyTensorFunctor , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz )
ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MultimatKineticEnergyTensorCombiner , mv2tensor , exaStamp::KineticEnergyTensorFunctor , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz , exanb::field::_type )

ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MonomatMomentumCombiner  , momentum , exaStamp::MomentumFunctor , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz )
ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MultimatMomentumCombiner , momentum , exaStamp::MomentumFunctor , exanb::field::_vx , exanb::field::_vy , exanb::field::_vz , exanb::field::_type )

ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MonomatMassCombiner  , mass , exaStamp::MassFunctor )
ONIKA_DECLARE_FIELD_COMBINER( exaStamp, MultimatMassCombiner , mass , exaStamp::MassFunctor , exanb::field::_type )

ONIKA_DECLARE_FIELD_COMBINER( exaStamp, AtomTypeChargeCombiner , type_charge , exaStamp::ChargeFunctor , exanb::field::_type )

