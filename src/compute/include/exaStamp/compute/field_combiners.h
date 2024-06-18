#pragma once

#include <exanb/fields.h>
#include <exanb/core/basic_types.h>
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

