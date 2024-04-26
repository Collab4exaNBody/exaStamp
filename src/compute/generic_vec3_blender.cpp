/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

//#pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/generic_vec3_blender.h>

namespace exanb
{
  template<class GridT> using BlendForceToFlatArray      = GenericVec3Blender< GridT, field::_fx, field::_fy, field::_fz , field::_flat_fx,field::_flat_fy,field::_flat_fz >;
  template<class GridT> using BlendForceFromFlatArray    = GenericVec3Blender< GridT, field::_flat_fx, field::_flat_fy, field::_flat_fz , field::_fx,field::_fy,field::_fz >;
  template<class GridT> using BlendPositionToFlatArray   = GenericVec3Blender< GridT, field::_rx, field::_ry, field::_rz , field::_flat_rx,field::_flat_ry,field::_flat_rz >;
  template<class GridT> using BlendPositionFromFlatArray = GenericVec3Blender< GridT, field::_flat_rx, field::_flat_ry, field::_flat_rz , field::_rx,field::_ry,field::_rz >;
  
 // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "blend_force_to_flat_array", make_grid_variant_operator< BlendForceToFlatArray > );
   OperatorNodeFactory::instance()->register_factory( "blend_position_to_flat_array", make_grid_variant_operator< BlendPositionToFlatArray > );

   OperatorNodeFactory::instance()->register_factory( "blend_force_from_flat_array", make_grid_variant_operator< BlendForceFromFlatArray > );
   OperatorNodeFactory::instance()->register_factory( "blend_position_from_flat_array", make_grid_variant_operator< BlendPositionFromFlatArray > );
  }

}

