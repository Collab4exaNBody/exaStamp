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

// // DO NOT REMOVE THIS LINE

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/generic_vec3_copy.h>

namespace exanb
{
  template<class GridT> using CopyForceToFlatArray      = GenericVec3Copy< GridT, field::_fx, field::_fy, field::_fz , field::_flat_fx,field::_flat_fy,field::_flat_fz >;
  template<class GridT> using CopyForceFromFlatArray    = GenericVec3Copy< GridT, field::_flat_fx, field::_flat_fy, field::_flat_fz , field::_fx,field::_fy,field::_fz >;
  template<class GridT> using CopyPositionToFlatArray   = GenericVec3Copy< GridT, field::_rx, field::_ry, field::_rz , field::_flat_rx,field::_flat_ry,field::_flat_rz >;
  template<class GridT> using CopyPositionFromFlatArray = GenericVec3Copy< GridT, field::_flat_rx, field::_flat_ry, field::_flat_rz , field::_rx,field::_ry,field::_rz >;
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(generic_vec3_copy)
  {
   OperatorNodeFactory::instance()->register_factory( "copy_force_to_flat_array", make_grid_variant_operator< CopyForceToFlatArray > );
   OperatorNodeFactory::instance()->register_factory( "copy_position_to_flat_array", make_grid_variant_operator< CopyPositionToFlatArray > );

   OperatorNodeFactory::instance()->register_factory( "copy_force_from_flat_array", make_grid_variant_operator< CopyForceFromFlatArray > );
   OperatorNodeFactory::instance()->register_factory( "copy_position_from_flat_array", make_grid_variant_operator< CopyPositionFromFlatArray > );
  }

}

