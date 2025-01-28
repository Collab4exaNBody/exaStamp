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

//_cuda_enable // DO NOT REMOVE THIS LINE

#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/generic_scalar_zero.h>

namespace exanb
{
  template<class GridT> using ZeroFlatForceEnergy    = GenericScalarZero< GridT, field::_flat_fx, field::_flat_fy, field::_flat_fz, field::_flat_ep >;
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(generic_scalar_zero)
  {
   OperatorNodeFactory::instance()->register_factory( "zero_flat_force_energy", make_grid_variant_operator< ZeroFlatForceEnergy > );
  }

}

