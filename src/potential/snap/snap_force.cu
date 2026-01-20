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

#include <md/snap/snap_force.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <onika/cpp_utils.h>

namespace exaStamp
{
# if defined(SNAP_FP32_MATH) || defined(SNAP_FP64_MATH)
  template<class GridT> using SnapForceXSTmpl = md::SnapForceGeneric<GridT,field::_ep,field::_virial>;
# else
  template<class GridT> using SnapForceXSTmpl = md::SnapForceGenericFP64<GridT,field::_ep,field::_virial>;
# endif

  // === register factories ===  
  ONIKA_AUTORUN_INIT(snap_force)
  {
    OperatorNodeFactory::instance()->register_factory( "snap_force" ,make_grid_variant_operator< SnapForceXSTmpl > );
  }

}


