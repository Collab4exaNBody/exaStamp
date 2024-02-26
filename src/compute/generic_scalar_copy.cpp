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

#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/generic_scalar_copy.h>

namespace exanb
{
  template<class GridT> using CopyTypeToFlatArray   = GenericScalarCopy< GridT, field::_type , field::_flat_type >;
  template<class GridT> using CopyTypeFromFlatArray = GenericScalarCopy< GridT, field::_flat_type , field::_type >;
  
 // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "copy_type_to_flat_array", make_grid_variant_operator< CopyTypeToFlatArray > );
   OperatorNodeFactory::instance()->register_factory( "copy_type_from_flat_array", make_grid_variant_operator< CopyTypeFromFlatArray > );
  }

}

