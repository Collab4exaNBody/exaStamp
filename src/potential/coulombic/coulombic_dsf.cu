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


#include <exanb/core/grid.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exaStamp/potential/coulombic/coulombic_pair_operator.h>
#include <exaStamp/potential/coulombic/dsf.h>

namespace exaStamp
{
  using namespace exanb;

  template< class GridT > using DsfPCTmpl = CoulombicPairPC< GridT , DsfParameters , DsfKernel >;

  template< class GridT > using DsfSelfTmpl = CoulombicSelfPC< GridT , DsfParameters , DsfKernel >;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_dsf)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_dsf" , make_grid_variant_operator<DsfPCTmpl> );
    OperatorNodeFactory::instance()->register_factory( "coulombic_dsf_self" , make_grid_variant_operator<DsfSelfTmpl> );
  }

}
