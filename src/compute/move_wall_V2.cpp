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
#include <onika/math/basic_types_stream.h>
#include <vector>
#include <iomanip>

namespace exaStamp
{
  using namespace exanb;

  class move_wall_v2 : public OperatorNode
  {
   
    ADD_SLOT( Vec3d  , init_normal , INPUT , Vec3d{1.0,0.0,0.0} );
    ADD_SLOT( double , init_offset , INPUT , 0.0 );
    ADD_SLOT( double , init_cutoff , INPUT , REQUIRED );
    ADD_SLOT( double , init_epsilon , INPUT , REQUIRED );
    ADD_SLOT( double , init_time , INPUT , REQUIRED );
    ADD_SLOT( double , init_velocity , INPUT , REQUIRED );
    ADD_SLOT( long , init_exponent , INPUT , 12 );
    ADD_SLOT( double , physical_time , INPUT , REQUIRED );

    // outputs for walll
    ADD_SLOT( Vec3d  , normal , OUTPUT );
    ADD_SLOT( double , offset , OUTPUT );
    ADD_SLOT( double , cutoff , OUTPUT );
    ADD_SLOT( double , epsilon , OUTPUT );
    ADD_SLOT( long , exponent , OUTPUT );
  public:

    inline void execute () override final
    {
      *normal = *init_normal;
      *cutoff = *init_cutoff;
      *epsilon = *init_epsilon;
      *offset = *init_offset;
      *exponent = *init_exponent;


      if( *physical_time >=  *init_time)
      {
        *offset = *init_offset + (*physical_time - *init_time) * (*init_velocity);
      }

      else
      {
        *offset = *init_offset;
      }
   	ldbg<<"offset="<<*offset<<std::endl;
    }

  };
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(move_wall_V2)
  {
    OperatorNodeFactory::instance()->register_factory( "move_wall_v2", make_simple_operator< move_wall_v2 > );
  }

}

