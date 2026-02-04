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

#include <exaStamp/particle_species/particle_specie_yaml.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/yaml/yaml_utils.h>
#include <exanb/core/particle_type_id.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>

namespace exaStamp
{
  using namespace exanb;

  struct SelectParticleSpecies : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( std::vector<std::string> , selection , INPUT , REQUIRED );
    ADD_SLOT(ParticleTypeMap , particle_type_map , INPUT_OUTPUT );

    inline void execute () override final
    {
      std::unordered_set<std::string> sel( selection->begin() , selection->end() );
      size_t count = 0;
      for(unsigned int i=0;i<species->size();i++)
      {
        if( sel.find(species->at(i).m_name) != sel.end() )
        {
          if(count!=i) species->at(count) = species->at(i);
          ++ count;
        }
      }
      ldbg << "selected "<<count<<" species among "<<species->size()<<std::endl;
      species->resize( count );

      // update particle type name map
      particle_type_map->clear();
      for(unsigned int a=0;a<species->size();a++)
      {
        (*particle_type_map) [ species->at(a).name() ] = a;
      }
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(select_species)
  {
    OperatorNodeFactory::instance()->register_factory( "select_species", make_simple_operator<SelectParticleSpecies> );
  }

}

