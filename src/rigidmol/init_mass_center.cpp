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
#include <onika/math/basic_types_stream.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_operators.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/yaml/yaml_utils.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>

namespace exaStamp
{
  using namespace exanb;

  struct InitMassCenter : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );

    inline void execute () override final
    {
      for(unsigned int a=0;a<species->size();a++)
      {
        ParticleSpecie& mol = species->at(a);
        if( mol.m_rigid_atom_count < 1 )
        {
          lerr << "rigid atom count must be >= 1"<<std::endl;
          std::abort();
        }
        /*if( mol.m_rigid_atom_count != mol.m_rigid_atom_names.size() )
        {
          lerr << "insconsistent rigid atom names count. count="<<mol.m_rigid_atom_count<<", names=";
          for(const auto& s:mol.m_rigid_atom_names) lerr<<s<<" ";
          lerr<<std::endl;
          std::abort();          
        }*/
        
        Vec3d masscenter_pos { 0.,0.,0. };
        double totalmass = 0.;
        for(unsigned int i=0;i<mol.m_rigid_atom_count;i++)
        {
          size_t site_atom_type = mol.m_rigid_atoms[i].m_atom_type;
          if( site_atom_type >= species->size() )
          {
            lerr << "bad rigid molecule site atom type "<<site_atom_type <<std::endl;
            std::abort();
          }
          double site_mass = species->at( site_atom_type ).m_mass;
          masscenter_pos += mol.m_rigid_atoms[i].m_pos * site_mass;
          totalmass += site_mass;
        }
        masscenter_pos /= totalmass;
        for(unsigned int i=0;i<mol.m_rigid_atom_count;i++)
        {
          //int site_atom_type = mol.m_rigid_atoms[i].m_atom_type;
          mol.m_rigid_atoms[i].m_pos -= masscenter_pos;
        }
      }
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(init_mass_center)
  {
    OperatorNodeFactory::instance()->register_factory( "init_mass_center", make_simple_operator< InitMassCenter > );
  }

}

