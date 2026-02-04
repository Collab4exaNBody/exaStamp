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

#include <onika/log.h>
#include <exaStamp/molecule/molecule_species.h>

namespace exaStamp
{
  using namespace exanb;
  
  class MoleculeSpeciesSetup : public OperatorNode
  {    
    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT_OUTPUT , DocString{"Molecule descriptions"} );

  public:
    inline bool is_sink() const override final { return true; }
    
    inline void execute ()  override final
    {
      for(const auto& mol:molecules->m_molecules)
      {
        ldbg << "+ molecule "<< mol.name() << std::endl;
      }
    }

  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_species)
  {
    OperatorNodeFactory::instance()->register_factory( "molecule_species", make_simple_operator< MoleculeSpeciesSetup > );
  }

}

