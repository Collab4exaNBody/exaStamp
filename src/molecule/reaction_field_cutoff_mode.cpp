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
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>

namespace exaStamp
{
  using namespace exanb;

  class ReactionFieldCutoffModeNode : public OperatorNode
  {
    ADD_SLOT( std::string        , mode                 , INPUT , REQUIRED , DocString{"truncated, shifted_force or switching"} );
    ADD_SLOT( double             , rs                   , INPUT , 0.0 , DocString{"switching start (mode switching only)"} );
    ADD_SLOT( LJExp6RFMultiParms , potentials_for_pairs , INPUT_OUTPUT , REQUIRED );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Sets the treatment of the reaction field at its cutoff rc for all pair potentials of the molecular force field
(potentials_for_pairs), whatever their origin (FAtomes and DONNEES files, molecule dump, molecule_force_field).
Use it in setup_system after the reader. Applies to ljexp6rf_pc, ljexp6rf_pc_fp32 and intramolecular_pair_force.
- truncated : plain reaction field (default), the force is not zero at rc.
- shifted_force : the force at rc is subtracted at every distance, E(r) - (r-rc) E'(rc).
- switching : E(r) is multiplied by a quintic switching function between rs and rc. Forces are unchanged below rs.
The cutoff mode is not stored in molecule dumps : set it again when restarting from a dump.
)EOF";
    }

    inline void execute () override final
    {
      ReactionFieldCutoffMode m;
      m.mode = reaction_field_cutoff_mode_from_name( *mode );
      m.rs = *rs;
      if( m.mode == RF_CUTOFF_SWITCHING && m.rs <= 0.0 )
      {
        fatal_error() << "reaction_field_cutoff_mode: mode switching needs the switching start rs" << std::endl;
      }
      if( m.mode != RF_CUTOFF_SWITCHING && m.rs != 0.0 )
      {
        fatal_error() << "reaction_field_cutoff_mode: rs is only used with mode switching" << std::endl;
      }

      size_t n_rf = 0;
      for(auto& pot : potentials_for_pairs->m_potentials)
      {
        if( pot.m_params.m_rf.is_null() ) continue;
        reaction_field_set_cutoff_mode( pot.m_params.m_rf , m );
        ++ n_rf;
      }
      if( n_rf == 0 )
      {
        lerr << "reaction_field_cutoff_mode: no pair potential has a reaction field, nothing to do" << std::endl;
      }
      lout << "Reaction field cutoff mode = " << reaction_field_cutoff_mode_name(m.mode);
      if( m.mode == RF_CUTOFF_SWITCHING ) lout << " , rs = " << m.rs;
      lout << " ("<<n_rf<<" pair potentials)" << std::endl;
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(reaction_field_cutoff_mode)
  {
   OperatorNodeFactory::instance()->register_factory( "reaction_field_cutoff_mode", make_simple_operator<ReactionFieldCutoffModeNode> );
  }

}
