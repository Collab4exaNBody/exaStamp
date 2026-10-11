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
#include <onika/string_utils.h>
#include <exaStamp/thermo_state/thermodynamic_state.h>
#include <onika/physics/constants.h>
#include <exanb/core/domain.h>

#include <exaStamp/io/thermodynamic_log_config.h>

#include <algorithm>

namespace exaStamp
{
  using namespace exanb;

  class PrintThermodynamicStateNode : public OperatorNode
  {  
    // thermodynamic state & physics data
    ADD_SLOT( long               , timestep            , INPUT , REQUIRED );
    ADD_SLOT( double             , physical_time       , INPUT , REQUIRED );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT , REQUIRED );

    // printing options
    ADD_SLOT( bool               , print_header        , INPUT, false );
    ADD_SLOT( bool               , internal_units      , INPUT, false );
    ADD_SLOT( std::string        , log_mode            , INPUT, "default" , DocString{"'default', 'thermo_basic', 'thermo_full', 'vol_fluct_ortho', 'vol_fluct_tricl', 'mechanical', or a list like 'stp;pht;prt;sta;toe;kie;poe;vol'"} );
    ADD_SLOT( std::string        , log_format          , INPUT, OPTIONAL , DocString{"optional ';'-separated printf-style format overrides applied positionally to the active columns, e.g. '%10.3f;%12.6e'"} );

    // LB and particle movement statistics
    ADD_SLOT( long               , lb_counter          , INPUT_OUTPUT );
    ADD_SLOT( long               , move_counter        , INPUT_OUTPUT );
    ADD_SLOT( long               , domain_ext_counter  , INPUT_OUTPUT );
    ADD_SLOT( double             , lb_inbalance_max    , INPUT_OUTPUT );

    // optional physics quantities
    ADD_SLOT( double             , total_electronic_energy , INPUT, OPTIONAL );
    ADD_SLOT( double             , ion_transfer_energy     , INPUT, OPTIONAL );

    // NEW
    ADD_SLOT(Domain              , domain              , INPUT , OPTIONAL, DocString{"Deformation box matrix"} );

    ADD_SLOT(ThermodynamicLogConfig, log_config        , PRIVATE );
  public:
    inline bool is_sink() const override final { return true; }
  
    inline void execute () override final
    {
      if( log_config->m_active_items.empty() )
      {
        *log_config = thermodynamic_log_config_default;
        ldbg << "log mode = "<< *log_mode << std::endl;
        thermodynamic_log_apply_mode( *log_config, *log_mode, "print_thermodynamic_state" );

        // one-time setup, same as the log_mode column list above: appending this on every
        // execute() call (as opposed to just once here) made m_active_items -- and so the
        // printed columns -- grow by one ELECTRON_E entry every single print step. Also guarded
        // against a user's own log_mode list already naming "ele"/"ite" explicitly.
        if( total_electronic_energy.has_value()
            && std::find(log_config->m_active_items.begin(), log_config->m_active_items.end(), ThermodynamicLogConfig::ELECTRON_E) == log_config->m_active_items.end() )
        {
          log_config->m_active_items.push_back( ThermodynamicLogConfig::ELECTRON_E );
        }
        if( ion_transfer_energy.has_value()
            && std::find(log_config->m_active_items.begin(), log_config->m_active_items.end(), ThermodynamicLogConfig::ION_TRANSFER_E) == log_config->m_active_items.end() )
        {
          log_config->m_active_items.push_back( ThermodynamicLogConfig::ION_TRANSFER_E );
        }

        if( log_format.has_value() )
        {
          thermodynamic_log_apply_format( *log_config, *log_format );
        }
      }

      double el_energy = total_electronic_energy.has_value() ? *total_electronic_energy : 0.0;
      double ion_energy = ion_transfer_energy.has_value() ? *ion_transfer_energy : 0.0;

      double conv_temperature, conv_energy, conv_pressure, conv_density;
      thermodynamic_log_conversion_factors( *internal_units, conv_temperature, conv_energy, conv_pressure, conv_density );

      bool lb_flag = (*lb_counter) > 0 ;
      long move_count = *move_counter ;
      long domext_count = *domain_ext_counter;
      double lb_inbalance = *lb_inbalance_max;

      *lb_counter = 0;
      *move_counter = 0;
      *domain_ext_counter = 0;
      *lb_inbalance_max = 0.0;

      double values[ThermodynamicLogConfig::LOG_ITEM_COUNT];
      thermodynamic_log_fill_values( values, *thermodynamic_state, *domain, *timestep, *physical_time, el_energy, ion_energy, lb_inbalance, conv_temperature, conv_energy, conv_pressure, conv_density );

      log_config->print_log( lout , values, *print_header , move_count , domext_count , lb_inbalance , lb_flag );
    }

  };
    
  // === register factories ===  
  ONIKA_AUTORUN_INIT(print_thermodynamic_state)
  {
   OperatorNodeFactory::instance()->register_factory( "print_thermodynamic_state", make_simple_operator<PrintThermodynamicStateNode> );
  }

}

