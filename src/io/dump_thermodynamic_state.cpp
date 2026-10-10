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
#include <onika/print_utils.h>
#include <exaStamp/thermo_state/thermodynamic_state.h>
#include <onika/physics/constants.h>
#include <exanb/core/domain.h>

#include <exaStamp/io/thermodynamic_log_config.h>

#include <algorithm>
#include <sstream>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct DumpThermodynamicStateNode : public OperatorNode
  {
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( long               , timestep            , INPUT, REQUIRED);
    ADD_SLOT( double             , physical_time       , INPUT );
    ADD_SLOT( bool               , print_header        , INPUT , true );
    ADD_SLOT( bool               , internal_units      , INPUT , false );
    ADD_SLOT( std::string        , log_mode            , INPUT , "dump_default" , DocString{"'dump_default' (historical column set), any of print_thermodynamic_state's presets ('thermo_full', 'vol_fluct_tricl', 'mechanical', ...), or a ';'-separated list like 'stp;pht;mas;vol'"} );
    ADD_SLOT( std::string        , log_format          , INPUT , OPTIONAL , DocString{"optional ';'-separated printf-style format overrides applied positionally to the active columns, e.g. '%10.3f;%12.6e'"} );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT, REQUIRED);
    ADD_SLOT( double             , total_electronic_energy , INPUT, OPTIONAL );
    ADD_SLOT( double             , ion_transfer_energy     , INPUT, OPTIONAL );
    ADD_SLOT( std::string        , thermostate_file    , INPUT , "thermodynamic_state.csv" );
    ADD_SLOT( bool               , force_flush_file    , INPUT , false );
    ADD_SLOT( bool               , force_append_thermo , INPUT , false );
    ADD_SLOT( bool               , is_dump_virial      , INPUT , false);
    // NEW
    ADD_SLOT(Domain              , domain              , INPUT , OPTIONAL, DocString{"Deformation box matrix"} );

    ADD_SLOT(ThermodynamicLogConfig, log_config        , PRIVATE );

    inline void execute () override final
    {
      // MPI Initialization
      int rank = 0;
      MPI_Comm_rank(*mpi, &rank);

      // initialisation : remove output.csv
      if(rank!=0) { return; }

      if( log_config->m_active_items.empty() )
      {
        *log_config = thermodynamic_dump_config_default;
        ldbg << "log mode = "<< *log_mode << std::endl;
        thermodynamic_log_apply_mode( *log_config, *log_mode, "dump_thermodynamic_state" );

        // same duplicate-column guard as print_thermodynamic_state: only auto-append if the
        // user's own log_mode list doesn't already name "ele"/"ite" explicitly.
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

      bool is_dump_virial = *(this->is_dump_virial);
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);

      double el_energy = total_electronic_energy.has_value() ? *total_electronic_energy : 0.0;
      double ion_energy = ion_transfer_energy.has_value() ? *ion_transfer_energy : 0.0;

      double conv_temperature, conv_energy, conv_pressure, conv_density;
      thermodynamic_log_conversion_factors( *internal_units, conv_temperature, conv_energy, conv_pressure, conv_density );

      double values[ThermodynamicLogConfig::LOG_ITEM_COUNT];
      thermodynamic_log_fill_values( values, sim_info, *domain, *timestep, *physical_time, el_energy, ion_energy, 0.0, conv_temperature, conv_energy, conv_pressure, conv_density );

      // virial columns are not part of the LogItemId column system (a separate, legacy per-9-
      // component raw stress tensor dump) -- spliced onto the same header/data line as before.
      std::string virial_header_suffix, virial_data_suffix;
      if( is_dump_virial )
      {
        virial_header_suffix = "  S11  S12  S13  S21  S22  S23  S31  S32  S33";
        virial_data_suffix = onika::format_string(" % .7e  % .7e  % .7e % .7e  % .7e  % .7e % .7e  % .7e  % .7e",
          sim_info.stress_tensor().m11 * conv_pressure,
          sim_info.stress_tensor().m12 * conv_pressure,
          sim_info.stress_tensor().m13 * conv_pressure,
          sim_info.stress_tensor().m21 * conv_pressure,
          sim_info.stress_tensor().m22 * conv_pressure,
          sim_info.stress_tensor().m23 * conv_pressure,
          sim_info.stress_tensor().m31 * conv_pressure,
          sim_info.stress_tensor().m32 * conv_pressure,
          sim_info.stress_tensor().m33 * conv_pressure);
      }

      std::ostringstream oss;
      log_config->print_log( oss, values, *print_header, 0, 0, 0.0, false, virial_header_suffix, virial_data_suffix );

      onika::FileAppendWriteBuffer::instance().append_to_file( *thermostate_file , oss.str(), *force_append_thermo );

      if( *force_flush_file )
      {
        onika::FileAppendWriteBuffer::instance().flush();
      }

    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
writes thermodynmic values to a log file that can later be used analyzed.
one can plot some time related value using commands like the following :
gnuplot -e 'plot "thermodynamic_state.csv" every ::1 using 2:4' # this plots total energy over time
)EOF";
    }



  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(dump_thermodynamic_state)
  {
   OperatorNodeFactory::instance()->register_factory( "dump_thermodynamic_state", make_simple_operator<DumpThermodynamicStateNode> );
  }

}

