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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <onika/log.h>
#include <onika/string_utils.h>
#include <onika/physics/constants.h>
#include <exaStamp/thermo_state/thermodynamic_state.h>
#include <exanb/core/domain.h>

namespace exaStamp
{
  using namespace exanb;

  struct ThermodynamicLogConfigItem
  {
    const char* m_name = "";
    std::string m_header = "";
    std::string m_format = "%g";
    long len = -1;
  };

  struct ThermodynamicLogConfig
  {
    enum LogItemId
    {
      TIME_STEP = 0 ,
      PHYS_TIME ,
      NB_PARTICLES ,
      SIM_STATUS ,
      TOTAL_E ,
      KINETIC_E ,
      POTENTIAL_E ,
      ELECTRON_E ,
      ION_TRANSFER_E ,
      TEMPERATURE ,
      Tx ,
      Ty ,
      Tz ,
      PRESSURE ,
      // full stress tensor (virial + kinetic), symmetric, 6 independent components
      Pxx ,
      Pyy ,
      Pzz ,
      Pxy ,
      Pxz ,
      Pyz ,
      // pure virial tensor (no kinetic term), symmetric, 6 independent components
      Vxx ,
      Vyy ,
      Vzz ,
      Vxy ,
      Vxz ,
      Vyz ,
      SMISES ,
      VOLUME ,
      MASS ,
      BOX_A ,
      BOX_B ,
      BOX_C ,
      BOX_ALPHA ,
      BOX_BETA ,
      BOX_GAMMA ,
      DENSITY,

      LOG_ITEM_COUNT
    };

    // header_suffix/data_suffix let a caller splice extra, hand-formatted columns (not part of
    // LogItemId) onto the same header/data line before its trailing newline -- used by
    // dump_thermodynamic_state's optional virial columns.
    template<class StreamT>
    inline void print_log( StreamT & out , const double values[LOG_ITEM_COUNT] , bool header, long move_count , long domext_count , double lb_inbalance, bool lb_flag, const std::string& header_suffix = "", const std::string& data_suffix = "" )
    {
      for(auto & p : m_avail_items)
      {
        if(p.second.len==-1)
        {
          std::string s = onika::format_string(p.second.m_format,0.0);
          size_t pfxlen = ( p.first == SIM_STATUS ) ? 3 : 0;
          p.second.len = std::max( s.length() + pfxlen ,  p.second.m_header.length() );
        }
      }
      if( header )
      {
        for(auto id : m_active_items) out << m_sep << onika::format_string("%*s", m_avail_items[id].len, m_avail_items[id].m_header );
        out << header_suffix;
        out << std::endl;
      }
      for(auto id : m_active_items)
      {
        if( id == SIM_STATUS )
        {
          char lb_move_char = ' ';
          if( move_count >= 1 )
          {
            if( move_count == 1 ) { lb_move_char = 'm'; }
            else if( move_count < 10 )  { lb_move_char = '0'+move_count; }
            else { lb_move_char = 'M'; }
          }
          char domext_char = ' ';
          if( domext_count >= 1 )
          {
            if( domext_count == 1 ) { domext_char = 'd'; }
            else if( domext_count < 10 )  { domext_char = '0'+domext_count; }
            else { domext_char = 'D'; }
          }
          std::string lb_value;
          if( lb_flag )
          {
            if( lb_inbalance == 0.0 )
            {
              lb_value = " N/A ";
            }
            else
            {
              lb_value = onika::format_string(m_avail_items[SIM_STATUS].m_format, values[SIM_STATUS] );
            }
          }

          out << m_sep << onika::format_string("%c %c%*s",lb_move_char,domext_char,m_avail_items[SIM_STATUS].len-3,lb_value);
        }
        else
        {
          out << m_sep << onika::format_string(m_avail_items[id].m_format,values[id]);
        }
      }
      out << data_suffix;
      out << std::endl;
    }

    std::map< LogItemId ,ThermodynamicLogConfigItem  > m_avail_items;
    std::vector<LogItemId> m_active_items;
    const char* m_sep = " ";
  };

  static const ThermodynamicLogConfig thermodynamic_log_config_default =
    {
      {
        { ThermodynamicLogConfig::TIME_STEP    , { "stp" , "Step"             , "% 9.0f" } } ,
        { ThermodynamicLogConfig::PHYS_TIME    , { "pht" , "Time (ps)"        , "% .6e" } } ,
        { ThermodynamicLogConfig::NB_PARTICLES , { "prt" , "Particles"        , "% 12.0f" } } ,
        { ThermodynamicLogConfig::SIM_STATUS   , { "sta" , "Mv/Ext/Imb."      , "%.1e" } } ,
        { ThermodynamicLogConfig::TOTAL_E      , { "toe" , "Tot. E. (eV/part)" , "% .10e" } } ,
        { ThermodynamicLogConfig::KINETIC_E    , { "kie" , "Kin. E. (eV/part)" , "% .10e" } } ,
        { ThermodynamicLogConfig::POTENTIAL_E  , { "poe" , "Pot. E. (eV/part)" , "% .10e" } } ,
        { ThermodynamicLogConfig::ELECTRON_E   , { "ele" , "Elec. E. (eV)"     , "% .10e" } } ,
        { ThermodynamicLogConfig::ION_TRANSFER_E, { "ite" , "Ion Transf. E. (eV)", "% .10e" } } ,
        { ThermodynamicLogConfig::TEMPERATURE  , { "tmp" , "Temp. (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::Tx           , { "tmx" , "   Tx (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::Ty           , { "tmy" , "   Ty (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::Tz           , { "tmz" , "   Tz (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::PRESSURE     , { "pre" , "Press. (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pxx          , { "pxx" , "   Pxx (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pyy          , { "pyy" , "   Pyy (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pzz          , { "pzz" , "   Pzz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pxy          , { "pxy" , "   Pxy (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pxz          , { "pxz" , "   Pxz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pyz          , { "pyz" , "   Pyz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Vxx          , { "vxx" , "   Vxx (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Vyy          , { "vyy" , "   Vyy (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Vzz          , { "vzz" , "   Vzz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Vxy          , { "vxy" , "   Vxy (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Vxz          , { "vxz" , "   Vxz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Vyz          , { "vyz" , "   Vyz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::SMISES       , { "smi" , "sMises (Pa)"      , "% .3e" } } ,
        { ThermodynamicLogConfig::VOLUME       , { "vol" , "Vol. (ang^3)"     , "% 13.6e" } } ,
        { ThermodynamicLogConfig::MASS         , { "mas" , "Mass"             , "% .3e" } } ,
        { ThermodynamicLogConfig::BOX_A        , { "bxa" , "A (ang)"          , "% .3f" } } ,
        { ThermodynamicLogConfig::BOX_B        , { "bxb" , "B (ang)"          , "% .3f" } } ,
        { ThermodynamicLogConfig::BOX_C        , { "bxc" , "C (ang)"          , "% .3f" } } ,
        { ThermodynamicLogConfig::BOX_ALPHA    , { "baa" , "alpha (deg)"      , "% 12.3f" } } ,
        { ThermodynamicLogConfig::BOX_BETA     , { "bab" , "beta (deg)"       , "% 12.3f" } } ,
        { ThermodynamicLogConfig::BOX_GAMMA    , { "bag" , "gamma (deg)"      , "% 12.3f" } } ,
        { ThermodynamicLogConfig::DENSITY      , { "rho" , "Rho (g/cm^3)"     , "% 12.6f" } }
      }
      ,
      {
        ThermodynamicLogConfig::TIME_STEP,
        ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::NB_PARTICLES,
        ThermodynamicLogConfig::SIM_STATUS,
        ThermodynamicLogConfig::TOTAL_E,
        ThermodynamicLogConfig::KINETIC_E,
        ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE,
        ThermodynamicLogConfig::PRESSURE,
        ThermodynamicLogConfig::VOLUME
      }
    };

  // Same column set/names as the screen config above, but with higher-precision formats suited to
  // a file meant for numeric post-processing (gnuplot, scripts) rather than terminal display --
  // used as dump_thermodynamic_state's own starting point instead of thermodynamic_log_config_default.
  static const ThermodynamicLogConfig thermodynamic_dump_config_default =
    {
      {
        { ThermodynamicLogConfig::TIME_STEP    , { "stp" , "Step"             , "% 9.0f" } } ,
        { ThermodynamicLogConfig::PHYS_TIME    , { "pht" , "Time (ps)"        , "% .6e" } } ,
        { ThermodynamicLogConfig::NB_PARTICLES , { "prt" , "Particles"        , "% 12.0f" } } ,
        { ThermodynamicLogConfig::SIM_STATUS   , { "sta" , "Mv/Ext/Imb."      , "%.1e" } } ,
        { ThermodynamicLogConfig::TOTAL_E      , { "toe" , "Tot. E. (eV/part)" , "% .10e" } } ,
        { ThermodynamicLogConfig::KINETIC_E    , { "kie" , "Kin. E. (eV/part)" , "% .10e" } } ,
        { ThermodynamicLogConfig::POTENTIAL_E  , { "poe" , "Pot. E. (eV/part)" , "% .10e" } } ,
        { ThermodynamicLogConfig::ELECTRON_E   , { "ele" , "Elec. E. (eV)"     , "% .10e" } } ,
        { ThermodynamicLogConfig::ION_TRANSFER_E, { "ite" , "Ion Transf. E. (eV)", "% .10e" } } ,
        { ThermodynamicLogConfig::TEMPERATURE  , { "tmp" , "Temp. (K)"        , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Tx           , { "tmx" , "   Tx (K)"        , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Ty           , { "tmy" , "   Ty (K)"        , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Tz           , { "tmz" , "   Tz (K)"        , "% 20.12f" } } ,
        { ThermodynamicLogConfig::PRESSURE     , { "pre" , "Press. (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Pxx          , { "pxx" , "   Pxx (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Pyy          , { "pyy" , "   Pyy (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Pzz          , { "pzz" , "   Pzz (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Pxy          , { "pxy" , "   Pxy (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Pxz          , { "pxz" , "   Pxz (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Pyz          , { "pyz" , "   Pyz (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Vxx          , { "vxx" , "   Vxx (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Vyy          , { "vyy" , "   Vyy (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Vzz          , { "vzz" , "   Vzz (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Vxy          , { "vxy" , "   Vxy (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Vxz          , { "vxz" , "   Vxz (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::Vyz          , { "vyz" , "   Vyz (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::SMISES       , { "smi" , "sMises (Pa)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::VOLUME       , { "vol" , "Vol. (ang^3)"     , "% 24.12f" } } ,
        { ThermodynamicLogConfig::MASS         , { "mas" , "Mass"             , "% 20.12f" } } ,
        { ThermodynamicLogConfig::BOX_A        , { "bxa" , "A (ang)"          , "% 20.12f" } } ,
        { ThermodynamicLogConfig::BOX_B        , { "bxb" , "B (ang)"          , "% 20.12f" } } ,
        { ThermodynamicLogConfig::BOX_C        , { "bxc" , "C (ang)"          , "% 20.12f" } } ,
        { ThermodynamicLogConfig::BOX_ALPHA    , { "baa" , "alpha (deg)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::BOX_BETA     , { "bab" , "beta (deg)"       , "% 20.12f" } } ,
        { ThermodynamicLogConfig::BOX_GAMMA    , { "bag" , "gamma (deg)"      , "% 20.12f" } } ,
        { ThermodynamicLogConfig::DENSITY      , { "rho" , "Rho (g/cm^3)"     , "% 20.12f" } }
      }
      ,
      {
        ThermodynamicLogConfig::TIME_STEP,
        ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::NB_PARTICLES,
        ThermodynamicLogConfig::TOTAL_E,
        ThermodynamicLogConfig::KINETIC_E,
        ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE,
        ThermodynamicLogConfig::Pxx,
        ThermodynamicLogConfig::Pyy,
        ThermodynamicLogConfig::Pzz,
        ThermodynamicLogConfig::Pxy,
        ThermodynamicLogConfig::Pxz,
        ThermodynamicLogConfig::Pyz,
        ThermodynamicLogConfig::BOX_A,
        ThermodynamicLogConfig::BOX_B,
        ThermodynamicLogConfig::BOX_C,
        ThermodynamicLogConfig::BOX_ALPHA,
        ThermodynamicLogConfig::BOX_BETA,
        ThermodynamicLogConfig::BOX_GAMMA,
        ThermodynamicLogConfig::VOLUME,
        ThermodynamicLogConfig::DENSITY
      }
    };

  inline std::vector<std::string> thermodynamic_log_split( const std::string& s )
  {
    std::vector<std::string> tokens;
    std::string::size_type pos = 0;
    while( pos != std::string::npos )
    {
      std::string::size_type next = s.find(';',pos);
      std::string token;
      if( next != std::string::npos ) { token = s.substr( pos , next-pos ); pos=next+1; }
      else { token = s.substr( pos ); pos = next; }
      tokens.push_back( token );
    }
    return tokens;
  }

  // Sets cfg.m_active_items from a named preset, or from a ';'-separated list of short keywords
  // (e.g. "stp;pht;mas;vol") looked up against cfg.m_avail_items. caller_name only decorates the
  // log/abort messages. Aborts with the full keyword table on an unrecognized keyword.
  inline void thermodynamic_log_apply_mode( ThermodynamicLogConfig& cfg, const std::string& mode, const char* caller_name )
  {
    if( mode == "default" || mode == "thermo_basic" ) {
      cfg.m_active_items = {
        ThermodynamicLogConfig::TIME_STEP, ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::TOTAL_E, ThermodynamicLogConfig::KINETIC_E, ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE, ThermodynamicLogConfig::PRESSURE, ThermodynamicLogConfig::SIM_STATUS };
    }
    else if( mode == "thermo" || mode == "thermo_full" ) {
      cfg.m_active_items = {
        ThermodynamicLogConfig::TIME_STEP, ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::TOTAL_E, ThermodynamicLogConfig::KINETIC_E, ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE, ThermodynamicLogConfig::Tx, ThermodynamicLogConfig::Ty, ThermodynamicLogConfig::Tz,
        ThermodynamicLogConfig::PRESSURE, ThermodynamicLogConfig::Pxx, ThermodynamicLogConfig::Pyy, ThermodynamicLogConfig::Pzz,
        ThermodynamicLogConfig::SIM_STATUS };
    }
    else if( mode == "vol_fluct_ortho_basic" ) {
      cfg.m_active_items = {
        ThermodynamicLogConfig::TIME_STEP, ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::TOTAL_E, ThermodynamicLogConfig::KINETIC_E, ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE, ThermodynamicLogConfig::PRESSURE, ThermodynamicLogConfig::VOLUME,
        ThermodynamicLogConfig::BOX_A, ThermodynamicLogConfig::BOX_B, ThermodynamicLogConfig::BOX_C,
        ThermodynamicLogConfig::DENSITY, ThermodynamicLogConfig::SIM_STATUS };
    }
    else if( mode == "vol_fluct_ortho" || mode == "vol_fluct_ortho_full" ) {
      cfg.m_active_items = {
        ThermodynamicLogConfig::TIME_STEP, ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::TOTAL_E, ThermodynamicLogConfig::KINETIC_E, ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE, ThermodynamicLogConfig::Tx, ThermodynamicLogConfig::Ty, ThermodynamicLogConfig::Tz,
        ThermodynamicLogConfig::PRESSURE, ThermodynamicLogConfig::Pxx, ThermodynamicLogConfig::Pyy, ThermodynamicLogConfig::Pzz,
        ThermodynamicLogConfig::VOLUME, ThermodynamicLogConfig::BOX_A, ThermodynamicLogConfig::BOX_B, ThermodynamicLogConfig::BOX_C,
        ThermodynamicLogConfig::DENSITY, ThermodynamicLogConfig::SIM_STATUS };
    }
    else if( mode == "vol_fluct_tricl_basic" ) {
      cfg.m_active_items = {
        ThermodynamicLogConfig::TIME_STEP, ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::TOTAL_E, ThermodynamicLogConfig::KINETIC_E, ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE, ThermodynamicLogConfig::PRESSURE, ThermodynamicLogConfig::VOLUME,
        ThermodynamicLogConfig::BOX_A, ThermodynamicLogConfig::BOX_B, ThermodynamicLogConfig::BOX_C,
        ThermodynamicLogConfig::BOX_ALPHA, ThermodynamicLogConfig::BOX_BETA, ThermodynamicLogConfig::BOX_GAMMA,
        ThermodynamicLogConfig::DENSITY, ThermodynamicLogConfig::SIM_STATUS };
    }
    else if( mode == "vol_fluct_tricl" || mode == "vol_fluct_tricl_full" ) {
      cfg.m_active_items = {
        ThermodynamicLogConfig::TIME_STEP, ThermodynamicLogConfig::PHYS_TIME,
        ThermodynamicLogConfig::TOTAL_E, ThermodynamicLogConfig::KINETIC_E, ThermodynamicLogConfig::POTENTIAL_E,
        ThermodynamicLogConfig::TEMPERATURE, ThermodynamicLogConfig::Tx, ThermodynamicLogConfig::Ty, ThermodynamicLogConfig::Tz,
        ThermodynamicLogConfig::PRESSURE, ThermodynamicLogConfig::Pxx, ThermodynamicLogConfig::Pyy, ThermodynamicLogConfig::Pzz,
        ThermodynamicLogConfig::Pxy, ThermodynamicLogConfig::Pxz, ThermodynamicLogConfig::Pyz,
        ThermodynamicLogConfig::VOLUME, ThermodynamicLogConfig::BOX_A, ThermodynamicLogConfig::BOX_B, ThermodynamicLogConfig::BOX_C,
        ThermodynamicLogConfig::BOX_ALPHA, ThermodynamicLogConfig::BOX_BETA, ThermodynamicLogConfig::BOX_GAMMA,
        ThermodynamicLogConfig::DENSITY, ThermodynamicLogConfig::SIM_STATUS };
    }
    else if( mode == "mechanical" ) {
      cfg.m_active_items = {
        ThermodynamicLogConfig::TIME_STEP, ThermodynamicLogConfig::PHYS_TIME, ThermodynamicLogConfig::NB_PARTICLES,
        ThermodynamicLogConfig::SIM_STATUS, ThermodynamicLogConfig::TOTAL_E, ThermodynamicLogConfig::KINETIC_E,
        ThermodynamicLogConfig::POTENTIAL_E, ThermodynamicLogConfig::TEMPERATURE, ThermodynamicLogConfig::PRESSURE,
        ThermodynamicLogConfig::SMISES, ThermodynamicLogConfig::VOLUME, ThermodynamicLogConfig::MASS };
    }
    else if( mode == "dump_default" ) {
      cfg.m_active_items = thermodynamic_dump_config_default.m_active_items;
    }
    else
    {
      cfg.m_active_items.clear();
      for( const std::string& token : thermodynamic_log_split(mode) )
      {
        if( token.empty() ) { continue; }
        bool found = false;
        for( const auto& item : cfg.m_avail_items )
        {
          if( item.second.m_name == token )
          {
            cfg.m_active_items.push_back( item.first );
            found = true;
            break;
          }
        }
        if( ! found )
        {
          lerr << caller_name << ": unrecognized log_mode keyword '"<<token<<"' in '"<<mode<<"'"<<std::endl;
          lerr << "Available keywords:" << std::endl;
          for( const auto& item : cfg.m_avail_items )
          {
            lerr << "  " << item.second.m_name << " : " << item.second.m_header << std::endl;
          }
          std::abort();
        }
      }
    }
  }

  // Overrides cfg.m_avail_items[...].m_format positionally for the currently active items, from a
  // ';'-separated list of printf-style format strings (e.g. "%9.3f;%12.6e"). Fewer entries than
  // m_active_items: only the first columns are overridden. Empty tokens leave that column's
  // current default format untouched. Must be called AFTER thermodynamic_log_apply_mode(), and
  // before the first print_log() call (which is when column widths get computed from m_format).
  inline void thermodynamic_log_apply_format( ThermodynamicLogConfig& cfg, const std::string& format_list )
  {
    auto tokens = thermodynamic_log_split(format_list);
    for( size_t k=0; k<cfg.m_active_items.size() && k<tokens.size(); k++ )
    {
      if( ! tokens[k].empty() )
      {
        cfg.m_avail_items[ cfg.m_active_items[k] ].m_format = tokens[k];
        cfg.m_avail_items[ cfg.m_active_items[k] ].len = -1;
      }
    }
  }

  inline void thermodynamic_log_conversion_factors( bool internal_units, double& conv_temperature, double& conv_energy, double& conv_pressure, double& conv_density )
  {
    conv_temperature = 1.e4 * onika::physics::atomicMass / onika::physics::boltzmann;       // internal units to Kelvin
    conv_energy      = 1.e4 * onika::physics::atomicMass / onika::physics::elementaryCharge; // internal units to eV
    conv_pressure    = 1.e4 * onika::physics::atomicMass * 1e30;                             // internal units to Pascal
    conv_density     = onika::physics::atomicMass * 1e3 * 1e24;                              // internal units to g/cm^3
    if( internal_units )
    {
      conv_temperature = 1.0;
      conv_energy = 1.0;
      conv_pressure = 1.0;
      conv_density = 1.0;
    }
  }

  // Fills the full LOG_ITEM_COUNT values array from a thermodynamic state snapshot -- shared by
  // print_thermodynamic_state and dump_thermodynamic_state so both can honor an arbitrary log_mode
  // column selection over the exact same set of quantities. el_energy/ion_energy are already
  // totals (not per-particle), matching total_electronic_energy/ion_transfer_energy's own units.
  inline void thermodynamic_log_fill_values(
    double values[ThermodynamicLogConfig::LOG_ITEM_COUNT],
    const ThermodynamicState& sim_info, const Domain& domain,
    long timestep, double physical_time, double el_energy, double ion_energy, double lb_inbalance,
    double conv_temperature, double conv_energy, double conv_pressure, double conv_density )
  {
    double total_energy_int_unit = sim_info.total_energy() + el_energy;

    Mat3d xform = domain.xform();
    Vec3d a = xform * Vec3d{domain.extent().x - domain.origin().x,0.,0.};
    Vec3d b = xform * Vec3d{0.,domain.extent().y - domain.origin().y,0.};
    Vec3d c = xform * Vec3d{0.,0.,domain.extent().z - domain.origin().z};
    double A = norm(a);
    double B = norm(b);
    double C = norm(c);
    double ALPHA = acos(dot(b,c)/(B*C))/acos(-1.)*180.;
    double BETA  = acos(dot(c,a)/(B*C))/acos(-1.)*180.;
    double GAMMA = acos(dot(a,b)/(B*C))/acos(-1.)*180.;

    values[ThermodynamicLogConfig::TIME_STEP]      = timestep;
    values[ThermodynamicLogConfig::PHYS_TIME]      = physical_time;
    values[ThermodynamicLogConfig::NB_PARTICLES]   = sim_info.particle_count();
    values[ThermodynamicLogConfig::SIM_STATUS]     = lb_inbalance;
    values[ThermodynamicLogConfig::TOTAL_E]        = total_energy_int_unit          / sim_info.particle_count() * conv_energy;
    values[ThermodynamicLogConfig::KINETIC_E]      = sim_info.kinetic_energy_scal() / sim_info.particle_count() * conv_energy;
    values[ThermodynamicLogConfig::POTENTIAL_E]    = sim_info.potential_energy()    / sim_info.particle_count() * conv_energy;
    values[ThermodynamicLogConfig::ELECTRON_E]     = el_energy  * conv_energy; // total, not per-particle
    values[ThermodynamicLogConfig::ION_TRANSFER_E] = ion_energy * conv_energy; // total, not per-particle
    values[ThermodynamicLogConfig::TEMPERATURE]    = sim_info.temperature_scal() / sim_info.particle_count() * conv_temperature;
    values[ThermodynamicLogConfig::Tx]             = sim_info.temperature().x   / sim_info.particle_count() * conv_temperature;
    values[ThermodynamicLogConfig::Ty]             = sim_info.temperature().y   / sim_info.particle_count() * conv_temperature;
    values[ThermodynamicLogConfig::Tz]             = sim_info.temperature().z   / sim_info.particle_count() * conv_temperature;
    values[ThermodynamicLogConfig::PRESSURE]       = sim_info.pressure_scal() * conv_pressure;
    // Pxx/Pyy/Pzz/Pxy/Pxz/Pyz are "the stress tensor" -- always include the kinetic contribution
    // (full_stress_tensor() = stress_tensor()+kinetic_tensor() is the true physical stress).
    // Vxx/Vyy/Vzz/Vxy/Vxz/Vyz are "the virial tensor" -- pure virial, no kinetic term. Both are
    // independently selectable via log_mode (dump_thermodynamic_state's legacy is_dump_virial
    // 9-component S11..S33 dump is a separate, older mechanism kept for backward compatibility).
    values[ThermodynamicLogConfig::Pxx]            = sim_info.full_stress_tensor().m11 * conv_pressure;
    values[ThermodynamicLogConfig::Pyy]            = sim_info.full_stress_tensor().m22 * conv_pressure;
    values[ThermodynamicLogConfig::Pzz]            = sim_info.full_stress_tensor().m33 * conv_pressure;
    values[ThermodynamicLogConfig::Pxy]            = sim_info.full_stress_tensor().m12 * conv_pressure;
    values[ThermodynamicLogConfig::Pxz]            = sim_info.full_stress_tensor().m13 * conv_pressure;
    values[ThermodynamicLogConfig::Pyz]            = sim_info.full_stress_tensor().m23 * conv_pressure;
    values[ThermodynamicLogConfig::Vxx]            = sim_info.stress_tensor().m11 * conv_pressure;
    values[ThermodynamicLogConfig::Vyy]            = sim_info.stress_tensor().m22 * conv_pressure;
    values[ThermodynamicLogConfig::Vzz]            = sim_info.stress_tensor().m33 * conv_pressure;
    values[ThermodynamicLogConfig::Vxy]            = sim_info.stress_tensor().m12 * conv_pressure;
    values[ThermodynamicLogConfig::Vxz]            = sim_info.stress_tensor().m13 * conv_pressure;
    values[ThermodynamicLogConfig::Vyz]            = sim_info.stress_tensor().m23 * conv_pressure;
    values[ThermodynamicLogConfig::SMISES]         = sim_info.vonmises_scal() * conv_pressure;
    values[ThermodynamicLogConfig::VOLUME]         = sim_info.volume();
    values[ThermodynamicLogConfig::MASS]           = sim_info.mass();
    values[ThermodynamicLogConfig::BOX_A]          = A;
    values[ThermodynamicLogConfig::BOX_B]          = B;
    values[ThermodynamicLogConfig::BOX_C]          = C;
    values[ThermodynamicLogConfig::BOX_ALPHA]      = ALPHA;
    values[ThermodynamicLogConfig::BOX_BETA]       = BETA;
    values[ThermodynamicLogConfig::BOX_GAMMA]      = GAMMA;
    values[ThermodynamicLogConfig::DENSITY]        = sim_info.mass()/sim_info.volume() * conv_density;
  }

}
