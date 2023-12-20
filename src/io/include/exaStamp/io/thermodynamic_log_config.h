#include <algorithm>
#include <string>
#include <vector>
#include <exanb/core/string_utils.h>

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
      TEMPERATURE ,
      Tx ,
      Ty ,
      Tz ,
      PRESSURE ,
      Px ,
      Py ,
      Pz ,
      Pxy ,
      Pxz ,
      Pyz ,
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
    
    template<class StreamT>
    inline void print_log( StreamT & out , const double values[LOG_ITEM_COUNT] , bool header, long move_count , long domext_count , double lb_inbalance, bool lb_flag )
    {
      for(auto & p : m_avail_items)
      {
        if(p.second.len==-1)
        {
          std::string s = format_string(p.second.m_format,0.0);
          size_t pfxlen = ( p.first == SIM_STATUS ) ? 3 : 0;
          p.second.len = std::max( s.length() + pfxlen ,  p.second.m_header.length() );
        }
      }
      if( header )
      {
        for(auto id : m_active_items) out << m_sep << format_string("%*s", m_avail_items[id].len, m_avail_items[id].m_header );
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
              lb_value = format_string(m_avail_items[SIM_STATUS].m_format, values[SIM_STATUS] );
            }
          }
          
          out << m_sep << format_string("%c %c%*s",lb_move_char,domext_char,m_avail_items[SIM_STATUS].len-3,lb_value);
        }
        else
        {
          out << m_sep << format_string(m_avail_items[id].m_format,values[id]);
        }
      }
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
        { ThermodynamicLogConfig::ELECTRON_E   , { "ele" , "Elec. E. (eV/part)", "% .10e" } } ,
        { ThermodynamicLogConfig::TEMPERATURE  , { "tmp" , "Temp. (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::Tx           , { "tmx" , "   Tx (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::Ty           , { "tmy" , "   Ty (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::Tz           , { "tmz" , "   Tz (K)"        , "% 10.3f" } } ,
        { ThermodynamicLogConfig::PRESSURE     , { "pre" , "Press. (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Px           , { "prx" , "    Px (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Py           , { "pry" , "    Py (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pz           , { "prz" , "    Pz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pxy          , { "pxy" , "   Pxy (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pxz          , { "pxz" , "   Pxz (Pa)"      , "% 11.3e" } } ,
        { ThermodynamicLogConfig::Pyz          , { "pyz" , "   Pyz (Pa)"      , "% 11.3e" } } ,
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

}

