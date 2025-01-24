#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/string_utils.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exanb/core/physics_constants.h>
#include <exanb/core/domain.h>

#include <exaStamp/io/thermodynamic_log_config.h>

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

    // LB and particle movement statistics
    ADD_SLOT( long               , lb_counter          , INPUT_OUTPUT );
    ADD_SLOT( long               , move_counter        , INPUT_OUTPUT );
    ADD_SLOT( long               , domain_ext_counter  , INPUT_OUTPUT );
    ADD_SLOT( double             , lb_inbalance_max    , INPUT_OUTPUT );

    // optional physics quantities
    ADD_SLOT( double             , electronic_energy   , INPUT, OPTIONAL );

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
        if( *log_mode == "default" || *log_mode == "thermo_basic") {
          log_config->m_active_items = {
            ThermodynamicLogConfig::TIME_STEP,
            ThermodynamicLogConfig::PHYS_TIME,
            ThermodynamicLogConfig::TOTAL_E,
            ThermodynamicLogConfig::KINETIC_E,
            ThermodynamicLogConfig::POTENTIAL_E,
            ThermodynamicLogConfig::TEMPERATURE,
            ThermodynamicLogConfig::PRESSURE,
            ThermodynamicLogConfig::SIM_STATUS };
	}
        else if( *log_mode == "thermo" || *log_mode == "thermo_full" )
        {
          log_config->m_active_items = {
            ThermodynamicLogConfig::TIME_STEP,
            ThermodynamicLogConfig::PHYS_TIME,
            ThermodynamicLogConfig::TOTAL_E,
            ThermodynamicLogConfig::KINETIC_E,
            ThermodynamicLogConfig::POTENTIAL_E,
            ThermodynamicLogConfig::TEMPERATURE,
            ThermodynamicLogConfig::Tx,
            ThermodynamicLogConfig::Ty,
            ThermodynamicLogConfig::Tz,
            ThermodynamicLogConfig::PRESSURE,
            ThermodynamicLogConfig::Px,
            ThermodynamicLogConfig::Py,
            ThermodynamicLogConfig::Pz,
            ThermodynamicLogConfig::SIM_STATUS };
	}
        else if( *log_mode == "vol_fluct_ortho_basic" )
        {
          log_config->m_active_items = {
            ThermodynamicLogConfig::TIME_STEP,
            ThermodynamicLogConfig::PHYS_TIME,
            ThermodynamicLogConfig::TOTAL_E,
            ThermodynamicLogConfig::KINETIC_E,
            ThermodynamicLogConfig::POTENTIAL_E,
            ThermodynamicLogConfig::TEMPERATURE,
            ThermodynamicLogConfig::PRESSURE,
            ThermodynamicLogConfig::VOLUME,
            ThermodynamicLogConfig::BOX_A,
            ThermodynamicLogConfig::BOX_B,
            ThermodynamicLogConfig::BOX_C,
            ThermodynamicLogConfig::DENSITY,
            ThermodynamicLogConfig::SIM_STATUS };
        }
        else if( *log_mode == "vol_fluct_ortho" || *log_mode == "vol_fluct_ortho_full" )
        {
          log_config->m_active_items = {
            ThermodynamicLogConfig::TIME_STEP,
            ThermodynamicLogConfig::PHYS_TIME,
            ThermodynamicLogConfig::TOTAL_E,
            ThermodynamicLogConfig::KINETIC_E,
            ThermodynamicLogConfig::POTENTIAL_E,
            ThermodynamicLogConfig::TEMPERATURE,
            ThermodynamicLogConfig::Tx,
            ThermodynamicLogConfig::Ty,
            ThermodynamicLogConfig::Tz,
            ThermodynamicLogConfig::PRESSURE,
            ThermodynamicLogConfig::Px,
            ThermodynamicLogConfig::Py,
            ThermodynamicLogConfig::Pz,
            ThermodynamicLogConfig::VOLUME,
            ThermodynamicLogConfig::BOX_A,
            ThermodynamicLogConfig::BOX_B,
            ThermodynamicLogConfig::BOX_C,
            ThermodynamicLogConfig::DENSITY,
            ThermodynamicLogConfig::SIM_STATUS };
        }
        else if( *log_mode == "vol_fluct_tricl_basic")
        {
          log_config->m_active_items = {
            ThermodynamicLogConfig::TIME_STEP,
            ThermodynamicLogConfig::PHYS_TIME,
            ThermodynamicLogConfig::TOTAL_E,
            ThermodynamicLogConfig::KINETIC_E,
            ThermodynamicLogConfig::POTENTIAL_E,
            ThermodynamicLogConfig::TEMPERATURE,
            ThermodynamicLogConfig::PRESSURE,
            ThermodynamicLogConfig::VOLUME,
            ThermodynamicLogConfig::BOX_A,
            ThermodynamicLogConfig::BOX_B,
            ThermodynamicLogConfig::BOX_C,
            ThermodynamicLogConfig::BOX_ALPHA,
            ThermodynamicLogConfig::BOX_BETA,
            ThermodynamicLogConfig::BOX_GAMMA,
            ThermodynamicLogConfig::DENSITY,
            ThermodynamicLogConfig::SIM_STATUS };
        }
        else if( *log_mode == "vol_fluct_tricl" || *log_mode == "vol_fluct_tricl_full")
        {
          log_config->m_active_items = {
            ThermodynamicLogConfig::TIME_STEP,
            ThermodynamicLogConfig::PHYS_TIME,
            ThermodynamicLogConfig::TOTAL_E,
            ThermodynamicLogConfig::KINETIC_E,
            ThermodynamicLogConfig::POTENTIAL_E,
            ThermodynamicLogConfig::TEMPERATURE,
            ThermodynamicLogConfig::Tx,
            ThermodynamicLogConfig::Ty,
            ThermodynamicLogConfig::Tz,
            ThermodynamicLogConfig::PRESSURE,
            ThermodynamicLogConfig::Px,
            ThermodynamicLogConfig::Py,
            ThermodynamicLogConfig::Pz,
            ThermodynamicLogConfig::Pxy,
            ThermodynamicLogConfig::Pxz,
            ThermodynamicLogConfig::Pyz,
            ThermodynamicLogConfig::VOLUME,
            ThermodynamicLogConfig::BOX_A,
            ThermodynamicLogConfig::BOX_B,
            ThermodynamicLogConfig::BOX_C,
            ThermodynamicLogConfig::BOX_ALPHA,
            ThermodynamicLogConfig::BOX_BETA,
            ThermodynamicLogConfig::BOX_GAMMA,
            ThermodynamicLogConfig::DENSITY,
            ThermodynamicLogConfig::SIM_STATUS };
        }
        else if( *log_mode == "mechanical" )
        {
          log_config->m_active_items = {
            ThermodynamicLogConfig::TIME_STEP,
            ThermodynamicLogConfig::PHYS_TIME,
            ThermodynamicLogConfig::NB_PARTICLES,
            ThermodynamicLogConfig::SIM_STATUS,
            ThermodynamicLogConfig::TOTAL_E,
            ThermodynamicLogConfig::KINETIC_E,
            ThermodynamicLogConfig::POTENTIAL_E,
            ThermodynamicLogConfig::TEMPERATURE,
            ThermodynamicLogConfig::PRESSURE,
            ThermodynamicLogConfig::SMISES,
            ThermodynamicLogConfig::VOLUME,
            ThermodynamicLogConfig::MASS };
        }
        else
        {
          // parse log_mode as a list of keywords
          // format like stp:pht:mas:vol ...
          std::string::size_type pos = 0;
          while( pos != std::string::npos )
          {
            std::string::size_type next = log_mode->find(';',pos);
            std::string token;
            if( next != std::string::npos ) { token = log_mode->substr( pos , next-pos ); pos=next+1; }
            else { token = log_mode->substr( pos ); pos = next; }
          }
        }
      }

      double el_energy = 0.0;
      if( electronic_energy.has_value() )
      {
        el_energy = *electronic_energy;
        log_config->m_active_items.push_back( ThermodynamicLogConfig::ELECTRON_E );
      }
    
      double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann ;       // internal units to Kelvin
      //double conv_energy = 1.e4 * legacy_constant::atomicMass / legacy_constant::elementaryCharge;    // internal units to Joule
      double conv_energy = 1.e4 * legacy_constant::atomicMass / legacy_constant::elementaryCharge;      // internal units to eV
//      double conv_pressure = legacy_constant::atomicMass * 1e20;					// ORIGINAL LINE - NO IDEA WHAT UNITS THIS IS (1e-14 Pascal)
      double conv_pressure = 1.e4 * legacy_constant::atomicMass * 1e30;                                 // internal units to Pascal
      double conv_density = legacy_constant::atomicMass*1e3*1e24;                                       // internal units to g/cm^3

            
      if( *internal_units )
      {
        conv_energy = 1.0;
        conv_temperature = 1.0;
        conv_pressure = 1.0;
        conv_density = 1.0;
      }

      bool lb_flag = (*lb_counter) > 0 ;
      long move_count = *move_counter ;
      long domext_count = *domain_ext_counter;
      double lb_inbalance = *lb_inbalance_max;
      
      *lb_counter = 0;
      *move_counter = 0;
      *domain_ext_counter = 0;
      *lb_inbalance_max = 0.0;
      
      const ThermodynamicState& sim_info = *thermodynamic_state;
      double total_energy_int_unit = sim_info.total_energy() + el_energy;
      
      Mat3d xform = domain->xform();
      Vec3d a = xform * Vec3d{domain->extent().x - domain->origin().x,0.,0.} ;
      Vec3d b = xform * Vec3d{0.,domain->extent().y - domain->origin().y,0.} ;
      Vec3d c = xform * Vec3d{0.,0.,domain->extent().z - domain->origin().z} ;
      double A = norm(a) ;
      double B = norm(b) ;
      double C = norm(c) ;
      double ALPHA = acos(dot(b,c)/(B*C))/acos(-1.)*180. ;
      double BETA  = acos(dot(c,a)/(B*C))/acos(-1.)*180. ;
      double GAMMA = acos(dot(a,b)/(B*C))/acos(-1.)*180. ;

      double values[ThermodynamicLogConfig::LOG_ITEM_COUNT];
      values[ThermodynamicLogConfig::TIME_STEP]    = *timestep;
      values[ThermodynamicLogConfig::PHYS_TIME]    = *physical_time;
      values[ThermodynamicLogConfig::NB_PARTICLES] = sim_info.particle_count();
      values[ThermodynamicLogConfig::SIM_STATUS]   = lb_inbalance;
      values[ThermodynamicLogConfig::TOTAL_E]      = total_energy_int_unit          / sim_info.particle_count() * conv_energy;
      values[ThermodynamicLogConfig::KINETIC_E]    = sim_info.kinetic_energy_scal() / sim_info.particle_count() * conv_energy;
      values[ThermodynamicLogConfig::POTENTIAL_E]  = sim_info.potential_energy()    / sim_info.particle_count() * conv_energy;
      values[ThermodynamicLogConfig::ELECTRON_E]   = el_energy                      / sim_info.particle_count() * conv_energy;
      values[ThermodynamicLogConfig::TEMPERATURE]  = sim_info.temperature_scal()    / sim_info.particle_count() * conv_temperature;
      values[ThermodynamicLogConfig::Tx]           = sim_info.temperature().x       / sim_info.particle_count() * conv_temperature;
      values[ThermodynamicLogConfig::Ty]           = sim_info.temperature().y       / sim_info.particle_count() * conv_temperature;
      values[ThermodynamicLogConfig::Tz]           = sim_info.temperature().z       / sim_info.particle_count() * conv_temperature;
      values[ThermodynamicLogConfig::PRESSURE]     = sim_info.pressure_scal()                                   * conv_pressure;
      values[ThermodynamicLogConfig::Px]           = sim_info.stress_tensor().m11                          * conv_pressure;
      values[ThermodynamicLogConfig::Py]           = sim_info.stress_tensor().m22                          * conv_pressure;
      values[ThermodynamicLogConfig::Pz]           = sim_info.stress_tensor().m33                          * conv_pressure;
      values[ThermodynamicLogConfig::Pxy]          = sim_info.full_stress_tensor().m12                          * conv_pressure;
      values[ThermodynamicLogConfig::Pxz]          = sim_info.full_stress_tensor().m13                          * conv_pressure;                  
      values[ThermodynamicLogConfig::Pyz]          = sim_info.full_stress_tensor().m23                          * conv_pressure;                  
      // values[ThermodynamicLogConfig::Px]           = sim_info.pressure().x                                      * conv_pressure;
      // values[ThermodynamicLogConfig::Py]           = sim_info.pressure().y                                      * conv_pressure;
      // values[ThermodynamicLogConfig::Pz]           = sim_info.pressure().z                                      * conv_pressure;
      // values[ThermodynamicLogConfig::Pxy]          = sim_info.stress_tensor().m12                               * conv_pressure;
      // values[ThermodynamicLogConfig::Pxz]          = sim_info.stress_tensor().m13                               * conv_pressure;
      // values[ThermodynamicLogConfig::Pyz]          = sim_info.stress_tensor().m23                               * conv_pressure;
      values[ThermodynamicLogConfig::SMISES]       = sim_info.vonmises_scal()                                   * conv_pressure;
      values[ThermodynamicLogConfig::VOLUME]       = sim_info.volume();
      values[ThermodynamicLogConfig::MASS]         = sim_info.mass();
      values[ThermodynamicLogConfig::BOX_A]        = A;
      values[ThermodynamicLogConfig::BOX_B]        = B;
      values[ThermodynamicLogConfig::BOX_C]        = C;
      values[ThermodynamicLogConfig::BOX_ALPHA]    = ALPHA;
      values[ThermodynamicLogConfig::BOX_BETA]     = BETA;
      values[ThermodynamicLogConfig::BOX_GAMMA]    = GAMMA;
      values[ThermodynamicLogConfig::DENSITY]      = sim_info.mass()/sim_info.volume()                          *conv_density;

      log_config->print_log( lout , values, *print_header , move_count , domext_count , lb_inbalance , lb_flag );
    }

  };
    
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "print_thermodynamic_state", make_simple_operator<PrintThermodynamicStateNode> );
  }

}

