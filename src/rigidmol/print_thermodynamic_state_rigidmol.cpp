#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/string_utils.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exanb/core/physics_constants.h>

namespace exaStamp
{
  using namespace exanb;

  struct PrintThermodynamicStateRigidmolNode : public OperatorNode
  {  
    // thermodynamic state & physics data
    ADD_SLOT( long               , timestep            , INPUT , REQUIRED );
    ADD_SLOT( double             , physical_time       , INPUT , REQUIRED );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT );

    // printing options
    ADD_SLOT( bool               , print_header        , INPUT, false );
    ADD_SLOT( bool               , internal_units      , INPUT, false );

    // LB and particle movement statistics
    ADD_SLOT( long               , lb_counter          , INPUT_OUTPUT );
    ADD_SLOT( long               , move_counter        , INPUT_OUTPUT );
    ADD_SLOT( long               , domain_ext_counter  , INPUT_OUTPUT );
    ADD_SLOT( double             , lb_inbalance_max    , INPUT_OUTPUT );

    // optional physics quantities
    ADD_SLOT( double             , electronic_energy   , INPUT, OPTIONAL );

    inline void execute () override final
    {
      double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann ;
      double conv_energy = 1.e4 * legacy_constant::atomicMass / legacy_constant::elementaryCharge;
      double conv_pressure = legacy_constant::atomicMass * 1e20;
      
      static const std::string header = "     Step     Time (ps)     Particles   Mv/Ext/Imb.  Tot. E. (eV/part)  Kin. E. (eV/part)  Rot. E. (eV/part)  Pot. E. (eV/part)  Temperature   Pressure     sMises     Volume       Mass";
      
      
      if( *internal_units )
      {
        conv_energy = 1.0;
        conv_temperature = 1.0;
        conv_pressure = 1.0;
      }

      bool lb_flag = (*lb_counter) > 0 ;
      long move_count = *move_counter ;
      long domext_count = *domain_ext_counter;
      double lb_inbalance = *lb_inbalance_max;
      
      //std::cout << "lb_counter = "<< *lb_counter << std::endl;
      
      *lb_counter = 0;
      *move_counter = 0;
      *domain_ext_counter = 0;
      *lb_inbalance_max = 0.0;
      
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);

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
          lb_value = "  N/A  ";
        }
        else
        {
          lb_value = format_string("%.1e", lb_inbalance);
        }
      }

      double total_energy_int_unit = sim_info.total_energy_rigidmol();
      if( electronic_energy.has_value() )
      {
        total_energy_int_unit += *electronic_energy;
      }

      if( *print_header )
      {
        lout << header;
        if( electronic_energy.has_value() ) { lout << "  Elect. Energy"; }
        lout << std::endl;
      }

//      lout<<format_string("%9ld % .6e %13ld  %c %c %8s  % .10e  % .10e  % .10e  % 11.3f % .3e % .3e % .3e % .3e",
      lout<<format_string("%9ld % .6e %13ld  %c %c %8s  % .10e  %  .10e  % .10e  % .10e  % 11.3f % .3e % .3e % .3e % .3e",
        *timestep,
	      *physical_time,
        sim_info.particle_count(),
        lb_move_char,domext_char,lb_value,
        total_energy_int_unit * conv_energy / sim_info.particle_count(),
        sim_info.kinetic_energy_scal() * conv_energy / sim_info.particle_count(),
        sim_info.rotational_energy_scal() * conv_energy / sim_info.particle_count(),
        sim_info.potential_energy() * conv_energy / sim_info.particle_count(),
        sim_info.temperature_rigidmol_scal() * conv_temperature,
        sim_info.pressure_scal() * conv_pressure,
        sim_info.vonmises_scal() * conv_pressure,
        sim_info.volume(),
        sim_info.mass()) ;
      
      if( electronic_energy.has_value() )
      {
        lout << format_string(" % .7e",(*electronic_energy) * conv_energy / sim_info.particle_count() );
      }
      lout << std::endl;
    }

  };
    
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "print_thermodynamic_state_rigidmol", make_simple_operator<PrintThermodynamicStateRigidmolNode> );
  }

}

