#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/string_utils.h>
#include <exanb/core/print_utils.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exanb/core/physics_constants.h>

#include <sstream>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct DumpThermodynamicStateRigidmolNode : public OperatorNode
  {  
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( long               , timestep            , INPUT, REQUIRED);
    ADD_SLOT( double             , physical_time       , INPUT );
    ADD_SLOT( bool               , print_header        , INPUT, true );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT, REQUIRED);
    ADD_SLOT( double             , electronic_energy   , INPUT, OPTIONAL );
    ADD_SLOT( std::string        , file                , INPUT , "thermodynamic_state.csv" );
    ADD_SLOT( bool               , force_flush_file    , INPUT , false );
    ADD_SLOT( bool               , force_append_thermo , INPUT , false );
    
    inline void execute () override final
    {
      static const double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann ;
      static const double conv_energy = 1.e4 * legacy_constant::atomicMass / legacy_constant::elementaryCharge;
      static const double conv_pressure = legacy_constant::atomicMass * 1e20;
      static const std::string header = "     Step     Time (ps)     Particles   Tot. E. (eV/part)  Kin. E. (eV/part)  Rot. E. (eV/part)  Pot. E. (eV/part)  Temperature   Pressure     sMises     Volume       Mass  Kin. Tx  Kin. Ty  Kin. Tz  Rot. Tx  Rot. Ty  Rot. Tz ";

      // MPI Initialization
      int rank = 0;
      MPI_Comm_rank(*mpi, &rank);

      // initialisation : remove output.csv
      if(rank!=0) { return; }

      const ThermodynamicState& sim_info = *(this->thermodynamic_state);

      std::ostringstream oss;

      if( *print_header )
      {
        oss << header;
        if( electronic_energy.has_value() ) { oss << "  Elect. Energy"; }
        oss << '\n';
      }

      double total_energy_int_unit = sim_info.total_energy_rigidmol();
      if( electronic_energy.has_value() )
      {
        total_energy_int_unit += *electronic_energy;
      }

      oss << format_string("%9ld % .6e %13ld  % .10e  % .10e  % .10e  % .10e  % 11.3f % .3e % .3e % .3e % .3e % 11.3f % 11.3f % 11.3f % 11.3f % 11.3f % 11.3f /",
                             *timestep,
                             *physical_time,
                             sim_info.particle_count(),
                             total_energy_int_unit * conv_energy / sim_info.particle_count(),
                             sim_info.kinetic_energy_scal() * conv_energy / sim_info.particle_count(),
                             sim_info.rotational_energy_scal() * conv_energy / sim_info.particle_count(),
                             sim_info.potential_energy() * conv_energy / sim_info.particle_count(),
                             sim_info.temperature_rigidmol_scal() * conv_temperature,
                             sim_info.pressure_scal() * conv_pressure,
                             sim_info.vonmises_scal() * conv_pressure,
                             sim_info.volume(),
                             sim_info.mass(),
                             sim_info.kinetic_temperature_x() * conv_temperature,
                             sim_info.kinetic_temperature_y() * conv_temperature,
                             sim_info.kinetic_temperature_z() * conv_temperature,
                             sim_info.rotational_temperature_x() * conv_temperature,
                             sim_info.rotational_temperature_y() * conv_temperature,
                             sim_info.rotational_temperature_z() * conv_temperature);
                             
      if( electronic_energy.has_value() )
      {
        oss << format_string(" % .7e",(*electronic_energy) * conv_energy / sim_info.particle_count() );
      }
      
      oss << "\n";

      FileAppendWriteBuffer::instance().append_to_file( *file , oss.str(), *force_append_thermo );

      if( *force_flush_file )
      {
        FileAppendWriteBuffer::instance().flush();
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
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "dump_thermodynamic_state_rigidmol", make_simple_operator<DumpThermodynamicStateRigidmolNode> );
  }

}

