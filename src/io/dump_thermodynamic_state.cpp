#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/string_utils.h>
#include <exanb/core/print_utils.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exanb/core/physics_constants.h>
#include <exanb/core/domain.h>

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
    ADD_SLOT( bool               , print_header        , INPUT, true );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT, REQUIRED);
    ADD_SLOT( double             , electronic_energy   , INPUT, OPTIONAL );
    ADD_SLOT( std::string        , file                , INPUT , "thermodynamic_state.csv" );
    ADD_SLOT( bool               , force_flush_file    , INPUT , false );
    ADD_SLOT( bool               , force_append_thermo , INPUT , false );    
    ADD_SLOT( bool               , is_dump_virial      , INPUT , false);
    // NEW
    ADD_SLOT(Domain              , domain              , INPUT , OPTIONAL, DocString{"Deformation box matrix"} );


    inline void execute () override final
    {
      static const double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann ;	// internal units to Kelvin
      //static const double conv_energy = 1.e4 * legacy_constant::atomicMass;					// internal units to Joule
      static const double conv_energy = 1.e4 * legacy_constant::atomicMass / legacy_constant::elementaryCharge;	// internal units to eV
      //static const double conv_pressure = legacy_constant::atomicMass * 1e20;					// ORIGINAL LINE - NO IDEA WHAT UNITS THIS IS (1e-14 Pascal)
      static const double conv_pressure = 1.e4 * legacy_constant::atomicMass * 1e30;				// internal units to Pascal
      static const double conv_density = legacy_constant::atomicMass*1e3*1e24; 					// internal units to g/cm^3
	
//     static const std::string header = "###  Step     Time (ps)     Particles  Tot. E. (eV/part)  Kin. E. (eV/part)  Pot. E. (eV/part)  Temp. (K)                     Tx/Ty/Tz (K)    Press. (Pa)                                            Pxx/Pyy/Pzz (Pa) Pxy/Pxz/Pyz (Pa)   sMises (Pa)                            A/B/C (ang)    alpha/beta/gamma (deg)     Vol. (ang^3)  Rho (g/cm^3)";
     static const std::string header = "# Step     Time (ps)     Particles  Tot. E. (eV/part)  Kin. E. (eV/part)  Pot. E. (eV/part)  Temp. (K) Pxx Pyy Pzz Pxy Pxz Pyz (Pa) A/B/C (ang)    alpha/beta/gamma (deg)     Vol. (ang^3)  Rho (g/cm^3)";     

      bool is_dump_virial = *(this->is_dump_virial);      

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
        if( is_dump_virial ) { oss << "  S11  S12  S13  S21  S22  S23  S31  S32  S33"; }
        oss << '\n';
      }

      double total_energy_int_unit = sim_info.total_energy();
      if( electronic_energy.has_value() )
      {
        total_energy_int_unit += *electronic_energy;
      }
                             
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

      oss <<format_string("%9ld % .6e %13ld  % .10e  % .10e  % .10e  % 9.12f  % 9.12f  % 9.12f  % 9.12f  % 9.12f  % 9.12f % 9.12f % 12.12f % 12.12f % 12.12f % 9.12f % 9.12f % 9.12f % 16.12f % 13.12f ",
        *timestep,
        *physical_time,
        sim_info.particle_count(),
        total_energy_int_unit                / sim_info.particle_count() * conv_energy,
        sim_info.kinetic_energy_scal()       / sim_info.particle_count() * conv_energy,
        sim_info.potential_energy()          / sim_info.particle_count() * conv_energy,
        sim_info.temperature_scal()          / sim_info.particle_count() * conv_temperature,
	sim_info.full_stress_tensor().m11                                            * conv_pressure,
	sim_info.full_stress_tensor().m22                                            * conv_pressure,
	sim_info.full_stress_tensor().m33                                            * conv_pressure,
	// sim_info.pressure().x                                            * conv_pressure,
	// sim_info.pressure().y                                            * conv_pressure,
	// sim_info.pressure().z                                            * conv_pressure,
	sim_info.full_stress_tensor().m12                                            * conv_pressure,
	sim_info.full_stress_tensor().m13                                            * conv_pressure,
	sim_info.full_stress_tensor().m23                                            * conv_pressure,
  //       sim_info.deviator().x                                            * conv_pressure,
  //       sim_info.deviator().y                                            * conv_pressure,
	// sim_info.deviator().z                                            * conv_pressure,
	A, B, C, ALPHA, BETA, GAMMA, sim_info.volume(), conv_density * sim_info.mass()/sim_info.volume()) ;

      if( electronic_energy.has_value() )
      {
        oss << format_string(" % .7e",(*electronic_energy) * conv_energy / sim_info.particle_count() );
      }

      if( is_dump_virial ) {
        oss << format_string(" % .7e  % .7e  % .7e % .7e  % .7e  % .7e % .7e  % .7e  % .7e",
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
   OperatorNodeFactory::instance()->register_factory( "dump_thermodynamic_state", make_simple_operator<DumpThermodynamicStateNode> );
  }

}

