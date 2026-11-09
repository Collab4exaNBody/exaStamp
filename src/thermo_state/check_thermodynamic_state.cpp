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
#include <onika/file_utils.h>
#include <onika/string_utils.h>
#include <exaStamp/thermo_state/thermodynamic_state.h>

#include <yaml-cpp/yaml.h>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <utility>

namespace exaStamp
{
  using namespace exanb;

  // Compares the thermodynamic state (energies, virial, volume, particle count) to a reference file,
  // or writes the reference file when it does not exist and create_if_missing is true.
  class CheckThermodynamicStateNode : public OperatorNode
  {
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT , REQUIRED );
    ADD_SLOT( std::string        , file                , INPUT , REQUIRED , DocString{"Reference file (YAML)"} );
    ADD_SLOT( double             , rel_threshold       , INPUT , 1.e-10 , DocString{"Maximum relative error for each quantity"} );
    ADD_SLOT( bool               , create_if_missing   , INPUT , true , DocString{"Write the reference file when it does not exist (otherwise fatal error)"} );
    ADD_SLOT( bool               , fatal               , INPUT , true , DocString{"Fatal error when a quantity exceeds the threshold"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Compares the final thermodynamic state to a reference file : particle count (exact), potential, kinetic, rotational
and total energies, virial tensor (error relative to the tensor norm) and volume, in internal units.
If the reference file does not exist, it is written from the current state when create_if_missing is true.

Usage example:

simulation_epilog:
  - force_thermo_state
  - compute_all_forces_energy
  - default_thermostate_compute
  - check_thermodynamic_state:
      file: ref_thermo.yaml
      rel_threshold: 1.e-10
)EOF";
    }

    inline void execute () override final
    {
      const ThermodynamicState& ts = *thermodynamic_state;
      const Mat3d& w = ts.virial();
      const std::vector< std::pair<std::string,double> > scalars = {
        { "potential_energy"  , ts.potential_energy() } ,
        { "kinetic_energy"    , ts.kinetic_energy_scal() } ,
        { "rotational_energy" , ts.rotational_energy_scal() } ,
        { "total_energy"      , ts.total_energy() } ,
        { "volume"            , ts.volume() } };
      const std::vector<double> virial = { w.m11, w.m12, w.m13, w.m21, w.m22, w.m23, w.m31, w.m32, w.m33 };

      int rank = 0;
      MPI_Comm_rank( *mpi, &rank );
      int status = 0; // 0 = ok, 1 = mismatch, 2 = reference written, 3 = missing reference, 4 = bad reference
      std::string message;

      if( rank == 0 )
      {
        const std::string file_name = onika::data_file_path( *file );
        if( ! std::ifstream( file_name ).good() )
        {
          std::string non_finite;
          for(const auto& v : scalars) if( ! std::isfinite(v.second) ) non_finite += " " + v.first;
          for(double v : virial) if( ! std::isfinite(v) ) { non_finite += " virial"; break; }
          if( ! non_finite.empty() )
          {
            status = 4;
            message = "non-finite value(s) for" + non_finite + " , reference file " + file_name + " not written";
          }
          else if( *create_if_missing )
          {
            std::ofstream fout( file_name );
            fout << "particle_count: " << ts.particle_count() << "\n";
            for(const auto& s : scalars) fout << s.first << ": " << onika::format_string("%.17g",s.second) << "\n";
            fout << "virial: [";
            for(size_t i=0;i<virial.size();i++) fout << (i?", ":" ") << onika::format_string("%.17g",virial[i]);
            fout << " ]\n";
            status = fout.good() ? 2 : 4;
            message = "reference file " + file_name + ( fout.good() ? " written" : " could not be written" );
          }
          else
          {
            status = 3;
            message = "reference file " + file_name + " not found";
          }
        }
        else
        {
          try
          {
            const YAML::Node ref = YAML::LoadFile( file_name );
            lout << "check_thermodynamic_state : reference " << file_name << " , threshold " << *rel_threshold << std::endl;
            if( ref["particle_count"].as<size_t>() != ts.particle_count() )
            {
              status = 1;
              lout << onika::format_string("  %-17s : %ld (reference %ld)", "particle_count", long(ts.particle_count()), long(ref["particle_count"].as<size_t>())) << std::endl;
            }
            for(const auto& s : scalars)
            {
              const double r = ref[s.first].as<double>();
              const double err = ( r != 0.0 ) ? std::fabs( s.second - r ) / std::fabs( r ) : std::fabs( s.second );
              const bool ok = err <= *rel_threshold;
              if( ! ok ) status = 1;
              lout << onika::format_string("  %-17s : % .17e (reference % .17e) rel. error %.3e %s", s.first, s.second, r, err, ok ? "" : "<- FAILED") << std::endl;
            }
            const auto rv = ref["virial"].as< std::vector<double> >();
            double norm = 0.0, diff = 0.0;
            for(size_t i=0;i<virial.size();i++) { norm += rv.at(i)*rv.at(i); diff += (virial[i]-rv.at(i))*(virial[i]-rv.at(i)); }
            norm = std::sqrt(norm); diff = std::sqrt(diff);
            const double err = ( norm != 0.0 ) ? diff / norm : diff;
            const bool ok = err <= *rel_threshold;
            if( ! ok ) status = 1;
            lout << onika::format_string("  %-17s : |W| % .17e , |W-Wref|/|Wref| %.3e %s", "virial", norm, err, ok ? "" : "<- FAILED") << std::endl;
            if( status == 1 ) message = "thermodynamic state differs from reference file " + file_name;
          }
          catch(const std::exception& e)
          {
            status = 4;
            message = "cannot read reference file " + file_name + " : " + e.what();
          }
        }
      }

      MPI_Bcast( &status, 1, MPI_INT, 0, *mpi );
      if( status == 2 ) { lout << "check_thermodynamic_state : " << message << std::endl; }
      else if( status == 3 || status == 4 || ( status == 1 && *fatal ) )
      {
        if( rank == 0 ) fatal_error() << "check_thermodynamic_state : " << message << std::endl;
        else fatal_error() << "check_thermodynamic_state failed (see rank 0)" << std::endl;
      }
      else if( status == 1 ) { lerr << "check_thermodynamic_state : " << message << std::endl; }
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(check_thermodynamic_state)
  {
    OperatorNodeFactory::instance()->register_factory( "check_thermodynamic_state", make_compatible_operator< CheckThermodynamicStateNode > );
  }
}
