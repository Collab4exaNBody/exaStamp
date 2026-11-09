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

// Names of former coulombic operators, so that older input files get a clear message instead of "unknown operator".

#include <exaStamp/operator_alias.h>
#include <onika/cpp_utils.h>

namespace exaStamp
{
  ONIKA_AUTORUN_INIT(coulombic_deprecated_names)
  {
    // same operator, same parameters
    register_deprecated_operator_alias( "reaction_field" , "coulombic_rf" );

    // replaced by operators with different parameters or behaviour
    register_removed_operator( "coul_wolf_pc" , "Use coulombic_wolf (same parameters, includes the self energy : remove coul_wolf_self)." );
    register_removed_operator( "coul_dsf_pc" , "Use coulombic_dsf (same parameters, includes the self energy)." );
    register_removed_operator( "coul_wolf_self" , "coulombic_wolf includes the self energy. With the coul_wolf pair style, use coulombic_wolf_self with per_atom_charge: false." );
    register_removed_operator( "ewald_init" , "Use coulombic_ewald_init (Ewald summation) or coulombic_pppm_init (PPPM), see the Electrostatics documentation." );
    register_removed_operator( "ewald_short_range_pc" , "Use coulombic_ewald_short_range, which reads its parameters from coulombic_ewald_init or coulombic_pppm_init." );
    register_removed_operator( "ewald_long_range" , "Use coulombic_ewald_long_range (or coulombic_pppm), set up by coulombic_ewald_init (or coulombic_pppm_init)." );
    register_removed_operator( "ewald_long_range_pc" , "Use coulombic_ewald_long_range (or coulombic_pppm), set up by coulombic_ewald_init (or coulombic_pppm_init)." );
    register_removed_operator( "ewald_potential_energy_shift" , "The self and neutralizing background energies are included in coulombic_ewald_long_range and coulombic_pppm." );
    for( const char* suffix : { "_compute_force" , "_compute_force_symetric" , "_multi_force" , "_plot" } )
    {
      register_removed_operator( std::string("ewald_short_range") + suffix , "Use coulombic_ewald_short_range, which reads its parameters from coulombic_ewald_init or coulombic_pppm_init." );
    }
  }
}
