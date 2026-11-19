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
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/impropers_potentials_parameters.h>
#include <exaStamp/molecule/intramolecular_pair_weight.h>
#include <exaStamp/molecule/iFFname.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>
#include <algorithm>

namespace exaStamp
{
  using namespace exanb;

  class MoleculeForceField : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies              , species                  , INPUT , REQUIRED );
    ADD_SLOT( BondsPotentialParameters     , potentials_for_bonds     , INPUT_OUTPUT , BondsPotentialParameters{} , DocString{"Bond potentials of the force field"} );
    ADD_SLOT( BendsPotentialParameters     , potentials_for_angles    , INPUT_OUTPUT , BendsPotentialParameters{} , DocString{"Bending angle potentials of the force field"} );
    ADD_SLOT( TorsionsPotentialParameters  , potentials_for_torsions  , INPUT_OUTPUT , TorsionsPotentialParameters{} , DocString{"Torsion potentials of the force field"} );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT_OUTPUT , ImpropersPotentialParameters{} , DocString{"Improper torsion potentials of the force field"} );
    ADD_SLOT( LJExp6RFMultiParms           , potentials_for_pairs     , INPUT_OUTPUT , LJExp6RFMultiParms{} , DocString{"Pair potentials (LJ or Exp6, plus reaction field) for each pair of species"} );
    ADD_SLOT( IntramolecularPairWeighting  , mol_pair_weights         , INPUT_OUTPUT , IntramolecularPairWeighting{} , DocString{"intramolecular pair weights, keyed by species pair 'typeA_typeB'"} );
    ADD_SLOT( bool                         , rf_self_pairs            , INPUT_OUTPUT , true , DocString{"Add the reaction field self energy of each atom"} );
    ADD_SLOT( bool                         , long_range_correction    , INPUT_OUTPUT , true , DocString{"Long range correction of the pair energy and virial"} );
    ADD_SLOT( IntraFFtypes                 , ffnameVector             , OUTPUT , DocString{"Force field type names, indexed by force field type id"} );
    ADD_SLOT( double                       , rcut_max                 , INPUT_OUTPUT , 0.0 , DocString{"Maximum cutoff, updated with the pair potential cutoffs"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Flexible molecule force field given in the input file, for systems read without their force field
(e.g. read_xyz_molecules) : intramolecular potentials, intramolecular pair weights and pair potentials. Use it in
setup_system after the reader. It provides the same data as read_fatomes_mol reads from FAtomes files. Force field
types of intramolecular potentials are the species force field types (species key 'ff', species name by default).
)EOF";
    }

    inline void execute () override final
    {
      // force field type names, indexed by the force field type id of species
      ffnameVector->clear();
      for(const auto& sp : *species)
      {
        if( sp.m_FFtypeId < 0 ) fatal_error() << "molecule_force_field: species "<<sp.name()<<" has no force field type" << std::endl;
        if( size_t(sp.m_FFtypeId) >= ffnameVector->size() ) ffnameVector->resize( sp.m_FFtypeId + 1 );
        (*ffnameVector)[ sp.m_FFtypeId ] = sp.FFname();
      }

      // pair potential cutoff, as computed by read_fatomes_mol
      for(const auto& pot : potentials_for_pairs->m_potentials)
      {
        *rcut_max = std::max( *rcut_max , std::max( pot.m_params.m_rcut , pot.m_params.m_rf.rc ) );
      }

      lout << "======== Molecule force field ========" << std::endl
           << "Force field types = "<<ffnameVector->size() << std::endl
           << "Bond potentials   = "<<potentials_for_bonds->m_bond_desc.size() << std::endl
           << "Angle potentials  = "<<potentials_for_angles->m_potentials.size() << std::endl
           << "Torsion potent.   = "<<potentials_for_torsions->m_potentials.size() << std::endl
           << "Improper potent.  = "<<potentials_for_impropers->m_potentials.size() << std::endl
           << "Pair potentials   = "<<potentials_for_pairs->m_potentials.size() << std::endl
           << "Pair weights      = "<<mol_pair_weights->m_molecule_weight.size() << std::endl
           << "rcut_max          = "<<*rcut_max << std::endl
           << "======================================" << std::endl;
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_force_field)
  {
   OperatorNodeFactory::instance()->register_factory( "molecule_force_field", make_simple_operator<MoleculeForceField> );
  }

}
