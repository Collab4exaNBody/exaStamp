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

//  // DO NOT REMOVE THIS LINE

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>

#include <exaStamp/molecule/pairs_force_functor.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/math/basic_types_yaml.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/log.h>

#include <string>
#include <iostream>
#include <sstream>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class IntramolecularPairForce : public OperatorNode
  {

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( Domain                  , domain                , INPUT );
    ADD_SLOT( ChemicalPairs           , chemical_pairs        , INPUT_OUTPUT, OPTIONAL , DocString{"Intramolecular pairs (atoms 1-2, 1-3 and 1-4 apart)"} );
    ADD_SLOT( ParticleSpecies         , species               , INPUT, REQUIRED );
    ADD_SLOT( BondsPotentialParameters, potentials_for_bonds  , INPUT_OUTPUT, REQUIRED , DocString{"Bond potentials of the force field"} );
    ADD_SLOT( GridT                   , grid                  , INPUT_OUTPUT );
    ADD_SLOT( GridParticleLocks       , particle_locks        , INPUT_OUTPUT , DocString{"Particle locks for concurrent force updates"} );

    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( IntramolecularParameterIndexLists, intramolecular_parameters , INPUT_OUTPUT, DocString{"Intramolecular functional parameter index lists"} );

    ADD_SLOT( bool                      , trigger_thermo_state, INPUT , OPTIONAL , DocString{"If set, energy and virial are only computed when true"} );
    ADD_SLOT( bool                      , per_atom_charge     , INPUT , true , DocString{"Use the charge field instead of species charges"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Computes the pair forces (pair potential and reaction field) of intramolecular pairs : atoms 1-2, 1-3 and 1-4 apart in the
bond graph (chemical_pairs, built by intramolecular_pair_list). These pairs are excluded from the neighbor lists
(extramolecular_neighbors). The pair potential and the reaction field of each pair are scaled by the bond, bend or torsion
weights of mol_pair_weights. Energy and virial are only computed on thermodynamic state steps (trigger_thermo_state).
)EOF";
    }

    inline void execute() override final
    {
      using onika::parallel::parallel_for;

      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_pairs.has_value() )
      {
        lerr << "chemical_pairs input missing" << std::endl;
        std::abort();
      }

      // energy and virial are computed on thermodynamic state steps only
      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }
      const bool need_virial = log_energy;
      // per particle energy is only needed when thermodynamic state is computed. if unknown, keep computing it
      const bool need_energy = trigger_thermo_state.has_value() ? *trigger_thermo_state : true;
      using ChargeFieldT = decltype( grid->field_accessor( field::charge ) );
      using VirialFieldT = decltype( grid->field_accessor( field::virial ) );
      ChargeFieldT charge_field = {};
      VirialFieldT virial_field = {};
      if( *per_atom_charge )
      {
        charge_field = grid->field_accessor( field::charge );
      }
      if( need_virial )
      {
        virial_field = grid->field_accessor( field::virial );
      }
      
      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};


      ldbg<<"n_pairs = "<<chemical_pairs->size()<<" , pair params = "<< molecule_compute_parameters->m_pair_params.size()<<std::endl;

      auto compute_pair_force_opt_virial = [&]( auto cpene, auto cpvir )
      {
        auto cells = grid->cells_accessor();
        IntramolecularPairForceOp<decltype(cells),ChargeFieldT,VirialFieldT,cpene.value,cpvir.value> pair_force_op =
          { cells
          , particle_locks->data()
          , molecule_compute_parameters->m_pair_params.data()
          , intramolecular_parameters->m_pair_param_idx.data()
          , chemical_pairs->data()
          , species->data()
          , intramolecular_parameters->m_pair_image.data() , size_box
          , domain->xform() , domain->xform_is_identity()
          , *per_atom_charge
          , charge_field, virial_field };

        /* auto parallel_op = */ parallel_for( chemical_pairs->size(), pair_force_op, parallel_execution_context("intramol_pair_force") );
      };

      if( need_virial )      compute_pair_force_opt_virial( onika::TrueType{}  , onika::TrueType{}  );
      else if( need_energy ) compute_pair_force_opt_virial( onika::TrueType{}  , onika::FalseType{} );
      else                   compute_pair_force_opt_virial( onika::FalseType{} , onika::FalseType{} );

      ldbg<<"intramolecular_pair_force done"<<std::endl;
    }

  };

  template<class GridT> using IntramolecularPairForceTmpl = IntramolecularPairForce<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(intramolecular_pair_force)
  {
    OperatorNodeFactory::instance()->register_factory( "intramolecular_pair_force", make_grid_variant_operator< IntramolecularPairForceTmpl > );
  }

}

