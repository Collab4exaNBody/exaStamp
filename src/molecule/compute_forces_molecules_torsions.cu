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

#include <yaml-cpp/yaml.h>
#include <memory>
#include <utility>// std::pair

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/math/basic_types_yaml.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/log.h>
#include <onika/thread.h>  // GridParticleLocks

#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exaStamp/molecule/torsions_force_functor.h>

namespace exaStamp
{

  ///**
  // * Signum function : give the sign of an expression (+1, -1 or 0)
  // */
  //template <typename T>
  //typename std::enable_if<std::is_unsigned<T>::value, int>::type
  //inline constexpr signum(T x) noexcept {
  //  return T(0) < x;
  //}

  //template <typename T>
  //typename std::enable_if<std::is_signed<T>::value, int>::type
  //inline constexpr signum(T x) noexcept {
  //  return (T(0) < x) - (x < T(0));
  //}


  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_ep >
    >
  class ComputeForcesTorsionsNode : public OperatorNode
  {
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( GridT                       , grid                    , INPUT_OUTPUT );
    ADD_SLOT( Domain                      , domain                  , INPUT );
    ADD_SLOT( ChemicalTorsions            , chemical_torsions       , INPUT, OPTIONAL );
    ADD_SLOT( TorsionsPotentialParameters , potentials_for_torsions , INPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies             , species                 , INPUT, REQUIRED );
    
    ADD_SLOT( GridParticleLocks        , particle_locks        , INPUT_OUTPUT);

    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( IntramolecularParameterIndexLists, intramolecular_parameters , INPUT_OUTPUT, DocString{"Intramolecular functional parmater index lists"} );

    ADD_SLOT( bool                      , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( bool                      , compute_virial      , INPUT , false );

    inline void execute ()  override final
    {
      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_torsions.has_value() )
      {
        fatal_error() << "chemical_torsions input missing" << std::endl;
      }
      
      // do we need to compute energy and virial
      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }
      const bool need_virial = log_energy && *compute_virial;
      using VirialFieldT = decltype( grid->field_accessor( field::virial ) );
      VirialFieldT virial_field = {};
      if( need_virial ) virial_field = grid->field_accessor( field::virial );

      ldbg<<"n_torsions = "<<chemical_torsions->size()<<std::endl;

      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 

#     ifndef NDEBUG
      if( ! domain->xform_is_identity() )
      {
        Vec3d tmp = domain->xform() * size_box;
        if( fabs(tmp.x) <= 1.5 || fabs(tmp.y) <= 1.5 || fabs(tmp.z) <= 1.5 )
        {
          lerr<<"xform="<<domain->xform()<<", size_box="<<size_box<<", tmp="<<tmp<<std::endl;
          std::abort();
        }
      }
#     endif

      auto compute_torsions_force_opt_virial = [&]( auto cpvir )
      {
        auto cells = grid->cells_accessor();
        TorsionForceOp<decltype(cells),VirialFieldT,cpvir.value> torsion_force_op =
          { cells
          , particle_locks->data()
          , molecule_compute_parameters->m_func_params.data()
          , intramolecular_parameters->m_torsion_param_idx.data()
          , chemical_torsions->data()
          , size_box , half_min_size_box
          , domain->xform() , domain->xform_is_identity()
          , virial_field };

        /* auto parallel_op = */ parallel_for( chemical_torsions->size(), torsion_force_op, parallel_execution_context("torsion_force") );
        // auto exec_ctrl_obj = parallel_execution_stream() << std::move(parallel_op) ;
      };

      if( need_virial ) compute_torsions_force_opt_virial( onika::TrueType{} );
      else              compute_torsions_force_opt_virial( onika::FalseType{} );

      ldbg<<"compute_force_torsion done"<<std::endl;
    }

  };


  template<class GridT> using ComputeForcesTorsionsNodeTmpl = ComputeForcesTorsionsNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_forces_molecules_torsions)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_force_torsion", make_grid_variant_operator< ComputeForcesTorsionsNodeTmpl> );
  }

}


