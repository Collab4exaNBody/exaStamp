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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/core/particle_type_id.h>

#include <memory>
#include <vector>
#include <mpi.h>

#include "airebo_params.h"
#include "airebo_force_op_cached.h"
#include "lj_force_op.h"
#include "torsion_force_op.h"

namespace exaStamp
{
  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class AireboForce : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                  , mpi                 , INPUT        , REQUIRED );
    ADD_SLOT( double                    , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                    , rebo_cutoff         , INPUT        , REQUIRED );
    ADD_SLOT( double                    , lj_cutoff           , INPUT        , REQUIRED );       
    ADD_SLOT( double                    , torsion_cutoff      , INPUT        , REQUIRED );           
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors     , INPUT        , exanb::GridChunkNeighbors{}, DocString{"neighbor list"} );
    ADD_SLOT( bool                      , ghost               , INPUT        , false );
    ADD_SLOT( bool                      , conv_coef_units     , INPUT        , true );
    ADD_SLOT( bool                      , trigger_thermo_state, INPUT        , OPTIONAL );
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain              , INPUT        , REQUIRED );
    ADD_SLOT( GridParticleLocks         , particle_locks      , INPUT        , OPTIONAL, DocString{"particle spin locks"} );
    ADD_SLOT( long                      , timestep            , INPUT        , REQUIRED, DocString{"Iteration number"} );
    ADD_SLOT( ParticleSpecies           , species             , INPUT        , REQUIRED );
    ADD_SLOT( AireboParams              , parameters          , INPUT        , REQUIRED );

    static constexpr bool UseWeights   = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBufferAiRebo = ComputePairBuffer2<UseWeights, UseNeighbors,
                                             AireboComputeBuffer, CopyParticleType,
                                             AIREBO_MAX_NEIGHBORS>;

    using CellParticles = typename GridT::CellParticles;
    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;

    using ComputeFieldsWithoutVirial = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type, field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field, ComputeFieldsWithVirial, ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};

  public:

    inline void execute() override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      if ( grid->number_of_cells() == 0 ) return;

      bool eflag = false;
      if (trigger_thermo_state.has_value()) {
        ldbg << "trigger_thermo_state = " << *trigger_thermo_state << std::endl;
        eflag = *trigger_thermo_state;
      }

      // Read Only AIREBO parameters
      const AireboParamsRO params_ro{ *parameters };

      // ------------------------------------------------------------------------------ //
      // First pass: we compute NijC and NijH bond order per atom
      // ------------------------------------------------------------------------------ //
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      // auto compute_buf_bondorder = make_compute_pair_buffer<ComputeBufferBondOrder>();     
      // Accessor for NijC and NijH
      auto nijc_acc  = grid->field_accessor( field::mk_generic_real( "NijC") );
      auto nijh_acc  = grid->field_accessor( field::mk_generic_real( "NijH") );
      auto nconj_acc = grid->field_accessor( field::mk_generic_real( "Nconj") );

      LinearXForm cp_xform { domain->xform() };
 
      ComputePairNullWeightIterator cp_weight{};
      auto compute_buf_rebo = make_compute_pair_buffer<ComputeBufferAiRebo>();

      const bool use_locks = omp_get_max_threads() > 1 && particle_locks.has_value();
      auto force_op_fields = make_field_tuple_from_field_set( ComputeFields{} , nijc_acc, nijh_acc, nconj_acc );
      // This maybe a lead for later 
      // auto force_nbh_fields = onika::make_flat_tuple(nconj_acc);
      // auto force_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks, ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{}, force_nbh_fields );
      // compute_cell_particle_pairs( *grid, *rcut, false,
      //                          force_optional, force_buf, force_op, force_op_fields , parallel_execution_context() );
        
      // REBO: functor receives nC_central, nH_central, Nconj_central as extra args.
      auto run_rebo = [&](auto force_op)
      {
        auto with_locks = [&](auto cp_locks)
        {
          compute_cell_particle_pairs(
              *grid, *rebo_cutoff, *ghost,
              make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks),
              compute_buf_rebo, force_op, force_op_fields,
              parallel_execution_context());
        };
        if (use_locks) with_locks(ComputePairOptionalLocks<true>{ particle_locks->data() });
        else           with_locks(ComputePairOptionalLocks<false>{});
      };

      // LJ: functor receives nC_central, nH_central, Nconj_central as extra args.
      auto run_lj = [&](auto force_op)
      {
        auto with_locks = [&](auto cp_locks)
        {
          compute_cell_particle_pairs(
              *grid, *lj_cutoff, *ghost,
              make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks),
              compute_buf_rebo, force_op, force_op_fields,
              parallel_execution_context());
        };
        if (use_locks) with_locks(ComputePairOptionalLocks<true>{ particle_locks->data() });
        else           with_locks(ComputePairOptionalLocks<false>{});
      };      

      // LJ / torsion: standard field set, no extra per-particle fields.
      auto run = [&](auto force_op)
      {
        auto with_locks = [&](auto cp_locks)
        {
          compute_cell_particle_pairs(
              *grid, *torsion_cutoff, *ghost,
              make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks),
              compute_buf_rebo, force_op, compute_force_field_set,
              parallel_execution_context());
        };
        if (use_locks) with_locks(ComputePairOptionalLocks<true>{ particle_locks->data() });
        else           with_locks(ComputePairOptionalLocks<false>{});
      };

      run_rebo( AIREBOForceOp<decltype(nijc_acc),decltype(nijh_acc),decltype(nconj_acc)> { &params_ro, nijc_acc, nijh_acc, nconj_acc } );
      if (parameters->ljflag)  run_lj( LJForceOp<decltype(nijc_acc),decltype(nijh_acc),decltype(nconj_acc)> { &params_ro, nijc_acc, nijh_acc, nconj_acc } );
      if (parameters->torflag) run( TorsionForceOp  { &params_ro } );
    }
  };

  template<class GridT> using AireboForceTmpl = AireboForce<GridT>;

  ONIKA_AUTORUN_INIT(airebo_force)
  {
    OperatorNodeFactory::instance()->register_factory(
        "airebo_force", make_grid_variant_operator<AireboForceTmpl>);
  }

}
