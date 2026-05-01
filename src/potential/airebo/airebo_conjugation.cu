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
#include "airebo_conjugation_op.h"

namespace exaStamp
{
  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_type >
    >
  class AireboConjugation : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                  , mpi                 , INPUT        , REQUIRED );
    ADD_SLOT( double                    , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                    , bondorder_cutoff    , INPUT        , REQUIRED );    
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors     , INPUT        , exanb::GridChunkNeighbors{}, DocString{"neighbor list"} );
    ADD_SLOT( bool                      , ghost               , INPUT        , false );
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain              , INPUT        , REQUIRED );
    ADD_SLOT( GridParticleLocks         , particle_locks      , INPUT        , OPTIONAL, DocString{"particle spin locks"} );
    ADD_SLOT( long                      , timestep            , INPUT        , REQUIRED, DocString{"Iteration number"} );
    ADD_SLOT( ParticleSpecies           , species             , INPUT        , REQUIRED );
    ADD_SLOT( AireboParams              , parameters          , INPUT        , REQUIRED );

    using ComputeBufferNconj     = ComputePairBuffer2<false,false,NconjOpExtStorage>;

  public:

    inline void execute() override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      if ( grid->number_of_cells() == 0 ) return;

      // Read Only AIREBO parameters
      const AireboParamsRO params_ro{ *parameters };
      // ------------------------------------------------------------------------------ //
      // First pass: we compute NijC and NijH bond order per atom
      // ------------------------------------------------------------------------------ //
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf_nconj = make_compute_pair_buffer<ComputeBufferNconj>();
      
      // Accessor for NijC, NijH and Nconj
      auto nijc_acc  = grid->field_accessor( field::mk_generic_real( "NijC") );
      auto nijh_acc  = grid->field_accessor( field::mk_generic_real( "NijH") );
      auto nconj_acc = grid->field_accessor( field::mk_generic_real( "Nconj") );

      // Conjugation operation declaration
      NconjOp<decltype(nijc_acc),decltype(nijh_acc),decltype(nconj_acc)> compute_op_nconj = { &params_ro , nijc_acc, nijh_acc, nconj_acc };
      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr std::integral_constant<bool,true> force_use_cells_accessor = {};
      static constexpr DefaultPositionFields posfields = {};
      compute_cell_particle_pairs2( *grid, *bondorder_cutoff, false, optional, compute_buf_nconj, compute_op_nconj, compute_field_set
                                  , posfields, parallel_execution_context(), force_use_cells_accessor );
      // ------------------------------------------------------------------------------ //      

    }
  };

  template<class GridT> using AireboConjugationTmpl = AireboConjugation<GridT>;

  ONIKA_AUTORUN_INIT(airebo_conjugation)
  {
    OperatorNodeFactory::instance()->register_factory(
        "airebo_conjugation", make_grid_variant_operator<AireboConjugationTmpl>);
  }

}
