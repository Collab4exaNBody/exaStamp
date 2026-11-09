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

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include "sw_force_op.h"
#include "sw_params.h"

namespace exaStamp
{
  using namespace exanb;

  template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz, field::_type>>
  class SWForce : public OperatorNode
  {
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT        , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                      , ghost           , INPUT        , false , DocString{"compute central atoms in the ghost layer too (not needed with ghost force fold-back)"} );
    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED );
    ADD_SLOT( SWContext                 , sw_ctx          , INPUT        , REQUIRED , DocString{"Stillinger-Weber tables built by sw_init"} );
    ADD_SLOT( GridParticleLocks         , particle_locks  , INPUT        , OPTIONAL , DocString{"particle spin locks (unused, updates are atomic)"} );

    static constexpr bool UseWeights   = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights, UseNeighbors, SWComputeBuffer, SWCopyParticleType, SW_MAX_NEIGHBORS>;

    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type, field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field, ComputeFieldsWithVirial, ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Stillinger-Weber many-body potential (two-body + three-body angular term), multi-species.
Parameters come from sw_init. Forces are atomically scattered to neighbors, including ghost
particles, so forces and energies must be zeroed in the ghost layer before the force computation
and ghost contributions folded back afterwards :

compute_force_prolog:
  - zero_force_energy: { ghost: true }
compute_force_epilog:
  - update_force_energy_from_ghost
  - force_to_accel

compute_force:
  - sw_force

Per-atom energy : half of each pair energy, and the whole energy of each triplet to its central atom.

)EOF";
    }

    inline void execute() override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      if( grid->number_of_cells() == 0 ) return;

      if( sw_ctx->ntypes <= 0 || sw_ctx->triplets.empty() )
      {
        fatal_error() << "sw_force: Stillinger-Weber tables are empty, sw_init must run first" << std::endl;
      }

      SWForceOp force_op{ sw_ctx->triplets.data(), sw_ctx->ntypes, sw_ctx->threebody };

      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      ComputePairNullWeightIterator cp_weight{};
      LinearXForm cp_xform{ domain->xform() };
      ComputePairOptionalLocks<false> cp_locks{};
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();

      compute_cell_particle_pairs( *grid, sw_ctx->rcut, *ghost,
                                   make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks ),
                                   force_buf, force_op, compute_force_field_set,
                                   parallel_execution_context() );
    }
  };

  template<class GridT> using SWForceTmpl = SWForce<GridT>;

  ONIKA_AUTORUN_INIT(sw_force)
  {
    OperatorNodeFactory::instance()->register_factory( "sw_force", make_grid_variant_operator<SWForceTmpl> );
  }

}
