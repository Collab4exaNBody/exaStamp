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

#include "pod_params.h"
#include "pod_config.h"
#include "pod_force_op.h"

namespace exaStamp
{

  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class PodForce : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                  , mpi                 , INPUT        , REQUIRED );
    ADD_SLOT( double                    , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors     , INPUT        , exanb::GridChunkNeighbors{}, DocString{"neighbor list"} );
    ADD_SLOT( bool                      , ghost               , INPUT        , false );
    ADD_SLOT( bool                      , conv_coef_units     , INPUT        , true );
    ADD_SLOT( bool                      , trigger_thermo_state, INPUT        , OPTIONAL );
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain              , INPUT        , REQUIRED );
    ADD_SLOT( GridParticleLocks         , particle_locks      , INPUT        , OPTIONAL, DocString{"particle spin locks"} );
    ADD_SLOT( long                      , timestep            , INPUT        , REQUIRED, DocString{"Iteration number"} );
    ADD_SLOT( ParticleSpecies           , species             , INPUT        , REQUIRED );
    ADD_SLOT( ParticleTypeMap           , particle_type_map   , INPUT        , REQUIRED );
    ADD_SLOT( PodContext                , pod_ctx             , INPUT );

    static constexpr bool UseWeights   = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights, UseNeighbors, PodComputeBuffer, CopyParticleType>;

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
      const size_t nt = omp_get_max_threads();

      // Thread contexts are pre-allocated in pod_init. Warn if the runtime
      // thread count exceeds what was available at init time.
      if (nt > pod_ctx->m_eapod.size()) {
        lerr << "POD: omp_get_max_threads() grew from " << pod_ctx->m_eapod.size()
             << " to " << nt << " after init — some threads lack an EAPOD context."
             << " Re-run with the correct OMP_NUM_THREADS set before launch." << std::endl;
        fatal_error() << "POD thread context size mismatch" << std::endl;
      }

      if (grid->number_of_cells() == 0) return;

      if (!particle_locks.has_value()) {
        fatal_error() << "No particle locks" << std::endl;
      }

      bool eflag = false;
      if (trigger_thermo_state.has_value()) {
        ldbg << "trigger_thermo_state = " << *trigger_thermo_state << std::endl;
        eflag = *trigger_thermo_state;
      } else {
        ldbg << "trigger_thermo_state missing" << std::endl;
      }

      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();
      LinearXForm cp_xform{ domain->xform() };

      auto compute_opt_locks = [&](auto cp_locks)
      {
        PodForceOp force_op{ pod_ctx->m_eapod, pod_ctx->type_map, !(*conv_coef_units), eflag };
        compute_cell_particle_pairs(
            *grid, *rcut_max, *ghost,
            make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks),
            force_buf, force_op, compute_force_field_set,
            parallel_execution_context());
      };

      if (omp_get_max_threads() > 1) compute_opt_locks(ComputePairOptionalLocks<true>{ particle_locks->data() });
      else                           compute_opt_locks(ComputePairOptionalLocks<false>{});
    }

  };

  template<class GridT> using PodForceTmpl = PodForce<GridT>;

  ONIKA_AUTORUN_INIT(pod)
  {
    OperatorNodeFactory::instance()->register_factory("pod_force", make_grid_variant_operator<PodForceTmpl>);
  }

}
