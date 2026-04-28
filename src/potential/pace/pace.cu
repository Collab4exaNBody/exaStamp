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

#include <vector>
#include <memory>
#include <mpi.h>

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

#include "pace_params.h"
#include "pace_config.h"
#include "pace_force_op.h"

namespace exaStamp
{

  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class PaceForce : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                   , mpi                 , INPUT        , REQUIRED );
    ADD_SLOT( double                     , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors  , chunk_neighbors     , INPUT        , exanb::GridChunkNeighbors{}, DocString{"neighbor list"} );
    ADD_SLOT( bool                       , ghost               , INPUT        , true );
    ADD_SLOT( bool                       , conv_coef_units     , INPUT        , true );
    ADD_SLOT( bool                       , trigger_thermo_state, INPUT        , OPTIONAL );
    ADD_SLOT( GridT                      , grid                , INPUT_OUTPUT );
    ADD_SLOT( Domain                     , domain              , INPUT        , REQUIRED );
    ADD_SLOT( GridParticleLocks          , particle_locks      , INPUT        , OPTIONAL, DocString{"particle spin locks"} );
    ADD_SLOT( long                       , timestep            , INPUT        , REQUIRED, DocString{"Iteration number"} );
    ADD_SLOT( ParticleSpecies            , species             , INPUT        , REQUIRED );
    ADD_SLOT( ParticleTypeMap            , particle_type_map   , INPUT        , REQUIRED );
    ADD_SLOT( PaceContext                , pace_ctx            , INPUT );

    static constexpr bool UseWeights  = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights, UseNeighbors, PaceComputeBuffer, CopyParticleType>;

    using CellParticles = typename GridT::CellParticles;
    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;

    using ComputeFieldsWithoutVirial = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type, field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field, ComputeFieldsWithVirial, ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};

  public:

    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      const size_t nt = omp_get_max_threads();

      if( nt > pace_ctx->m_test.size() )
      {
        const size_t old_nt = pace_ctx->m_test.size();
        ldbg << "Resizing PACE thread context: " << old_nt << " -> " << nt << std::endl;
        pace_ctx->m_test.resize( nt );
        for(size_t i = old_nt; i < nt; i++)
        {
          assert( pace_ctx->m_test[i] == nullptr );
          pace_ctx->m_test[i] = std::make_shared<ACEImpl>();
          pace_ctx->m_test[i]->basis_set = std::make_unique<ACECTildeBasisSet>( *pace_ctx->aceimpl->basis_set );
          pace_ctx->m_test[i]->ace       = std::make_unique<ACERecursiveEvaluator>();
          pace_ctx->m_test[i]->ace->set_recursive(pace_ctx->recursive);
          pace_ctx->m_test[i]->ace->element_type_mapping.init(pace_ctx->nspecies + 1);
          for (int j = 1; j <= pace_ctx->nspecies; j++)
            pace_ctx->m_test[i]->ace->element_type_mapping(j) = pace_ctx->aceimpl->ace->element_type_mapping(j);
          pace_ctx->m_test[i]->ace->set_basis( *pace_ctx->m_test[i]->basis_set, 1 );
        }
      }

      if( grid->number_of_cells() == 0 ) return;

      if( !particle_locks.has_value() )
      {
        fatal_error() << "No particle locks" << std::endl;
      }

      bool eflag = false;
      if( trigger_thermo_state.has_value() )
      {
        ldbg << "trigger_thermo_state = " << *trigger_thermo_state << std::endl;
        eflag = *trigger_thermo_state;
      }
      else
      {
        ldbg << "trigger_thermo_state missing" << std::endl;
      }

      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();
      LinearXForm cp_xform{ domain->xform() };

      auto compute_opt_locks = [&](auto cp_locks)
      {
        PaceForceOp force_op{ pace_ctx->m_test, !(*conv_coef_units), eflag };
        compute_cell_particle_pairs( *grid, *rcut_max, *ghost, make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks), force_buf, force_op, compute_force_field_set, parallel_execution_context() );
      };

      if( omp_get_max_threads() > 1 ) compute_opt_locks( ComputePairOptionalLocks<true>{ particle_locks->data() } );
      else                            compute_opt_locks( ComputePairOptionalLocks<false>{} );
    }

  };

  template<class GridT> using PaceForceTmpl = PaceForce<GridT>;

  ONIKA_AUTORUN_INIT(pace)
  {
    OperatorNodeFactory::instance()->register_factory( "pace_force", make_grid_variant_operator<PaceForceTmpl> );
  }

}
