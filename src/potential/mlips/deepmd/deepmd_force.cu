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
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/file_utils.h>

#include "DeepPot.h"
#include <exaStamp/potential/mlips/deepmd/deepmd.h>
#include <exaStamp/potential/mlips/deepmd/deepmd_force_op.h>

#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid_fields.h>

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>

#include <onika/memory/allocator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/source_term.h>

#include <mpi.h>
#include <iomanip>
#include <memory>

namespace exaStamp
{
  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz, field::_type >
    >
  class DeepMDForce : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi          , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( ParticleSpecies, species      , INPUT , REQUIRED );
    ADD_SLOT( GridChunkNeighbors    , chunk_neighbors   , INPUT        , GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT        , false    );    
    ADD_SLOT( double         , rcut_max     , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( GridT  , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain , domain       , INPUT , REQUIRED );
    ADD_SLOT( long           , grid_subdiv  , INPUT , 3 );
    ADD_SLOT( GridCellValues , grid_cell_values      , INPUT_OUTPUT );
    ADD_SLOT( deepmd::DeepPot , deep_pot     , INPUT );
    ADD_SLOT( std::string , model   , INPUT_OUTPUT , REQUIRED );    
    ADD_SLOT( DPMDContext , dpmd_ctx , PRIVATE );
    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT , OPTIONAL , DocString{"particle spin locks"} );
    
    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors>;
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type, field::_virial >;
    using ComputeFields              = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};
        
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      ldbg << "DeepMD force computation" << std::endl;

      if( grid->number_of_cells() == 0 ) { return; }
      
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      size_t nt = omp_get_max_threads();
      if (nt > dpmd_ctx->m_thread_ctx.size()) {
        size_t old_nt = dpmd_ctx->m_thread_ctx.size();
        dpmd_ctx->m_thread_ctx.resize( nt );

        for(size_t j=old_nt;j<nt;j++)
          {
            assert( dpmd_ctx->m_thread_ctx[j].dpmd_model == nullptr );
            dpmd_ctx->m_thread_ctx[j].dpmd_model = new deepmd::DeepPot (*model);
          }
      }

      *rcut_max = std::max( *rcut_max, (*deep_pot).cutoff() );

      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
		
      if( ! particle_locks.has_value() )
      {
        fatal_error() << "No particle locks" << std::endl;
      }

      ComputePairNullWeightIterator cp_weight{};
      GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();
      LinearXForm cp_xform { domain->xform() };

      using TypeFieldT = decltype( grid->field_accessor( field::type ) );
      TypeFieldT   type_field   = {};
      type_field = grid->field_accessor( field::type );
      
      auto compute_opt_locks = [&](auto cp_locks)
      {
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
        DPMDForceOp<TypeFieldT> force_op { dpmd_ctx->m_thread_ctx.data() , type_field };
        compute_cell_particle_pairs( *grid, (*deep_pot).cutoff(), *ghost, optional, force_buf, force_op , compute_force_field_set , parallel_execution_context() );
      };

      if( omp_get_max_threads() > 1 ) {
        compute_opt_locks( ComputePairOptionalLocks<true>{ particle_locks->data() } );
      } else {
        compute_opt_locks( ComputePairOptionalLocks<false>{} );
      }
      
    }

  };

  template<class GridT> using DeepMDForceTmpl = DeepMDForce<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(deepmd_force)
  {
   OperatorNodeFactory::instance()->register_factory("deepmd_force", make_grid_variant_operator< DeepMDForceTmpl > );
  }

}
