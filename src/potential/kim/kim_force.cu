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

#include "kim.h"
#include "kim_force_op.h"

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/cuda/cuda_context.h>

#include <vector>
#include <memory>
#include <iostream>

namespace exaStamp
{
  using onika::memory::DEFAULT_ALIGNMENT;
  
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type >
    >
  class KIMComputeForce : public OperatorNode
  {

    using CellParticles = typename GridT::CellParticles;
    
    // ========= I/O slots =======================
    ADD_SLOT( KIMParams             , parameters        , INPUT        , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0      );
    ADD_SLOT( ParticleSpecies       , species           , INPUT        , REQUIRED );
    ADD_SLOT( GridChunkNeighbors    , chunk_neighbors   , INPUT        , GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT        , false    );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT            );
    ADD_SLOT( Domain                , domain            , INPUT        , REQUIRED );
    ADD_SLOT( KIMContext            , kim_ctx           , INPUT );    
    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT , OPTIONAL , DocString{"particle spin locks"} );
    ADD_SLOT( std::vector<int>      , kim_particle_codes, INPUT);
    
    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    //    using ComputeBuffer = ComputePairBuffer2<false,false>;
    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors,KimComputeBuffer,CopyParticleType>;
    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type, field::_virial >;
    using ComputeFields              = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};
    
  public:
    // Operator execution
    inline void execute () override final
    {
      KIM::Log::PushDefaultVerbosity(KIM::LOG_VERBOSITY::silent);

      ldbg << "KIM model = " << parameters->model << std::endl;
      ldbg << "KIM rcut  = " << parameters->rcut << std::endl;      
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      size_t nt = omp_get_max_threads();
      if (nt > kim_ctx->m_thread_ctx.size()) {
        size_t old_nt = kim_ctx->m_thread_ctx.size();
        kim_ctx->m_thread_ctx.resize( nt );
        int requestedUnitsAccepted;
        int error;
        for(size_t j=old_nt;j<nt;j++)
          {
            assert( kim_ctx->m_thread_ctx[j].kim_model == nullptr );
            error = KIM::Model::Create(KIM::NUMBERING::zeroBased,
                                       KIM::LENGTH_UNIT::A,
                                       KIM::ENERGY_UNIT::eV,
                                       KIM::CHARGE_UNIT::e,
                                       KIM::TEMPERATURE_UNIT::K,
                                       KIM::TIME_UNIT::ps,
                                       parameters->model,
                                       &requestedUnitsAccepted,
                                       &kim_ctx->m_thread_ctx[j].kim_model);
            if (error) { MY_ERROR("KIM::Model::Create()"); }
            // Check for compatibility with the model
            if (!requestedUnitsAccepted) { MY_ERROR("Must Adapt to model units"); }

            // This is done in kim_init so if it does get through there no need to do it here.
            // If an error occur, it should crash at kim_init
            
            // // Check that we know about all required routines
            // int numberOfModelRoutineNames;
            // KIM::MODEL_ROUTINE_NAME::GetNumberOfModelRoutineNames(&numberOfModelRoutineNames);
      
            // for (int i = 0; i < numberOfModelRoutineNames; ++i)
            //   {
            //     KIM::ModelRoutineName modelRoutineName;
            //     int error = KIM::MODEL_ROUTINE_NAME::GetModelRoutineName(i, &modelRoutineName);
            //     if (error) { MY_ERROR("Unable to get ModelRoutineName."); }
            //     int present;
            //     int required;
            //     error = kim_ctx->m_thread_ctx[j].kim_model->IsRoutinePresent(modelRoutineName, &present, &required);
            //     if (error) { MY_ERROR("Unable to get routine present/required."); }

            //     ldbg << "Model routine name \"" << modelRoutineName.ToString()
            //          << "\" has present = " << present
            //          << " and required = " << required << "." << std::endl;

            //     if ((present == true) && (required == true))
            //       {
            //         using namespace KIM::MODEL_ROUTINE_NAME;
            //         if (!((modelRoutineName == Create)
            //               || (modelRoutineName == ComputeArgumentsCreate)
            //               || (modelRoutineName == Compute) || (modelRoutineName == Refresh)
            //               || (modelRoutineName == ComputeArgumentsDestroy)
            //               || (modelRoutineName == Destroy)))
            //           {
            //             MY_ERROR("Unknown Routine \"" + modelRoutineName.ToString()
            //                      + "\" is required by model.");
            //           }
            //       }
            //   }
            
            // // print model units
            // KIM::LengthUnit lengthUnit;
            // KIM::EnergyUnit energyUnit;
            // KIM::ChargeUnit chargeUnit;
            // KIM::TemperatureUnit temperatureUnit;
            // KIM::TimeUnit timeUnit;

            // kim_ctx->m_thread_ctx[j].kim_model->GetUnits(&lengthUnit, &energyUnit, &chargeUnit, &temperatureUnit, &timeUnit);

            // ldbg << "LengthUnit is \"" << lengthUnit.ToString() << "\"" << std::endl
            //      << "EnergyUnit is \"" << energyUnit.ToString() << "\"" << std::endl
            //      << "ChargeUnit is \"" << chargeUnit.ToString() << "\"" << std::endl
            //      << "TemperatureUnit is \"" << temperatureUnit.ToString() << "\""
            //      << std::endl
            //      << "TimeUnit is \"" << timeUnit.ToString() << "\"" << std::endl;

            // // check species
            // int speciesIsSupported;
            // int modelCode;
            // error = kim_ctx->m_thread_ctx[j].kim_model->GetSpeciesSupportAndCode(KIM::SPECIES_NAME::Ta, &speciesIsSupported, &modelCode);
            // if ((error) || (!speciesIsSupported))
            //   {
            //     MY_ERROR("Species Ta not supported");
            //   }

            KIM::ComputeArguments * computeArguments;
            error = kim_ctx->m_thread_ctx[j].kim_model->ComputeArgumentsCreate(&computeArguments);
            if (error) { MY_ERROR("Unable to create a ComputeArguments object."); }
          }
      }
          
      *rcut_max = std::max( *rcut_max , parameters->rcut );

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

      auto compute_opt_locks = [&](auto cp_locks)
      {
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
        KimForceOp force_op { kim_ctx->m_thread_ctx.data(), *kim_particle_codes };
        compute_cell_particle_pairs( *grid, parameters->rcut, *ghost, optional, force_buf, force_op , compute_force_field_set , parallel_execution_context() );
      };
      if( omp_get_max_threads() > 1 ) {
        compute_opt_locks( ComputePairOptionalLocks<true>{ particle_locks->data() } );
      } else {
        compute_opt_locks( ComputePairOptionalLocks<false>{} );
      }
      
    }

  };

  template<class GridT> using KIMComputeForceTmpl = KIMComputeForce<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(kim_force)
  {
    OperatorNodeFactory::instance()->register_factory( "kim_force" ,make_grid_variant_operator< KIMComputeForceTmpl > );
  }

}


