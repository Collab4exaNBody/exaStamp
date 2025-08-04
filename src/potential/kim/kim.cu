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

#include <onika/physics/constants.h>
#include <onika/cuda/cuda_context.h>

#include <vector>
#include <memory>
#include <iostream>
#define MY_ERROR(message)                                                \
  {                                                                      \
    std::cout << "* Error : \"" << message << "\" : " << __LINE__ << ":" \
              << __FILE__ << std::endl;                                  \
    exit(1);                                                             \
  }

namespace exaStamp
{
  using onika::memory::DEFAULT_ALIGNMENT;


  /* Define neighborlist structure */
  typedef struct
  {
    double cutoff;
    int numberOfParticles;
    int * NNeighbors;
    int * neighborList;
  } NeighList;
  
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_id >
    >
  class KIMComputeForce : public OperatorNode
  {
    // ========= I/O slots =======================
    //    ADD_SLOT( KIMParams             , parameters        , INPUT        , REQUIRED );
    ADD_SLOT( std::string           , kim_model_name    , INPUT );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0      );
    ADD_SLOT( ParticleSpecies       , species           , INPUT        , REQUIRED );
    ADD_SLOT( int64_t               , timestep          , INPUT        , REQUIRED );
    ADD_SLOT( GridChunkNeighbors    , chunk_neighbors   , INPUT        , GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT        , false    );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT            );
    ADD_SLOT( Domain                , domain            , INPUT        , REQUIRED );
    ADD_SLOT( KIMContext            , kim_ctx           , INPUT );

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    using ComputeBuffer = ComputePairBuffer2<false,false>;
    using CellParticles = typename GridT::CellParticles;
    //    using ParticleLock  = decltype( ComputePairOptionalLocks<false>{}[0][0] );

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_id >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_id, field::_virial >;
    using ComputeFields              = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};
    
  public:
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      size_t nt = omp_get_max_threads();
      if (nt > kim_ctx->m_thread_ctx.size()) {
        size_t old_nt = kim_ctx->m_thread_ctx.size();
        std::cout << "resizing thread context " << std::endl;
        std::cout << "\told size = " << old_nt << ", new size = " << nt << std::endl;
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
                                       *kim_model_name,
                                       &requestedUnitsAccepted,
                                       &kim_ctx->m_thread_ctx[j].kim_model);
            if (error) { MY_ERROR("KIM::Model::Create()"); }
            // Check for compatibility with the model
            if (!requestedUnitsAccepted) { MY_ERROR("Must Adapt to model units"); }


            // Check that we know about all required routines
            int numberOfModelRoutineNames;
            KIM::MODEL_ROUTINE_NAME::GetNumberOfModelRoutineNames(&numberOfModelRoutineNames);
      
            for (int i = 0; i < numberOfModelRoutineNames; ++i)
              {
                KIM::ModelRoutineName modelRoutineName;
                int error = KIM::MODEL_ROUTINE_NAME::GetModelRoutineName(i, &modelRoutineName);
                if (error) { MY_ERROR("Unable to get ModelRoutineName."); }
                int present;
                int required;
                error = kim_ctx->m_thread_ctx[j].kim_model->IsRoutinePresent(modelRoutineName, &present, &required);
                if (error) { MY_ERROR("Unable to get routine present/required."); }

                std::cout << "Model routine name \"" << modelRoutineName.ToString()
                          << "\" has present = " << present
                          << " and required = " << required << "." << std::endl;

                if ((present == true) && (required == true))
                  {
                    using namespace KIM::MODEL_ROUTINE_NAME;
                    if (!((modelRoutineName == Create)
                          || (modelRoutineName == ComputeArgumentsCreate)
                          || (modelRoutineName == Compute) || (modelRoutineName == Refresh)
                          || (modelRoutineName == ComputeArgumentsDestroy)
                          || (modelRoutineName == Destroy)))
                      {
                        MY_ERROR("Unknown Routine \"" + modelRoutineName.ToString()
                                 + "\" is required by model.");
                      }
                  }
              }

            // print model units
            KIM::LengthUnit lengthUnit;
            KIM::EnergyUnit energyUnit;
            KIM::ChargeUnit chargeUnit;
            KIM::TemperatureUnit temperatureUnit;
            KIM::TimeUnit timeUnit;

            kim_ctx->m_thread_ctx[j].kim_model->GetUnits(&lengthUnit, &energyUnit, &chargeUnit, &temperatureUnit, &timeUnit);

            std::cout << "LengthUnit is \"" << lengthUnit.ToString() << "\"" << std::endl
                      << "EnergyUnit is \"" << energyUnit.ToString() << "\"" << std::endl
                      << "ChargeUnit is \"" << chargeUnit.ToString() << "\"" << std::endl
                      << "TemperatureUnit is \"" << temperatureUnit.ToString() << "\""
                      << std::endl
                      << "TimeUnit is \"" << timeUnit.ToString() << "\"" << std::endl;

            // check species
            int speciesIsSupported;
            int modelTaCode;
            error = kim_ctx->m_thread_ctx[j].kim_model->GetSpeciesSupportAndCode(KIM::SPECIES_NAME::Ta, &speciesIsSupported, &modelTaCode);
            if ((error) || (!speciesIsSupported))
              {
                MY_ERROR("Species Ta not supported");
              }

            KIM::ComputeArguments * computeArguments;
            error = kim_ctx->m_thread_ctx[j].kim_model->ComputeArgumentsCreate(&computeArguments);
            if (error) { MY_ERROR("Unable to create a ComputeArguments object."); }
          }
      }
          
      //      if (nt > pace
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
		
      ForceOp force_op { *rcut_max , kim_ctx->m_thread_ctx.data() };
      ComputePairNullWeightIterator          cp_weight{};
      ComputePairOptionalLocks<false>        cp_locks {};
      GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();
      ComputePairTrivialCellFiltering cpu_cell_filter = {};

      if( domain->xform_is_identity() )
        {
          NullXForm cp_xform;
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks, cpu_cell_filter );
          compute_cell_particle_pairs( *grid, *rcut_max, *ghost, optional, force_buf, force_op , compute_force_field_set, parallel_execution_context() );
        }
      else
        {
          LinearXForm cp_xform { domain->xform() };
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks, cpu_cell_filter );
          compute_cell_particle_pairs( *grid, *rcut_max, *ghost, optional, force_buf, force_op , compute_force_field_set, parallel_execution_context());
        }
    }
    
    private:
    
    struct alignas(DEFAULT_ALIGNMENT) ForceOp 
    {
      const double m_rcut;

      KIMThreadContext* m_thread_ctx = nullptr;
      
      inline void operator ()
        (
        size_t n,
        ComputeBuffer& buf,
        double& en,
        double& fx,
        double& fy,
        double& fz,
        unsigned int type, // On a besoin du type de l'atome courant
        unsigned int id, // idem pour l'identifiant de l'atome courant
        CellParticles* unused
        ) const
      {
        Mat3d virial;
        this->operator () ( n,buf,en,fx,fy,fz,type,id,virial, unused );
      }

      inline void operator ()
        (
        size_t n,
        ComputeBuffer& buf,
        double& en,
        double& fx,
        double& fy,
        double& fz,
        unsigned int type,
        unsigned int id,
        Mat3d& virial ,
        CellParticles*
        ) const
      {
        
        size_t tid = omp_get_thread_num();
        assert(tid < m_thread_ctx.size());
        KIMThreadContext & kim_ctx = m_thread_ctx[tid];
        auto kimptr = kim_ctx.kim_model;
        
        // energy and force contributions to the particle
        double _fx = 0.;	
        double _fy = 0.;
        double _fz = 0.;
        int numberOfParticles_cluster = n;
        double energy_cluster_model;
        double forces_cluster[n * 3];

        /* Setup local neighborhood coordinates */
        double coords_cluster[n][3];
        for (int i = 0; i < n; ++i) {
          coords_cluster[i][0] = buf.drx[i];
          coords_cluster[i][1] = buf.dry[i];
          coords_cluster[i][2] = buf.drz[i];
        }
        
        // Flag to decide whether neighboring particle contributes to central particle's
        int particleContributing_cluster_model[n];
        for (int i = 0; i < n; ++i)
          particleContributing_cluster_model[i] = 1; /* every particle contributes */

        /* setup particleSpecies */
        int particleSpecies_cluster_model[n];
        int isSpeciesSupported;
        int error = kimptr->GetSpeciesSupportAndCode(KIM::SPECIES_NAME::Ta,
                                                     &isSpeciesSupported,
                                                     &(particleSpecies_cluster_model[0]));
        if (error) MY_ERROR("get_species_code");
        for (int i = 1; i < n; ++i)
          particleSpecies_cluster_model[i] = particleSpecies_cluster_model[0];
        std::vector<Vec3d> centerParticleCoordinates;
        centerParticleCoordinates.resize(n);

        //        MatrixXd localForces(subconfigOfParticle.numberOfParticles, DIM);
        // localForces.setZero();
        KIM::ComputeArguments * computeArguments;
        error = kimptr->ComputeArgumentsCreate(&computeArguments);

        error = kimptr->Compute(computeArguments);

        int numberOfComputeArgumentNames;
        KIM::COMPUTE_ARGUMENT_NAME::GetNumberOfComputeArgumentNames(&numberOfComputeArgumentNames);
        for (int i = 0; i < numberOfComputeArgumentNames; ++i)
          {
            KIM::ComputeArgumentName computeArgumentName;
            KIM::SupportStatus supportStatus;
            KIM::COMPUTE_ARGUMENT_NAME::GetComputeArgumentName(i, &computeArgumentName);
            KIM::DataType dataType;
            KIM::COMPUTE_ARGUMENT_NAME::GetComputeArgumentDataType(computeArgumentName,
                                                                   &dataType);
            error = computeArguments->GetArgumentSupportStatus(computeArgumentName,
                                                               &supportStatus);
            if (error) MY_ERROR("unable to get ComputeArgument SupportStatus");

            std::cout << "ComputeArgument Name \"" << computeArgumentName.ToString()
                      << "\""
                      << " is of type \"" << dataType.ToString() << "\""
                      << " and has supportStatus \"" << supportStatus.ToString() << "\""
                      << std::endl;
          }
        
        int numberOfParameters;
        kimptr->GetNumberOfParameters(&numberOfParameters);
        for (int i = 0; i < numberOfParameters; ++i)
          {
            KIM::DataType dataType;
            std::string const * strName;
            std::string const * strDesc;
            int extent;
            kimptr->GetParameterMetadata(i, &dataType, &extent, &strName, &strDesc);
            std::cout << "Parameter No. " << i << " has" << std::endl
                      << " data type   : \"" << dataType.ToString() << "\"" << std::endl
                      << " extent      : " << extent << std::endl
                      << " name        : " << *strName << std::endl
                      << " description : " << *strDesc << std::endl;
          }

        // Check supported extensions, if any
        int present;
        error = kimptr->IsRoutinePresent(KIM::MODEL_ROUTINE_NAME::Extension, &present, NULL);
        if (error) { MY_ERROR("Unable to get Extension present/required."); }
        if (present)
          {
            KIM::SupportedExtensions supportedExtensions;
            error = kimptr->Extension(KIM_SUPPORTED_EXTENSIONS_ID,&supportedExtensions);
            if (error) { MY_ERROR("Error returned from KIM::Model::Extension()."); }
            std::cout << "Model Supports "
                      << supportedExtensions.numberOfSupportedExtensions
                      << " Extensions:" << std::endl;
            for (int i = 0; i < supportedExtensions.numberOfSupportedExtensions; ++i)
              {
                std::cout << " spportedExtensionID[" << std::setw(2) << i << "] = \""
                          << supportedExtensions.supportedExtensionID[i] << "\" "
                          << "which has required = "
                          << supportedExtensions.supportedExtensionRequired[i] << "."
                          << std::endl;
              }
          }

        error = computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles,
                                                     (int *) &numberOfParticles_cluster)
          || computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,
                                                  particleSpecies_cluster_model)
          || computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::particleContributing,
                                                  particleContributing_cluster_model)
          || computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::coordinates,
                                                  (double *) coords_cluster)
          || computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
                                                  &energy_cluster_model)
          || computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialForces,
                                                  (double *) forces_cluster);

        // if (error) MY_ERROR("KIM_API_set_data");
        // error = computeArguments->SetCallbackPointer(KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
        //                                              KIM::LANGUAGE_NAME::cpp,
        //                                              (KIM::Function *) &get_cluster_neigh,
        //                                              &nl_cluster_model);
        // if (error) MY_ERROR("set_call_back");
        
        error = kimptr->Compute(computeArguments);
        if (error) MY_ERROR("compute");
        
        error = kimptr->ComputeArgumentsDestroy(&computeArguments);
        
                // // broadcast to model
                // p_kimLocal->broadcastToModel(&subconfigOfParticle,
                //                      subconfigOfParticle.particleContributing,
                //                      &localForces,
                //                      nlOfParticle,
                //                      (KIM::Function *) &nbl_get_neigh,
                //                      nullptr,
                //                      nullptr);
                // // compute partial forces
                // p_kimLocal->compute();
        
      }
    };

  };

  template<class GridT> using KIMComputeForceTmpl = KIMComputeForce<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(kim)
  {
    OperatorNodeFactory::instance()->register_factory( "kim_force" ,make_grid_variant_operator< KIMComputeForceTmpl > );
  }

}


