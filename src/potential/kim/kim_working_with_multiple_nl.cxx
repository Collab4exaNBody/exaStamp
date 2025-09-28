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
//#include "Eigen/Eigen/Dense"
#include <onika/physics/units.h>

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
        
        // FILL BELOW with the appropriate code to compute forces and energy on central particles.
        // number of particles in this local cluster: 1 (center) + n neighbors
        const int np = static_cast<int>(n) + 1;

        // coordinates layout: [x0,y0,z0, x1,y1,z1, ...]
        // put central at origin; neighbors are already r_ij = (drx, dry, drz)
        std::vector<double> coords(3 * static_cast<size_t>(np), 0.0);
        for (int i = 0; i < static_cast<int>(n); ++i) {
          coords[3 * (i + 1) + 0] = buf.drx[i];
          coords[3 * (i + 1) + 1] = buf.dry[i];
          coords[3 * (i + 1) + 2] = buf.drz[i];
        }

        int particleSpecies[np];
        int particleContributing[np];
        
        for (int i = 0; i < np; ++i)
          particleContributing[i] = 1;

        /* setup particleSpecies */
        int isSpeciesSupported;
        int error = kimptr->GetSpeciesSupportAndCode(KIM::SPECIES_NAME::Ta,
                                                     &isSpeciesSupported,
                                                     &(particleSpecies[0]));
        if (error) MY_ERROR("get_species_code");
        for (int i = 1; i < np; ++i)
          particleSpecies[i] = particleSpecies[0];
        
        // species: reuse the code you already queried into particleSpecies_cluster_model[0]
        std::vector<int> species_codes(np, particleSpecies[0]);

        // contributing flags: only the central needs to contribute here
        std::vector<int> contributing(np, 0);
        contributing[0] = 1;

        // outputs
        double localEnergy = 0.0;
        std::vector<double> forces(3 * static_cast<size_t>(np), 0.0);

        // lightweight neighbor list object for the callback
        struct LocalNeighList {
          double cutoff = 0.0;
          int    numberOfParticles = 0;
          std::vector<int> NNeighbors;    // size = np
          std::vector<int> neighborList;  // size = np*np, row-major; row i starts at i*np
        } nl;

        nl.numberOfParticles = np;
        nl.NNeighbors.assign(np, 0);
        nl.neighborList.assign(np * np, 0);

        // provide neighbors for the central particle (index 0): 1..n
        nl.NNeighbors[0] = static_cast<int>(n);
        for (int j = 0; j < static_cast<int>(n); ++j) {
          nl.neighborList[0 * np + j] = j + 1;  // central's j-th neighbor is particle (j+1)
        }

        // prepare ComputeArguments
        KIM::ComputeArguments* computeArguments = nullptr;
        {
          int err = kimptr->ComputeArgumentsCreate(&computeArguments);
          if (err) { MY_ERROR("KIM::ComputeArgumentsCreate() failed."); }
        }

        // (optional) read model neighbor-list hints to set an appropriate cutoff
        int num_nlists = 0;
        double const* cutoffs = nullptr;
        int const* will_not_request_noncontrib = nullptr;
        kimptr->GetNeighborListPointers(&num_nlists, &cutoffs, &will_not_request_noncontrib);
        if (num_nlists != 1) { kimptr->ComputeArgumentsDestroy(&computeArguments); MY_ERROR("Unexpected number of neighbor lists."); }
        nl.cutoff = (cutoffs != nullptr) ? cutoffs[0] : m_rcut;

        // wire required argument pointers
        {
          int np_local = np;
          int err =
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &np_local) ||
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes, species_codes.data()) ||
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::particleContributing, contributing.data()) ||
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::coordinates, coords.data()) ||
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, &localEnergy) ||
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::partialForces, forces.data());
          if (err) { kimptr->ComputeArgumentsDestroy(&computeArguments); MY_ERROR("Error in SetArgumentPointer."); }
        }

        // install a GetNeighborList callback (non-capturing lambda decays to function ptr)
        using GetNeighSig = int(void*, int, double const*, int, int, int*, int const**);
        GetNeighSig* get_neigh_cb = +[](void* dataObject,
                                        int numberOfNeighborLists,
                                        double const* /*cutoffs*/,
                                        int neighborListIndex,
                                        int particleNumber,
                                        int* numberOfNeighbors,
                                        int const** neighborsOfParticle) -> int
        {
          auto* nl = static_cast<LocalNeighList*>(dataObject);
          if (numberOfNeighborLists != 1) return 1;
          if (neighborListIndex != 0)     return 1;
          if (particleNumber < 0 || particleNumber >= nl->numberOfParticles) return 1;

          *numberOfNeighbors   = nl->NNeighbors[particleNumber];
          *neighborsOfParticle = &(nl->neighborList[particleNumber * nl->numberOfParticles]);
          return 0;
        };

        {
          int err = computeArguments->SetCallbackPointer(
                                                         KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
                                                         KIM::LANGUAGE_NAME::cpp,
                                                         reinterpret_cast<KIM::Function*>(get_neigh_cb),
                                                         &nl);
          if (err) { kimptr->ComputeArgumentsDestroy(&computeArguments); MY_ERROR("Error in SetCallbackPointer(GetNeighborList)."); }
        }

        // call the model (use the 1-argument overload as in the official example)
        {
          int err = kimptr->Compute(computeArguments);
          if (err) { kimptr->ComputeArgumentsDestroy(&computeArguments); MY_ERROR("KIM::Model::Compute() failed."); }
        }

        double conv_energy_factor = ONIKA_CONST_QUANTITY( 1. * eV ).convert();
        
        // accumulate central-particle results (index 0)
        for (int i = 1; i < np; ++i)
          {
            fx += ( forces[3*i+0] * conv_energy_factor );
            fy += ( forces[3*i+1] * conv_energy_factor );
            fz += ( forces[3*i+2] * conv_energy_factor );
          }            
        en = (localEnergy * conv_energy_factor);

        // cleanup
        {
          int err = kimptr->ComputeArgumentsDestroy(&computeArguments);
          if (err) { MY_ERROR("KIM::ComputeArgumentsDestroy() failed."); }
        }
        // -----------------------------------------------------------------------------

        //
        
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


