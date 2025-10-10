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
#include "KIM_Log.hpp"
#include "KIM_LogVerbosity.hpp"
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
      KIM::Log::PushDefaultVerbosity(KIM::LOG_VERBOSITY::silent);
      
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

                ldbg << "Model routine name \"" << modelRoutineName.ToString()
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

            ldbg << "LengthUnit is \"" << lengthUnit.ToString() << "\"" << std::endl
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
          
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
		
      ForceOp force_op { *rcut_max , kim_ctx->m_thread_ctx.data() };
      ComputePairNullWeightIterator          cp_weight{};
      ComputePairOptionalLocks<false>        cp_locks {};
      GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      ComputePairTrivialCellFiltering        cpu_cell_filter = {};
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();

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
      KIM::Log::PopDefaultVerbosity(); 
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
        unsigned int type,
        unsigned int id,
        CellParticles* unused
        ) const
      {
        Mat3d virial;
        this->operator () ( n,buf,en,fx,fy,fz,type,id,virial, unused);
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
        //        assert(tid < (*m_thread_ctx).size());
        KIMThreadContext & kim_ctx = m_thread_ctx[tid];
        auto kimptr = kim_ctx.kim_model;
        
        // number of particles in this local cluster: 1 (center) + n neighbors
        const int np = static_cast<int>(n) + 1;
        
        // put central at origin; neighbors are already r_ij = (drx, dry, drz)
        std::vector<double> coords(3 * static_cast<size_t>(np), 0.0);
        for (int i = 0; i < static_cast<int>(n); ++i) {
          coords[3 * ( i + 1 ) + 0] = buf.drx[i];
          coords[3 * ( i + 1 ) + 1] = buf.dry[i];
          coords[3 * ( i + 1 ) + 2] = buf.drz[i];
          //          std::cout << "cx,cy,cz = " << coords[3 * i + 0] << ","<< coords[3 * i + 1] << ","<< coords[3 * i + 2] << std::endl;
        }
        
        // contributing: only central particle is a contributing particle. Other particles just serve to compute energy and force on central particle.
        std::vector<int> contributing(np, 0);
        contributing[0] = 1;
        
        // species: reuse the code you already queried into particleSpecies_cluster_model[0]
        int isSpeciesSupported;
        std::vector<int> species_codes(np, 0);        
        int error = kimptr->GetSpeciesSupportAndCode(KIM::SPECIES_NAME::Ta,
                                                     &isSpeciesSupported,
                                                     &(species_codes[0]));
        if (error) MY_ERROR("get_species_code");

        // Defining outputs
        double localEnergy = 0.0;
        std::vector<double> energies(static_cast<size_t>(np), 0.0);
        std::vector<double> forces(3 * static_cast<size_t>(np), 0.0);
        std::vector<double> virials(6 * static_cast<size_t>(np), 0.0);

        // lightweight neighbor-list payload: ONLY central has neighbors
        struct CentralOnlyNL {
          int np;
          int n;                              // number of neighbors of central
          std::vector<int> neighbors_indices; // length n, values 1..n
        } nl;
        
        nl.np = np;
        nl.n  = static_cast<int>(n);
        nl.neighbors_indices.resize(nl.n);
        for (int j = 0; j < nl.n; ++j) nl.neighbors_indices[j] = j + 1; // neighbors of central

        // prepare ComputeArguments
        KIM::ComputeArguments* computeArguments = nullptr;
        {
          int err = kimptr->ComputeArgumentsCreate(&computeArguments);
          if (err) { MY_ERROR("KIM::ComputeArgumentsCreate() failed."); }
        }

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
                                                 KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, energies.data());
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::partialForces, forces.data());
            computeArguments->SetArgumentPointer(
                                                 KIM::COMPUTE_ARGUMENT_NAME::partialVirial, virials.data());          
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
          //          auto* nl = static_cast<LocalNeighList*>(dataObject);
          auto* nl = static_cast<CentralOnlyNL*>(dataObject);          
          if (numberOfNeighborLists != 1) return 1;
          if (neighborListIndex != 0)     return 1;
          
          // if (particleNumber < 0 || particleNumber >= nl->numberOfParticles) return 1;

          // *numberOfNeighbors   = nl->NNeighbors[particleNumber];
          // *neighborsOfParticle = &(nl->neighborList[particleNumber * nl->numberOfParticles]);

          if (particleNumber == 0) {
            *numberOfNeighbors   = nl->n;
            *neighborsOfParticle = nl->neighbors_indices.data();
          } else {
            // No neighbor list for non-central particles
            *numberOfNeighbors   = 0;
            *neighborsOfParticle = nullptr;
          }
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

        //std::cout << "###################################" << std::endl;
        //std::cout << "list of forces = " << std::endl;
        //std::cout << "fx,fy,fz central atom = " << Vec3d{forces[0],forces[1],forces[2]} << std::endl;
        // Vec3d sumforces {0.,0.,0.};
        // for (int i = 1; i < static_cast<int>(n+1); ++i) {
        //   sumforces += Vec3d{forces[3*i+0],forces[3*i+1],forces[3*i+2] };
        //   //          //std::cout << "fx,fy,fz for atom " << i << " = " << Vec3d{forces[3*i+0],forces[3*i+1],forces[3*i+2] } << std::endl;
        // }
        //std::cout << "sumforces on neihgors = " << sumforces << std::endl;
        //        ldbg << "###################################" << std::endl;

        //        std::cout << "###################################" << std::endl;
        //        std::cout << "list of energies = " << std::endl;
        //        std::cout << "fx,fy,fz central atom = " << Vec3d{forces[0],forces[1],forces[2]} << std::endl;
        //        std::cout << "energy central atom = " << energies[0] << std::endl;        
        // for (int i = 1; i < static_cast<int>(n+1); ++i) {
        //   std::cout << "energy for atom " << i << " = " << energies[i] << std::endl;
        // }
        //std::cout << "sumforces on neihgors = " << sumforces << std::endl;
        //        ldbg << "###################################" << std::endl;                
        double conv_energy_factor = ONIKA_CONST_QUANTITY( 1. * eV ).convert();
        Vec3d localForce = Vec3d{forces[0],forces[1],forces[2]} * conv_energy_factor * 2.0;
        fx += localForce.x;
        fy += localForce.y;
        fz += localForce.z;
        //        en = localEnenergies[0] * conv_energy_factor;
        en = (localEnergy * conv_energy_factor);

        // cleanup
        {
          int err = kimptr->ComputeArgumentsDestroy(&computeArguments);
          if (err) { MY_ERROR("KIM::ComputeArgumentsDestroy() failed."); }
        }
        // -----------------------------------------------------------------------------

        //        KIM::Log::PopDefaultVerbosity();
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


