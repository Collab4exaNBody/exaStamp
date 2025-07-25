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

#pragma once

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <onika/log.h>
#include <onika/cpp_utils.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include "meam_parameters_yaml.h"
#include "meam_parameters_stream.h"

#include <iostream>


#ifdef XSTAMP_MEAM_LJ_MULTIMAT

#include "meam_lj_compute_func.h"
#define MeamName1 meam_lj_force
#define MEAM_FORCE_REQUIRED_FIELDS field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type

#else

#include "meam_compute_func.h"
#define MeamName1 meam_force
#define MEAM_FORCE_REQUIRED_FIELDS field::_ep ,field::_fx ,field::_fy ,field::_fz

#endif


#define MeamPotentialOperatorName MeamName1
#define MeamPotentialStr USTAMP_STR(MeamName1)

#define MEAM_REGISTER_INIT() _MEAM_REGISTER_INIT( MEAM_CONSTRUCTOR_FUNC_NAME )
#define MEAM_CONSTRUCTOR_FUNC_NAME USTAMP_CONCAT(MeamPotentialOperatorName,_init)
#define _MEAM_REGISTER_INIT(name) CONSTRUCTOR_ATTRIB void MAKE_UNIQUE_NAME(name,_,__LINE__,ONIKA_CURRENT_PACKAGE_NAME) ()


namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, MEAM_FORCE_REQUIRED_FIELDS >
    >
  class MeamPotentialOperatorName : public OperatorNode
  {
    using CellParticles = typename GridT::CellParticles;

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    static constexpr FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz > compute_field_set = {};

#   ifdef XSTAMP_MEAM_LJ_MULTIMAT
    // functor that populate compute buffer's extended storage's S field with particle types
    using ComputeBufferT = ComputePairBuffer2<false,true, MeamLJExtraStorageT<MeamParameters::MAX_PARTICLE_NEIGHBORS> , MeamLJNeighborFilter , MeamParameters::MAX_PARTICLE_NEIGHBORS >;
    using ForceOp = MeamLJForceComputeFunctor<CellParticles,ComputeBufferT>;
#   else
    using ComputeBufferT = ComputePairBuffer2<false,true, MeamExtraStorageT<MeamParameters::MAX_PARTICLE_NEIGHBORS> , DefaultComputePairBufferAppendFunc , MeamParameters::MAX_PARTICLE_NEIGHBORS >;;
    using ForceOp = MeamForceComputeFunctor<CellParticles,ComputeBufferT>;
#   endif

    struct ComputeMeamScratch
    {
      GridT * gridp = nullptr;
      Mat3d xform;
      ForceOp cp_force;
#   ifdef XSTAMP_MEAM_LJ_MULTIMAT
      std::vector<MeamLJParms> lj_parameters;
#   endif
      bool xform_is_identity = false;
    };


    // ========= I/O slots =======================
    ADD_SLOT( MeamParameters        , parameters        , INPUT        , REQUIRED , DocString{"set of MEAM parameters"} );
    ADD_SLOT( double                , rcut              , INPUT        , REQUIRED , DocString{"base rcut"} );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0      , DocString{"maximum rcut among all potentials"} );
    ADD_SLOT( double                , ghost_dist_max    , INPUT_OUTPUT , 0.0 );   
    ADD_SLOT( bool                  , ghost             , INPUT , true , DocString{"if true, computations are performed in ghost area also."} );
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT , REQUIRED , DocString{"particle grid"} );
    ADD_SLOT( Domain                , domain            , INPUT        , REQUIRED , DocString{"Domain description, including deformation box"} );

#   ifdef XSTAMP_MEAM_LJ_MULTIMAT
    ADD_SLOT( ParticleSpecies       , species           , INPUT        , REQUIRED );
    ADD_SLOT( std::string           , meam_type         , INPUT        , REQUIRED , DocString{"'Sees' only particles of this type for the MEAM part"} );
    ADD_SLOT( UserMeamLJParameters  , lj_parameters     , INPUT        , REQUIRED , DocString{"LJ parameters for non-MEAM types"} );
#   endif 

    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT        , OPTIONAL , DocString{"particle spin locks"} );

    ADD_SLOT( ComputeMeamScratch    , compute_meam_scratch, PRIVATE );
    // ===========================================

  public:
    // -----------------------------------------------
    // ----------- Operator documentation ------------
    inline std::string documentation() const override final
    {
      return R"EOF(
        MEAM (Modified Embedded Atom Model) potential implemntation, with screening.
        No additional mpi communication required.
        2x neighbor distance in ghost area required.
        global parameter ghost_dist_max must be 2x bigger than potential's rcut to have correct results
        )EOF";
    }

    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      // compute maximum rcut needed
      MeamPotential meam( *rcut, *parameters );
      double compute_rcut = meam.p.Rcut;
      if( *ghost )
      {
        *ghost_dist_max = std::max( *ghost_dist_max , meam.p.Rcut * 2.0 );
      }
      else
      {
        // if we don't rely on computation on ghosts,
        // then we don't need correct neighborhood for particles
        // surrounding owned subdomain, thus we don't need ghost_dist_max to be scaled by 2 anymore
        *ghost_dist_max = std::max( *ghost_dist_max , meam.p.Rcut );
      }

#     ifdef XSTAMP_MEAM_LJ_MULTIMAT
      for(const auto & ljp : *lj_parameters)
      {
        compute_rcut = std::max( compute_rcut , ljp.rcut );
        *ghost_dist_max = std::max( *ghost_dist_max , ljp.rcut );
      }
#     endif
      *rcut_max = std::max( compute_rcut , *rcut_max );

      ldbg << "MEAM : rcut="<<meam.p.Rcut<<" , compute_rcut="<<compute_rcut<<", rcut_max="<<(*rcut_max) << std::endl;

      if( grid->number_of_cells() == 0 ) { return ; }

#     ifdef XSTAMP_MEAM_LJ_MULTIMAT
      ldbg << "MEAM/LJ : rcut="<<meam.p.Rcut<<", rcut_max="<< (*rcut_max) << std::endl;
      int ftype = -1;
      auto & spvec = *species;
      auto find_type_id = [&spvec](const std::string& name) -> int
        {
          int ftype = -1;
          for(unsigned int i=0;i<spvec.size();i++)
          {
            if( spvec[i].name() == name ) ftype=i;
          }
          return ftype;
        };
      ftype = find_type_id(*meam_type);
      if( ftype == -1 )
      {
        lerr << "MEAM type '"<< *meam_type << "' not found" << std::endl;
        std::abort();
      }
      
      compute_meam_scratch->cp_force = ForceOp { { std::move(meam) } };
      compute_meam_scratch->cp_force.m_meam_lj_multimat.m_meam_type = ftype;
      compute_meam_scratch->cp_force.m_meam_lj_multimat.m_meam_pair_id = unique_pair_id( ftype, ftype );
      compute_meam_scratch->cp_force.m_meam_lj_multimat.m_meam_rcut2 = meam.p.Rcut * meam.p.Rcut;
      
      if( ! lj_parameters->empty() && compute_meam_scratch->lj_parameters.empty() )
      {
        compute_meam_scratch->lj_parameters.assign( unique_pair_count(species->size()) , MeamLJParms{} );
        for(const auto & ljp : *lj_parameters)
        {
          int ta = find_type_id(ljp.type_a);
          if( ta<0 ) { fatal_error() << "type '"<<ljp.type_a<<"' not found in species"<<std::endl; }
          int tb = find_type_id(ljp.type_b);
          if( tb<0 ) { fatal_error() << "type '"<<ljp.type_b<<"' not found in species"<<std::endl; }
          const unsigned int pair_id = unique_pair_id(ta,tb);
          compute_meam_scratch->lj_parameters[pair_id] = ljp.lj;
          ldbg << "MEAM/LJ: add LJ pair "<<ljp.type_a<<"/"<<ljp.type_b<<" : rcut="<<ljp.rcut<<", rcut2="<<ljp.lj.rcut2<<", ecut="<<ljp.lj.ecut<< std::endl;
        }
      }
      if( compute_meam_scratch->lj_parameters.size() > MeamLJMultiMatParms::MAX_TYPE_PAIRS )
      {
        fatal_error()<<"MEAM/LJ pairs failure : maximum allowed materials overflow ("<<species->size()<<">"<< MeamLJMultiMatParms::MAX_TYPES<<")"<<std::endl;
      }
      for(unsigned int i=0;i<compute_meam_scratch->lj_parameters.size();i++)
      {
        compute_meam_scratch->cp_force.m_meam_lj_multimat.m_lj_parms[i] = compute_meam_scratch->lj_parameters[i];
      }
#     else
      compute_meam_scratch->cp_force = ForceOp { std::move(meam) };
#     endif

      compute_meam_scratch->xform_is_identity = domain->xform_is_identity();
      compute_meam_scratch->xform = domain->xform();

      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      LinearXForm cp_xform{ domain->xform() };
      ComputePairNullWeightIterator cp_weight{};

      bool allow_cuda_exec = ( global_cuda_ctx() != nullptr );
      if( allow_cuda_exec ) allow_cuda_exec = global_cuda_ctx()->has_devices();
      if( allow_cuda_exec )
      {
        ldbg << "MEAM : GPU version" << std::endl;
        ComputePairOptionalLocks<false> cp_locks {}; // no locks needed for GPU version, it uses atomicAdd
        auto cp_optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks );
        auto cp_force_buf = make_compute_pair_buffer< ComputeBufferT >();
        compute_cell_particle_pairs( *grid, compute_rcut, *ghost, cp_optional, cp_force_buf, compute_meam_scratch->cp_force, compute_field_set, parallel_execution_context() );
      }
      else
      {
        ldbg << "MEAM : CPU version" << std::endl;
        if( ! particle_locks.has_value() )
        {
          lerr << "Fatal: no input particle_locks provided" << std::endl;
          std::abort();
        }
        ComputePairNullWeightIterator cp_weight{};
        ComputePairOptionalLocks<true> cp_locks { particle_locks->data() };
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks );
        auto cp_force_buf = make_compute_pair_buffer< ComputeBufferT >();
        compute_cell_particle_pairs( *grid, compute_rcut, *ghost, optional, cp_force_buf, compute_meam_scratch->cp_force, compute_field_set, parallel_execution_context() );
      }
    }

  };

  namespace tplhelper
  {
    template<class GridT> using MeamPotentialOperatorName = ::exaStamp::MeamPotentialOperatorName<GridT>;
  }

  // === register factories ===  
  MEAM_REGISTER_INIT()
  {
    OperatorNodeFactory::instance()->register_factory( MeamPotentialStr , make_grid_variant_operator< tplhelper::MeamPotentialOperatorName > );
  }

}

