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

//  // DO NOT REMOVE THIS LINE
//  // DO NOT REMOVE THIS LINE

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

#include <exaStamp/potential/reaction_field/reaction_field.h>

#include <exanb/core/config.h> // for MAX_PARTICLE_NEIGHBORS constant
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/core/concurent_add_contributions.h>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  template<bool _ComputeEnergy, bool _ComputeVirial>
  struct ReactionFieldComputeContext;

  template<> struct ReactionFieldComputeContext<false,false>
  {
    static inline constexpr bool ComputeEnergy = false;
    static inline constexpr bool ComputeVirial = false;
    double charge_a = 0.0;
    Vec3d f = {0.,0.,0.};
  };

  template<> struct ReactionFieldComputeContext<true,false>
  {
    static inline constexpr bool ComputeEnergy = true;
    static inline constexpr bool ComputeVirial = false;
    double charge_a = 0.0;
    Vec3d f = {0.,0.,0.};
    double ep = 0.0;
  };

  template<> struct ReactionFieldComputeContext<true,true>
  {
    static inline constexpr bool ComputeEnergy = true;
    static inline constexpr bool ComputeVirial = true;
    double charge_a = 0.0;
    Vec3d f = {0.,0.,0.};
    double ep = 0.0;
    Mat3d virial = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  };

  // Reaction Field Compute functor
  template<class CPLocksT, class ChargeFieldT, class TypeFieldT , class VirialFieldT, bool _PerAtomCharge=false, bool _UseSymetry=false, bool _ComputeEnergy = false, bool _ComputeVirial = false>
  struct ReactioFieldForceOp
  {
    static inline constexpr bool PerAtomCharge = _PerAtomCharge;
    static inline constexpr bool ComputeEnergy = _ComputeEnergy;
    static inline constexpr bool ComputeVirial = _ComputeVirial;
    static inline constexpr bool UseSymetry = _UseSymetry;

    static_assert( !ComputeVirial || ComputeEnergy );
    
    // poetential parameters
    const ReactionFieldParms m_params;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    CPLocksT & m_locks;
    ChargeFieldT m_charge_field;
    TypeFieldT m_type_field;
    VirialFieldT m_virial_field;

    using ParticleLockT = decltype( m_locks[0][0] );

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a , size_t p_a, exanb::ComputePairParticleContextStart ) const
    {
      ctx.ext.f = Vec3d{0.,0.,0.};
      if constexpr ( ComputeEnergy )
      {
        ctx.ext.ep = 0.0;
        if constexpr ( ComputeVirial ) ctx.ext.virial = Mat3d{0.,0.,0.,0.,0.,0.,0.,0.,0.};
      }
      if constexpr (  PerAtomCharge ) ctx.ext.charge_a = cells[cell_a][m_charge_field][p_a];
      if constexpr ( !PerAtomCharge ) ctx.ext.charge_a = m_species[ cells[cell_a][field::type][p_a] ].m_charge;
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
    {
      static constexpr bool CPAA = UseSymetry &&   gpu_device_execution();
      static constexpr bool LOCK = UseSymetry && ! gpu_device_execution();

      if constexpr ( ComputeEnergy && ComputeVirial )
      {
        concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
            m_locks[cell_a][p_a]
          , cells[cell_a][field::fx][p_a], cells[cell_a][field::fy][p_a], cells[cell_a][field::fz][p_a], cells[cell_a][field::ep][p_a], cells[cell_a][m_virial_field][p_a]
          , ctx.ext.f.x, ctx.ext.f.y, ctx.ext.f.z, ctx.ext.ep, ctx.ext.virial );
      }
      if constexpr ( ComputeEnergy && !ComputeVirial )
      {
        concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
            m_locks[cell_a][p_a]
          , cells[cell_a][field::fx][p_a], cells[cell_a][field::fy][p_a], cells[cell_a][field::fz][p_a], cells[cell_a][field::ep][p_a]
          , ctx.ext.f.x, ctx.ext.f.y, ctx.ext.f.z, ctx.ext.ep );
      }
      if constexpr ( !ComputeEnergy && !ComputeVirial )
      {
        concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double> (
            m_locks[cell_a][p_a]
          , cells[cell_a][field::fx][p_a], cells[cell_a][field::fy][p_a], cells[cell_a][field::fz][p_a]
          , ctx.ext.f.x, ctx.ext.f.y, ctx.ext.f.z );
      }
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () (
      ComputeBufferT& ctx, Vec3d dr,double d2,
      CellParticlesT cells,size_t cell_b,size_t p_b, double weight ) const
    {
      static constexpr bool CPAA = UseSymetry &&   gpu_device_execution();
      static constexpr bool LOCK = UseSymetry && ! gpu_device_execution();

      double charge_b = 0.0;
      if constexpr (  PerAtomCharge ) charge_b = cells[cell_b][m_charge_field][p_b];
      if constexpr ( !PerAtomCharge ) charge_b = m_species[ cells[cell_b][field::type][p_b] ].m_charge;

      const double r = std::sqrt(d2);
      double e=0.0, de=0.0;
      reaction_field_compute_energy( m_params, ctx.ext.charge_a * charge_b, r, e, de );
      e *= weight; de *= weight; // weighting function
      de /= r;
      const Vec3d dr_fe = de * dr;
      ctx.ext.f += dr_fe;
      [[maybe_unused]] Mat3d virial;
      
      if constexpr ( ComputeEnergy )
      {
        ctx.ext.ep += .5 * e;
        if constexpr ( ComputeVirial )
        {
          virial = tensor( dr_fe, dr ) * -0.5;
          ctx.ext.virial += virial;
        }
      }

      if constexpr ( UseSymetry )
      {
        if constexpr ( ComputeEnergy && ComputeVirial )
        {
          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              m_locks[cell_b][p_b]
            , cells[cell_b][field::fx][p_b], cells[cell_b][field::fy][p_b], cells[cell_b][field::fz][p_b], cells[cell_b][field::ep][p_b], cells[cell_b][m_virial_field][p_b]
            , -dr_fe.x, -dr_fe.y , -dr_fe.z, .5*e, virial );
        }
        if constexpr ( ComputeEnergy && !ComputeVirial )
        {
          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              m_locks[cell_b][p_b]
            , cells[cell_b][field::fx][p_b], cells[cell_b][field::fy][p_b], cells[cell_b][field::fz][p_b], cells[cell_b][field::ep][p_b]
            , -dr_fe.x, -dr_fe.y , -dr_fe.z, .5*e );
        }
        if constexpr ( !ComputeEnergy && !ComputeVirial )
        {
          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double> (
              m_locks[cell_b][p_b]
            , cells[cell_b][field::fx][p_b], cells[cell_b][field::fy][p_b], cells[cell_b][field::fz][p_b]
            , -dr_fe.x, -dr_fe.y , -dr_fe.z );
        }
      }

    }
    
  };

}

namespace exanb
{
  template<class CPLocksT, class ChargeFieldT, class TypeFieldT , class VirialFieldT, bool _PerAtomCharge, bool _UseSymetry, bool _ComputeEnergy, bool _ComputeVirial>
  struct ComputePairTraits< exaStamp::ReactioFieldForceOp<CPLocksT,ChargeFieldT,TypeFieldT,VirialFieldT,_PerAtomCharge,_UseSymetry,_ComputeEnergy,_ComputeVirial> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible      = false;
    static inline constexpr bool BufferLessCompatible         = true;
    static inline constexpr bool HasParticleContextStart      = true;    
    static inline constexpr bool HasParticleContext           = true;
    static inline constexpr bool HasParticleContextStop       = true;
    static inline constexpr bool CudaCompatible               = true;
  };
}


namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class ReactionFieldPC : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( ReactionFieldParms        , parameters          , INPUT , REQUIRED );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors     , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( CompactGridPairWeights    , compact_nbh_weight  , INPUT , OPTIONAL );
    ADD_SLOT( bool                      , enable_pair_weights , INPUT, true );
    ADD_SLOT( bool                      , per_atom_charge     , INPUT, true );
    ADD_SLOT( bool                      , use_symmetry        , INPUT, false );
    ADD_SLOT( bool                      , compute_virial      , INPUT, false );
    ADD_SLOT( bool                      , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( Domain                    , domain              , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies           , species             , INPUT , REQUIRED );    

    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( double                    , rcut_max            , INPUT_OUTPUT , 0.0 );

    ADD_SLOT( GridParticleLocks         , particle_locks      , INPUT_OUTPUT , OPTIONAL , DocString{"particle spin locks"} );

    // ========= Internal types =======================

    // cell particles array type
    using CellParticles = typename GridT::CellParticles;

  public:
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      const double rcut = parameters->rc;
      *rcut_max = std::max( *rcut_max , rcut );
      
      size_t n_cells = grid->number_of_cells();

      // in this case, nothing to compute.
      // this is usefull case where compute_force is called at the very first to initialize rcut_max
      if( n_cells==0 ) return ;

      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }
       
      const bool need_particle_locks = ( omp_get_max_threads() > 1 ) && ( *use_symmetry ) ;
      const bool need_virial = log_energy ; // && *compute_virial;
      const bool pair_weights = compact_nbh_weight.has_value() && ( *enable_pair_weights );
      const bool need_ghost = *use_symmetry;

      if( need_particle_locks && ! particle_locks.has_value() )
      {
        fatal_error() << "missing particle_locks value"<<std::endl;
      }

      ldbg << std::boolalpha
           <<"Reaction field: rc="<< rcut
           <<" , pair_weights="<< pair_weights
           <<" , log_energy="<< log_energy
           <<" , need_virial="<< need_virial
           <<" , use_symmetry="<< *use_symmetry
           <<" , ghost="<< need_ghost
           <<" , per_atom_charge="<< *per_atom_charge
           <<" , need_locks="<< need_particle_locks << std::endl;


      using ChargeFieldT = decltype( grid->field_accessor( field::charge ) );
      using VirialFieldT = decltype( grid->field_accessor( field::virial ) );
      using TypeFieldT   = decltype( grid->field_accessor( field::type   ) );      

      ChargeFieldT charge_field = {};
      TypeFieldT   type_field   = {};
      VirialFieldT virial_field = {};

      if( *per_atom_charge )
      {
        charge_field = grid->field_accessor( field::charge );
      }
      else
      {
        type_field = grid->field_accessor( field::type );
      }      
      if( need_virial )
      {
        virial_field = grid->field_accessor( field::virial );
      }
            
      auto compute_force_energy = [&](auto & cp_locks, const auto & cp_weight, const auto & force_op, auto && cpbuf_factory )
      {
        exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
        LinearXForm cp_xform { domain->xform() };
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
//        [[maybe_unused]] static constexpr onika::parallel::AssertFunctorSizeFitIn< alignof(force_op) , 1 , decltype(force_op) > _check_functor_size = {};
        static constexpr std::true_type use_cells_accessor = {};
        compute_cell_particle_pairs2( *grid, rcut, need_ghost, optional, cpbuf_factory, force_op, onika::FlatTuple<>{}, DefaultPositionFields{}, parallel_execution_context(), use_cells_accessor );
      };
      
      auto compute_force_energy_opt_weights = [&](auto & cp_locks, auto && cp_weight)
      {
        // template<class CPLocksT, class ChargeFieldT, class TypeFieldT , class VirialFieldT, bool _PerAtomCharge=false, bool _UseSymetry=false, bool _ComputeEnergy = false, bool _ComputeVirial = false>
        if( log_energy )
        {
          using ForceOp = ReactioFieldForceOp<decltype(cp_locks),ChargeFieldT,TypeFieldT,VirialFieldT,true,false,true,true>;
          using CPBufT = ComputePairBuffer2<false,false, ReactionFieldComputeContext<true,true> >;
          compute_force_energy( cp_locks, cp_weight, ForceOp{*parameters, species->data(), cp_locks, charge_field, type_field, virial_field} , make_compute_pair_buffer<CPBufT>() );
        }
        else
        {
          using ForceOp = ReactioFieldForceOp<decltype(cp_locks),ChargeFieldT,TypeFieldT,VirialFieldT,true,false,false,false>;
          using CPBufT = ComputePairBuffer2<false,false, ReactionFieldComputeContext<false,false> >;
          compute_force_energy( cp_locks, cp_weight, ForceOp{*parameters, species->data(), cp_locks, charge_field, type_field, virial_field} , make_compute_pair_buffer<CPBufT>() );
        }
      };
      
      auto compute_force_energy_opt_locks = [&](auto && cp_locks)
      {
        if( pair_weights ) compute_force_energy_opt_weights( cp_locks, CompactPairWeightIterator{ compact_nbh_weight->m_cell_weights.data() } );
        else               compute_force_energy_opt_weights( cp_locks, ComputePairNullWeightIterator{} );
      };

      if( need_particle_locks ) compute_force_energy_opt_locks( ComputePairOptionalLocks<false>{} );
      else                      compute_force_energy_opt_locks( ComputePairOptionalLocks<true>{ particle_locks->data() } );      
    }

  };

  template<class GridT> using ReactionFieldPCTmpl = ReactionFieldPC<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(reaction_field)
  {  
    OperatorNodeFactory::instance()->register_factory( "reaction_field" , make_grid_variant_operator<ReactionFieldPCTmpl> );
  }

}


