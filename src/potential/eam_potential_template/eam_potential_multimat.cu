// _cuda_enable // DO NOT REMOVE THIS LINE

//  // DO NOT REMOVE THIS LINE

#include <exanb/core/grid.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/particle_type_pair.h>
#include <exanb/core/domain.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/cpp_utils.h>

#include <exaStamp/potential/eam/eam_buffer.h>
#include <exaStamp/potential/eam/eam_yaml.h>
#include "potential.h"

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_cell_particles.h>
#include <onika/parallel/memset.h>

#ifdef USTAMP_POTENTIAL_EAM_MM // operator compiled only if potential is multimaterial

#include "eam_force_op_multimat.h"

#ifndef USTAMP_POTENTIAL_EAM_MM_INIT_TYPES
#define USTAMP_POTENTIAL_EAM_MM_INIT_TYPES(p,nt,pe) /**/
#endif

namespace exaStamp
{
  using namespace exanb;

  class EamParameterInitName : public OperatorNode
  {  
    ADD_SLOT( USTAMP_POTENTIAL_PARMS, parameters       , OUTPUT , REQUIRED );
  public:
    inline void execute () override final {}
  };

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz, field::_type >
    >
  class EamPotentialOperatorName : public OperatorNode
  {  
    using CellParticles = typename GridT::CellParticles;

    using EamScratch = EamMultimatPotentialScratch< USTAMP_POTENTIAL_PARMS >;
    using StringVector = std::vector< std::string >;
    template<bool Sym,class CPLocksT> using SymRhoOp = PRIV_NAMESPACE_NAME::SymRhoOp<Sym,CPLocksT>;
    template<bool Sym,class CPLocksT, class VirFieldT> using SymForceOp = PRIV_NAMESPACE_NAME::SymForceOp<Sym,CPLocksT,VirFieldT>;
    using Rho2EmbOp = PRIV_NAMESPACE_NAME::Rho2EmbOp;
    using ForceOpExt = PRIV_NAMESPACE_NAME::ForceOpExtStorage;
    using ForceOpEnergyExt = PRIV_NAMESPACE_NAME::ForceEnergyOpExtStorage;
    using RhoOpExtStorage = PRIV_NAMESPACE_NAME::RhoOpExtStorage;

    // attributes processed during computation
    //static inline constexpr FieldSet< field::_type , field::_virial > cp_force_fields_v{};
    static inline constexpr FieldSet< field::_type > cp_emb_fields_v{};
    static inline constexpr FieldSet< field::_ep , field::_type > cp_emb_fields_energy_v{};

    // ========= I/O slots =======================
    ADD_SLOT( ParticleSpecies       , species          , INPUT , REQUIRED );
    ADD_SLOT( USTAMP_POTENTIAL_PARMS, parameters       , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( StringVector          , types            , INPUT , StringVector{} , DocString{"Empty list means all types are used, otherwise list types handled by this potential"} );
    ADD_SLOT( double                , rcut             , INPUT );
    ADD_SLOT( double                , rcut_max         , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                , ghost_dist_max   , INPUT_OUTPUT , 0.0 );   
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors  , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( GridT                 , grid             , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain           , INPUT , REQUIRED );

    ADD_SLOT( bool                  , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( bool                  , compute_virial   , INPUT , false );

    ADD_SLOT( bool                  , eam_rho      , INPUT , false );
    ADD_SLOT( bool                  , eam_rho2emb  , INPUT , false );
    ADD_SLOT( bool                  , eam_ghost    , INPUT , true );
    ADD_SLOT( bool                  , eam_force    , INPUT , true );
    ADD_SLOT( bool                  , eam_symmetry , INPUT , false );

    ADD_SLOT( GridParticleLocks     , particle_locks      , INPUT_OUTPUT , OPTIONAL , DocString{"particle spin locks"} );

    ADD_SLOT( EamScratch            , eam_scratch      , PRIVATE );
  public:

    // Operator execution
    inline void execute () override final
    {
      //MeamPotential meam( *rcut, *parameters );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( ( *eam_rho || *eam_force ) && *eam_ghost )
      {
        *ghost_dist_max = std::max( *ghost_dist_max , (*rcut) * 2.0 );
      }

      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 ) { return ; } // short cut to avoid errors in pre-initialization step

      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }

      const bool need_particle_locks = ( omp_get_max_threads() > 1 ) && ( *eam_symmetry ) ;
      const bool need_virial = log_energy && *compute_virial;

      ldbg << "EAM Multimat: rho="<<std::boolalpha<< *eam_rho
           <<" , rho2emb="<< *eam_rho 
           <<" , rho2emb="<< *eam_rho2emb 
           <<" , force="<< *eam_force
           <<" , ghost="<< *eam_ghost 
           <<" , sym="<< *eam_symmetry
           <<" , eflag="<< log_energy
           <<" , virflag="<<need_virial
           <<" , need_locks="<< need_particle_locks << std::endl;
      
      // we use eflag also to trigger virial computation
      log_energy = log_energy || need_virial;

      size_t n_species = species->size();
      size_t n_type_pairs = unique_pair_count( n_species );
      
      bool initialize_scratch = eam_scratch->m_pair_enabled.empty();
      if( initialize_scratch )
      {
        eam_scratch->m_pair_enabled.assign( n_type_pairs , false );
        for(size_t i=0;i<n_type_pairs;i++)
        {
          unsigned int a=0, b=0;
          pair_id_to_type_pair(i,a,b);
          const bool a_enabled = types->empty() || ( std::find( types->begin() , types->end() , species->at(a).name() ) != types->end() );
          const bool b_enabled = types->empty() || ( std::find( types->begin() , types->end() , species->at(b).name() ) != types->end() );
          eam_scratch->m_pair_enabled[i] = ( a_enabled && b_enabled );
        }
      }
      USTAMP_POTENTIAL_EAM_MM_INIT_TYPES( *parameters , n_species , eam_scratch->m_pair_enabled.data() );

      auto rho_emb_field = grid->field_accessor( field::rho_dEmb );
      auto * rho_emb_ptr = rho_emb_field.m_flat_array_ptr;

      auto c_rho_emb_field = grid->field_const_accessor( field::rho_dEmb );
      const auto * c_rho_emb_ptr = c_rho_emb_field.m_flat_array_ptr;

      // execute the 2 passes
      auto compute_eam_force = [&]( const auto& cp_xform , const auto& cp_locks , auto cp_emb_fields , auto force_buf )
      {
        using CPLocksT = std::remove_reference_t<decltype(cp_locks)>;
        // common bricks for both compute passes
        exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
        exanb::GridChunkNeighborsLightWeightIt<true> nbh_it_sym{ *chunk_neighbors };
        ComputePairNullWeightIterator cp_weight {};
        const auto * __restrict__ c_particle_offset = grid->cell_particle_offset_data();

        // 1st pass (new) computes rho then emb, without compute buffer
        if( *eam_rho )
        {
          onika::parallel::parallel_memset( rho_emb_ptr , grid->number_of_particles() , 0.0 , parallel_execution_context() );

          auto rho_op_fields = make_field_tuple_from_field_set( FieldSet< field::_type >{} );
          auto rho_buf_factory = ComputePairBufferFactory< ComputePairBuffer2<false,false,RhoOpExtStorage> > {}; // make_empty_pair_buffer<RhoOpExtStorage>();
          if( *eam_symmetry )
          {
            auto rho_optional = make_compute_pair_optional_args( nbh_it_sym, cp_weight , cp_xform , cp_locks, ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{} );
            SymRhoOp<true,CPLocksT> rho_op { *parameters, eam_scratch->m_pair_enabled.data(), c_particle_offset, rho_emb_ptr, cp_locks };
            compute_cell_particle_pairs( *grid, *rcut, *eam_ghost, rho_optional, rho_buf_factory, rho_op, rho_op_fields, parallel_execution_context() );
          }
          else
          {
            auto rho_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks, ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{} );
            SymRhoOp<false,CPLocksT> rho_op { *parameters, eam_scratch->m_pair_enabled.data(), c_particle_offset, rho_emb_ptr, cp_locks };
            compute_cell_particle_pairs( *grid, *rcut, *eam_ghost, rho_optional, rho_buf_factory, rho_op, rho_op_fields, parallel_execution_context() );
          }
        }
        
        if( *eam_rho2emb )
        {
          Rho2EmbOp rho2emb_op { *parameters , eam_scratch->m_pair_enabled.data() };
          auto rho2emb_op_fields = make_field_tuple_from_field_set( cp_emb_fields , rho_emb_field );
          compute_cell_particles( *grid , *eam_ghost , rho2emb_op , rho2emb_op_fields , parallel_execution_context() );
        }

        // 2nd pass parameters: compute final force using the emb term, only for non ghost particles (but reading EMB terms from neighbor particles)
        if( *eam_force )
        {
          auto vir_field = grid->field_accessor(field::virial);
          using VirFieldT = decltype(vir_field);
          auto cp_force_fields_v = onika::make_flat_tuple( grid->field_accessor(field::type) , vir_field );
          
          if( *eam_symmetry )
          {
            auto force_optional = make_compute_pair_optional_args( nbh_it_sym, cp_weight , cp_xform , cp_locks );
            SymForceOp<true,CPLocksT,VirFieldT> force_op { *parameters, eam_scratch->m_pair_enabled.data(), c_particle_offset, c_rho_emb_ptr , cp_locks , vir_field };
            compute_cell_particle_pairs( *grid, *rcut, false, force_optional, force_buf, force_op, cp_force_fields_v, parallel_execution_context() );
          }
          else
          {
            auto force_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks );
            SymForceOp<false,CPLocksT,VirFieldT> force_op { *parameters, eam_scratch->m_pair_enabled.data(), c_particle_offset, c_rho_emb_ptr, cp_locks , vir_field };
            compute_cell_particle_pairs( *grid, *rcut, false, force_optional, force_buf, force_op, cp_force_fields_v, parallel_execution_context() );
          }          
        }
      };

      if( need_particle_locks && ! particle_locks.has_value() )
      {
        fatal_error()<<"particle_locks is needed, but corresponding slot has no value" << std::endl;
      }
      
      auto compute_eam_xform_locks = [&](const auto& cp_xform, const auto& cp_locks)
      {
        if( log_energy ) compute_eam_force( cp_xform, cp_locks, cp_emb_fields_energy_v , ComputePairBufferFactory< ComputePairBuffer2<false,false,ForceOpEnergyExt> >{} );
        else             compute_eam_force( cp_xform, cp_locks, cp_emb_fields_v        , ComputePairBufferFactory< ComputePairBuffer2<false,false,ForceOpExt> >{} );
      };

      if( domain->xform_is_identity() )
      {
        if( need_particle_locks ) compute_eam_xform_locks( exanb::NullXForm{} , exanb::ComputePairOptionalLocks<true>{particle_locks->data()} );
        else                      compute_eam_xform_locks( exanb::NullXForm{} , exanb::ComputePairOptionalLocks<false>{} );
      }
      else
      {
        if( need_particle_locks ) compute_eam_xform_locks( exanb::LinearXForm{domain->xform()} , exanb::ComputePairOptionalLocks<true>{particle_locks->data()} );
        else                      compute_eam_xform_locks( exanb::LinearXForm{domain->xform()} , exanb::ComputePairOptionalLocks<false>{} );
      }

    }

  };

  namespace tmplhelper
  {
    template<class GridT> using EamPotentialOperatorName  = ::exaStamp::EamPotentialOperatorName<GridT>;
  }

  // === register factories ===  
  ONIKA_AUTORUN_INIT(eam_potential_multimat)
  {
    OperatorNodeFactory::instance()->register_factory( EamPotentialStr , make_grid_variant_operator< tmplhelper::EamPotentialOperatorName > );
    OperatorNodeFactory::instance()->register_factory( EamParameterInitStr , make_simple_operator< EamParameterInitName > );
    
  }

}

#endif // only compiled if potential supports multimaterial

