//  // DO NOT REMOVE THIS LINE

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
#include <onika/cpp_utils.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <onika/thread.h>

#include <exaStamp/potential/eam/eam_buffer.h>
#include <exaStamp/potential/eam/eam_yaml.h>
#include "potential.h"

#include <exanb/particle_neighbors/flat_neighbor_lists.h>
#include <onika/parallel/parallel_for.h>
#include <onika/parallel/memset.h>

#ifdef USTAMP_POTENTIAL_EAM_MM // operator compiled only if potential is multimaterial

#include "eam_force_op_multimat_flat.h"

namespace exaStamp
{
  using namespace exanb;

  template< class GridT >
  class EamPotentialFlatName : public OperatorNode
  {  
    using CellParticles = typename GridT::CellParticles;

    using EamScratch = EamMultimatPotentialScratch< USTAMP_POTENTIAL_PARMS >;
    using StringVector = std::vector< std::string >;

    template<bool NewtonSym, class XFormT> using FlatSymRhoOp = PRIV_NAMESPACE_NAME::FlatSymRhoOp<NewtonSym,XFormT>;
    template<bool EnergyFlag> using FlatRho2EmbOp = PRIV_NAMESPACE_NAME::FlatRho2EmbOp<EnergyFlag>;
    template<bool NewtonSym, class XFormT> using FlatSymForceOp = PRIV_NAMESPACE_NAME::FlatSymForceOp<NewtonSym,XFormT>;

    // ========= I/O slots =======================
    ADD_SLOT( ParticleSpecies       , species          , INPUT , REQUIRED );
    ADD_SLOT( USTAMP_POTENTIAL_PARMS, parameters       , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( StringVector          , types            , INPUT , StringVector{} , DocString{"Empty list means all types are used, otherwise list types handled by this potential"} );
    ADD_SLOT( double                , rcut             , INPUT );
    ADD_SLOT( double                , rcut_max         , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                , ghost_dist_max   , INPUT_OUTPUT , 0.0 );   
    ADD_SLOT( FlatPartNbhList       , flat_nbh_list    , INPUT_OUTPUT );
    ADD_SLOT( GridT                 , grid             , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain           , INPUT , REQUIRED );

    ADD_SLOT( bool                  , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( bool                  , compute_virial   , INPUT , false );

    ADD_SLOT( bool                  , eam_rho      , INPUT , true );
    ADD_SLOT( bool                  , eam_rho2emb  , INPUT , true );
    ADD_SLOT( bool                  , eam_ghost    , INPUT , true );
    ADD_SLOT( bool                  , eam_force    , INPUT , true );
    ADD_SLOT( bool                  , eam_symmetry , INPUT , false );

    ADD_SLOT( spin_mutex_array      , flat_particle_locks       , INPUT_OUTPUT , OPTIONAL , DocString{"particle spin locks"} );

    ADD_SLOT( EamScratch            , eam_scratch          , PRIVATE );
    
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

      const size_t n_particles = grid->number_of_particles();
      if( n_particles == 0 ) { return ; } // short cut to avoid errors in pre-initialization step

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

      ldbg << "EAM Multimat Flat:"
           <<  " rho="<<std::boolalpha<< *eam_rho
           <<" , rho2emb="<< *eam_rho 
           <<" , rho2emb="<< *eam_rho2emb 
           <<" , force="<< *eam_force
           <<" , ghost="<< *eam_ghost 
           <<" , sym="<< *eam_symmetry
           <<" , eflag="<< log_energy
           <<" , virflag="<<need_virial
           <<" , need_locks="<< need_particle_locks << std::endl;
      
      if( need_particle_locks )
      {
        if( ! flat_particle_locks.has_value() )
        {
          fatal_error()<<"flat_particle_locks is needed, but corresponding slot has no value" << std::endl;
        }
        if( flat_particle_locks->size() != grid->number_of_particles() )
        {
          fatal_error()<<"flat_particle_locks has wrong size : "<<flat_particle_locks->size()<<" <> "<< grid->number_of_particles() << std::endl;
        }
      }

      const size_t n_species = species->size();
      const size_t n_type_pairs = unique_pair_count( n_species );
      const bool initialize_scratch = eam_scratch->m_pair_enabled.empty();
      
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

      auto rho_emb = grid->field_accessor( field::rho_dEmb );
      auto energy = grid->field_accessor( field::flat_ep );
      auto fx = grid->field_accessor( field::flat_fx );
      auto fy = grid->field_accessor( field::flat_fy );
      auto fz = grid->field_accessor( field::flat_fz );
      auto rx = grid->field_accessor( field::flat_rx );
      auto ry = grid->field_accessor( field::flat_ry );
      auto rz = grid->field_accessor( field::flat_rz );
      auto types = grid->field_accessor( field::flat_type );

      const double rc = *rcut;
      const double rc2 = rc*rc;

      // execute the 2 passes
      auto compute_eam_force = [&]( const auto& cp_xform , auto newtonSym )
      {
        using XFormT = std::remove_reference_t<decltype(cp_xform)>;
        
        if( *eam_rho )
        {          
          onika::parallel::parallel_memset( rho_emb.m_flat_array_ptr , n_particles , 0.0 , parallel_execution_context() );

          FlatSymRhoOp< newtonSym.value , XFormT > rho_op =
            { *parameters, rc2
            , flat_nbh_list->m_neighbor_offset.data(), flat_nbh_list->m_neighbor_list.data(), flat_nbh_list->m_half_count.data()
            , (*eam_ghost) ? nullptr : grid->particle_ghost_flag_data() , eam_scratch->m_pair_enabled.data()
            , rho_emb.m_flat_array_ptr, types.m_flat_array_ptr, rx.m_flat_array_ptr, ry.m_flat_array_ptr, rz.m_flat_array_ptr
            , cp_xform , flat_particle_locks->data() };

          parallel_for( n_particles , rho_op , parallel_execution_context() );
        }
        
        if( *eam_rho2emb )
        {
          if( log_energy )
          {
            std::cout << "Log energy true" << std::endl;
            FlatRho2EmbOp<true> rho2emb_op{ *parameters, (*eam_ghost) ? nullptr : grid->particle_ghost_flag_data(), eam_scratch->m_pair_enabled.data()
                                          , types.m_flat_array_ptr, rho_emb.m_flat_array_ptr, energy.m_flat_array_ptr };
            parallel_for( n_particles , rho2emb_op , parallel_execution_context() );
          }
          else
          {
            FlatRho2EmbOp<false> rho2emb_op{ *parameters, (*eam_ghost) ? nullptr : grid->particle_ghost_flag_data(), eam_scratch->m_pair_enabled.data()
                                          , types.m_flat_array_ptr, rho_emb.m_flat_array_ptr };
            parallel_for( n_particles , rho2emb_op , parallel_execution_context() );
          }
        }

        if( *eam_force )
        {
          FlatSymForceOp< newtonSym.value, XFormT > force_op =
            { *parameters, rc2
            , flat_nbh_list->m_neighbor_offset.data(), flat_nbh_list->m_neighbor_list.data(), flat_nbh_list->m_half_count.data()
            , (*eam_ghost) ? nullptr : grid->particle_ghost_flag_data() , eam_scratch->m_pair_enabled.data()
            , rho_emb.m_flat_array_ptr, types.m_flat_array_ptr
            , rx.m_flat_array_ptr, ry.m_flat_array_ptr, rz.m_flat_array_ptr
            , fx.m_flat_array_ptr, fy.m_flat_array_ptr, fz.m_flat_array_ptr
            , energy.m_flat_array_ptr, cp_xform, flat_particle_locks->data() };
          parallel_for( n_particles , force_op , parallel_execution_context() );            
        }

      };
  
      auto compute_eam_xform = [&]( const auto& cp_xform )
      {
        if( *eam_symmetry ) compute_eam_force( cp_xform , std::true_type{} );
        else                compute_eam_force( cp_xform , std::false_type{} );
      };

      if( domain->xform_is_identity() ) compute_eam_xform( exanb::NullXForm{} );
      else                              compute_eam_xform( exanb::LinearXForm{domain->xform()} );
    }

  };

  namespace tmplhelper
  {
    template<class GridT> using EamPotentialFlatName  = ::exaStamp::EamPotentialFlatName<GridT>;
  }

  // === register factories ===  
  ONIKA_AUTORUN_INIT(eam_potential_multimat_flat)
  {
    OperatorNodeFactory::instance()->register_factory( EamPotentialFlatStr , make_grid_variant_operator< tmplhelper::EamPotentialFlatName > );
  }

}

#endif // only compiled if potential supports multimaterial

