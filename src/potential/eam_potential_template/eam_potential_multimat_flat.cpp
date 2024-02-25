// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

// #pragma xstamp_grid_variant // DO NOT REMOVE THIS LINE

#include <exanb/core/grid.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/particle_type_pair.h>
#include <exanb/core/domain.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>

#include <exaStamp/potential/eam/eam_buffer.h>
#include <exaStamp/potential/eam/eam_yaml.h>
#include "potential.h"

#include <exanb/particle_neighbors/flat_neighbor_lists.h>
#include <onika/parallel/parallel_for.h>
#include <onika/parallel/memset.h>

#ifdef USTAMP_POTENTIAL_EAM_MM // operator compiled only if potential is multimaterial

#include "eam_force_op_multimat.h"

#ifndef USTAMP_POTENTIAL_EAM_MM_INIT_TYPES
#define USTAMP_POTENTIAL_EAM_MM_INIT_TYPES(p,nt,pe) /**/
#endif

namespace exaStamp
{
  using namespace exanb;

  template< class GridT >
  class EamPotentialFlatName : public OperatorNode
  {  
    using CellParticles = typename GridT::CellParticles;

    using EamScratch = EamMultimatPotentialScratch< USTAMP_POTENTIAL_PARMS >;
    using StringVector = std::vector< std::string >;

    template<bool NewtonSym, bool Ghost, class ParticleLocksT> using FlatSymRhoOp = PRIV_NAMESPACE_NAME::FlatSymRhoOp<NewtonSym,Ghost,ParticleLocksT>;

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

    ADD_SLOT( GridParticleLocks     , particle_locks       , INPUT_OUTPUT , OPTIONAL , DocString{"particle spin locks"} );

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

      size_t n_particles = grid->number_of_particles();
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

      auto rho_emb = grid->field_accessor( field::rho_dEmb );
      //auto energy = grid->field_accessor( field::flat_ep );
      //auto fx = grid->field_accessor( field::flat_fx );
      //auto fy = grid->field_accessor( field::flat_fy );
      //auto fz = grid->field_accessor( field::flat_fz );
      auto rx = grid->field_accessor( field::flat_rx );
      auto ry = grid->field_accessor( field::flat_ry );
      auto rz = grid->field_accessor( field::flat_rz );
      auto types = grid->field_accessor( field::flat_type );

      if( *eam_rho )
      {
        const size_t n_particles = grid->number_of_particles();
        eam_scratch->m_particle_locks.resize( n_particles );
        const double rc = *rcut;
        const double rc2 = rc*rc;
        FlatSymRhoOp< true , false , spin_mutex_array > func =
          { *parameters, eam_scratch->m_pair_enabled.data(), rc2
          , flat_nbh_list->m_neighbor_offset.data(), flat_nbh_list->m_neighbor_list.data(), flat_nbh_list->m_half_count.data(), grid->particle_ghost_flag_data()
          , rho_emb.m_func.m_data_array, types.m_func.m_data_array, rx.m_func.m_data_array, ry.m_func.m_data_array, rz.m_func.m_data_array, eam_scratch->m_particle_locks };
        parallel_for( n_particles , func , parallel_execution_context() );
      }
  
    }

  };

  namespace tmplhelper
  {
    template<class GridT> using EamPotentialFlatName  = ::exaStamp::EamPotentialFlatName<GridT>;
  }

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( EamPotentialFlatStr , make_grid_variant_operator< tmplhelper::EamPotentialFlatName > );
  }

}

#endif // only compiled if potential supports multimaterial

