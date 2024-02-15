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

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_cell_particles.h>

#ifdef USTAMP_POTENTIAL_EAM_MM // operator compiled only if potential is multimaterial

#include "eam_force_op_multimat.h"

#ifndef USTAMP_POTENTIAL_EAM_MM_INIT_TYPES
#define USTAMP_POTENTIAL_EAM_MM_INIT_TYPES(p,nt,pe) /**/
#endif

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    bool ComputeRho2Emb = false ,
    bool ComputeEmb = true ,
    bool ComputeGhostEmbRho = true ,
    bool ComputeForce = true ,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type >
    >
  class EamPotentialOperatorName : public OperatorNode
  {  
    using CellParticles = typename GridT::CellParticles;

    using EamScratch = EamMultimatPotentialScratch< USTAMP_POTENTIAL_PARMS >;
    using StringVector = std::vector< std::string >;
    using EmbOp = PRIV_NAMESPACE_NAME::EmbOp;
    using ForceOp = PRIV_NAMESPACE_NAME::ForceOp;

    // compile time constant indicating if grid has virial field
    static inline constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type, field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static inline constexpr ComputeFields force_compute_fields_v{};

    // ========= I/O slots =======================
    ADD_SLOT( ParticleSpecies       , species          , INPUT , REQUIRED );
    ADD_SLOT( USTAMP_POTENTIAL_PARMS, parameters       , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( StringVector          , types            , INPUT , StringVector{} , DocString{"Empty list means all types are used, otherwise list types handled by this potential"} );
    ADD_SLOT( double                , rcut             , INPUT );
    ADD_SLOT( double                , rcut_max         , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                , ghost_dist_max   , INPUT_OUTPUT , 0.0 );   
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors  , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( GridT                 , grid             , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain           , INPUT , REQUIRED );
        
    ADD_SLOT( EamAdditionalFields   , eam_extra_fields , INPUT_OUTPUT );
    ADD_SLOT( EamScratch            , eam_scratch      , PRIVATE );

  public:

    // Operator execution
    inline void execute () override final
    {
      //MeamPotential meam( *rcut, *parameters );
      *rcut_max = std::max( *rcut , *rcut_max );
      if constexpr ( ComputeGhostEmbRho )
      {
        *ghost_dist_max = std::max( *ghost_dist_max , (*rcut) * 2.0 );
      }
 
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
      
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

      // execute the 2 passes
      auto compute_eam_force = [&]( const auto& cp_xform )
      {
        // common bricks for both compute passes
        exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
        ComputePairNullWeightIterator cp_weight {};
        ComputePairOptionalLocks<false> cp_locks;

        // 1st pass (new) computes rho then emb, without compute buffer
        if constexpr ( ComputeRho2Emb )
        {
          eam_extra_fields->m_rho.clear();
          eam_extra_fields->m_rho.resize( grid->number_of_particles() );
          double * rho_ptr = eam_extra_fields->m_emb.data();
          auto rho_field = make_external_field_flat_array_accessor( *grid , rho_ptr , field::rho );

          eam_extra_fields->m_emb.clear();
          eam_extra_fields->m_emb.resize( grid->number_of_particles() );
          double * emb_ptr = eam_extra_fields->m_emb.data();
          auto emb_field = make_external_field_flat_array_accessor( *grid , emb_ptr , field::dEmb );

          auto rho_op_fields = make_field_tuple_from_field_set( FieldSet< field::_ep , field::_type >{} , emb_field );
          auto rho_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks, ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{} );
          ComputePairBufferFactory< ComputePairBuffer2<> > rho_buf_factory = {}
          SymRhoOp<false> rho_op { *parameters , eam_scratch->m_pair_enabled.data() };
          compute_cell_particle_pairs( *grid, *rcut, ComputeGhostEmbRho, rho_optional, rho_buf_factory, rho_op, rho_op_fields, parallel_execution_context() );
          
          Rho2EmbOp rho2emb_op { *parameters , eam_scratch->m_pair_enabled.data() };
          auto rho2emb_op_fields = make_field_tuple_from_field_set( FieldSet< field::_ep , field::_type >{} , rho_field , emb_field );
          compute_cell_particles( *grid , ComputeGhostEmbRho , rho2emb_op , rho2emb_op_fields , parallel_execution_context() );
        }

        // 1st (legacy) pass parameters : compute per particle EMB term, using compute buffer
        if constexpr ( ComputeEmb )
        {
          eam_extra_fields->m_emb.clear();
          eam_extra_fields->m_emb.resize( grid->number_of_particles() );

          using EmbCPBuf = SimpleNbhComputeBuffer< FieldSet<field::_type> >; 
          ComputePairBufferFactory< EmbCPBuf > emb_buf;
          double * emb_ptr = eam_extra_fields->m_emb.data();
          auto emb_field = make_external_field_flat_array_accessor( *grid , emb_ptr , field::dEmb );
          auto emb_op_fields = make_field_tuple_from_field_set( FieldSet< field::_ep , field::_type >{} , emb_field );
          auto emb_nbh_fields = onika::make_flat_tuple( field::type ); // we want internal type field for each neighbor atom
          auto emb_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks, ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{}, emb_nbh_fields );
          EmbOp emb_op { *parameters , eam_scratch->m_pair_enabled.data() };
          compute_cell_particle_pairs( *grid, *rcut, ComputeGhostEmbRho, emb_optional, emb_buf, emb_op, emb_op_fields, parallel_execution_context() );
        }

        // 2nd pass parameters: compute final force using the emb term, only for non ghost particles (but reading EMB terms from neighbor particles)
        if constexpr ( ComputeForce )
        {
          using ForceCPBuf = SimpleNbhComputeBuffer< FieldSet< field::_type , field::_dEmb > >; //ComputePairBuffer2<false,false,EamComputeBufferEmbTypeExt,EamCopyParticleEmbType>;
          ComputePairBufferFactory<ForceCPBuf> force_buf = {};
          double * c_emb_ptr = eam_extra_fields->m_emb.data();
          auto c_emb_field = make_external_field_flat_array_accessor( *grid , c_emb_ptr , field::dEmb );
          auto force_op_fields = make_field_tuple_from_field_set( force_compute_fields_v , c_emb_field );
          auto force_nbh_fields = onika::make_flat_tuple( field::type , c_emb_field); // we want type and emb for each neighbor atom
          auto force_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks, ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{}, force_nbh_fields );
          ForceOp force_op { *parameters , eam_scratch->m_pair_enabled.data() };
          compute_cell_particle_pairs( *grid, *rcut, false, force_optional, force_buf, force_op, force_op_fields, parallel_execution_context() );
        }
      };

      if( domain->xform_is_identity() ) compute_eam_force( exanb::NullXForm{} );
      else                              compute_eam_force( exanb::LinearXForm{ domain->xform() } );
    }

  };

  namespace tmplhelper
  {
    template<class GridT> using EamPotentialOperatorName  = ::exaStamp::EamPotentialOperatorName<GridT,false,true,true,true>;
    template<class GridT> using EamPotentialRhoOnlyName   = ::exaStamp::EamPotentialOperatorName<GridT,true,false,false,false>;
    template<class GridT> using EamPotentialEmbOnlyName   = ::exaStamp::EamPotentialOperatorName<GridT,false,true,false,false>;
    template<class GridT> using EamPotentialForceOnlyName = ::exaStamp::EamPotentialOperatorName<GridT,false,false,false,true>;
  }

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( EamPotentialStr          , make_grid_variant_operator< tmplhelper::EamPotentialOperatorName  > );
    OperatorNodeFactory::instance()->register_factory( EamPotentialEmbOnlyStr   , make_grid_variant_operator< tmplhelper::EamPotentialEmbOnlyName   > );
    OperatorNodeFactory::instance()->register_factory( EamPotentialRhoOnlyStr   , make_grid_variant_operator< tmplhelper::EamPotentialEmbOnlyName   > );
    OperatorNodeFactory::instance()->register_factory( EamPotentialForceOnlyStr , make_grid_variant_operator< tmplhelper::EamPotentialForceOnlyName > );
  }

}

#endif // only compiled if potential supports multimaterial

