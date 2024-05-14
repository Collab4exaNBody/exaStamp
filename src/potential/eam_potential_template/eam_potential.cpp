// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

// #pragma xstamp_grid_variant // DO NOT REMOVE THIS LINE

#include <exaStamp/potential/eam/eam_buffer.h>

#include <exanb/core/grid.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/domain.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/xform.h>

#include "potential.h"

#ifndef USTAMP_POTENTIAL_EAM_MM // only for monomaterial EAM potentials

#include "eam_force_op_singlemat.h"

#define EamPotentialOperatorName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force)
#define EamPotentialComputeEmbNoGhostName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_emb)
#define EamPotentialComputeForceOnlyName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force_reuse_emb)

#define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)
#define EamPotentialEmbNoGhostStr USTAMP_STR(EamPotentialComputeEmbNoGhostName)
#define EamPotentialForceOnlyStr USTAMP_STR(EamPotentialComputeForceOnlyName)

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    bool ComputeEmb = true ,
    bool ComputeGhostEmb = true ,
    bool ComputeForce = true ,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class EamPotentialOperatorName : public OperatorNode
  {  
    using CellParticles = typename GridT::CellParticles;
    using EmbOp = PRIV_NAMESPACE_NAME::EmbOp;
    using ForceOp = PRIV_NAMESPACE_NAME::ForceOp;
//    using ParticleLock = decltype( ComputePairOptionalLocks<false>{}[0][0] );

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;

    // ========= I/O slots =======================
    ADD_SLOT( USTAMP_POTENTIAL_PARMS, parameters        , INPUT_OUTPUT );
    ADD_SLOT( double                , rcut              , INPUT );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                , ghost_dist_max    , INPUT_OUTPUT , 0.0 );   
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
 
    ADD_SLOT( EamAdditionalFields   , eam_extra_fields , INPUT_OUTPUT );

  public:

    // Operator execution
    inline void execute () override final
    {
      ///MeamPotential meam( *rcut, *parameters );
      *rcut_max = std::max( *rcut , *rcut_max );
      if constexpr ( ComputeGhostEmb )
      {
        *ghost_dist_max = std::max( *ghost_dist_max , (*rcut) * 2.0 );
      }
 
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
            
      auto compute_eam_force = [&]( const auto& cp_xform )
      {
        // common building blocks for both compute passes
        exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
        ComputePairNullWeightIterator cp_weight {};
        ComputePairOptionalLocks<false> cp_locks;

        // 1st pass parameters : compute per particle EMB term, including ghost particles
        if constexpr ( ComputeEmb )
        {
          // reset emb field to zero
          eam_extra_fields->m_rho_emb.clear();
          eam_extra_fields->m_rho_emb.resize( grid->number_of_particles() );
          //eam_extra_fields->m_emb.assign( grid->number_of_particles() , 0.0 );

          // build compute buffer
          using EmbCPBuf = ComputePairBuffer2<>;
          ComputePairBufferFactory< EmbCPBuf > emb_buf;
          EmbOp emb_op { *parameters };

          // Emb term computation will access, for each central atom, potential energy (internal field) and emb_field (externally stored extra field)
          auto emb_field = grid->field_accessor( field::rho_dEmb );
          auto emb_op_fields = make_field_tuple_from_field_set( FieldSet<field::_ep>{} , emb_field );
          auto emb_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks /* no additional fields required for neighbors */ );
          compute_cell_particle_pairs( *grid, *rcut, ComputeGhostEmb, emb_optional, emb_buf, emb_op, emb_op_fields , parallel_execution_context() );      
        }

        // 2nd pass parameters: compute final force using the emb term, only for non ghost particles (but reading EMB terms from neighbor particles)
        if constexpr ( ComputeForce )
        {
          using ForceCPBuf = SimpleNbhComputeBuffer< FieldSet<field::_rho_dEmb> >; /* we want extra neighbor storage space to store these fields */
          ComputePairBufferFactory< ForceCPBuf > force_buf;  
          ForceOp force_op { *parameters };
          const double * c_emb_ptr = eam_extra_fields->m_rho_emb.data();
          auto c_emb_field = grid->field_const_accessor( field::rho_dEmb );

          // force computation will access, for each central atom, fields defined in ComputeFields plus external constant field c_emb_field
          auto force_op_fields = make_field_tuple_from_field_set( ComputeFields{} , c_emb_field );
          auto force_nbh_fields = onika::make_flat_tuple(c_emb_field); // we want external field c_emb_field to be populated for each neighbor
          auto force_optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform , cp_locks, ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{}, force_nbh_fields );
          compute_cell_particle_pairs( *grid, *rcut, false, force_optional, force_buf, force_op, force_op_fields , parallel_execution_context() );
        }
      };

      if( domain->xform_is_identity() ) compute_eam_force( exanb::NullXForm{} );
      else                              compute_eam_force( exanb::LinearXForm{ domain->xform() } );
    }

  };

  namespace tmplhelper
  {
    template<class GridT> using EamPotentialOperatorName          = ::exaStamp::EamPotentialOperatorName<GridT,true,true,true>;
    template<class GridT> using EamPotentialComputeEmbNoGhostName = ::exaStamp::EamPotentialOperatorName<GridT,true,false,false>;
    template<class GridT> using EamPotentialComputeForceOnlyName  = ::exaStamp::EamPotentialOperatorName<GridT,false,false,true>;
  }

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( EamPotentialStr           , make_grid_variant_operator< tmplhelper::EamPotentialOperatorName          > );
    OperatorNodeFactory::instance()->register_factory( EamPotentialEmbNoGhostStr , make_grid_variant_operator< tmplhelper::EamPotentialComputeEmbNoGhostName > );
    OperatorNodeFactory::instance()->register_factory( EamPotentialForceOnlyStr  , make_grid_variant_operator< tmplhelper::EamPotentialComputeForceOnlyName  > );
  }

}

#endif // only mono material EAM potentials will compile this code

