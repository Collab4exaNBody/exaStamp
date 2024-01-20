#pragma xstamp_cuda_enable

#pragma xstamp_grid_variant

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

#include "eam_buffer.h"
#include "potential.h"
#include "eam_force_op_singlemat.h"

// this allows for parallel compilation of templated operator for each available field set



#define EamPotentialOperatorName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force)
#define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
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
    ADD_SLOT( USTAMP_POTENTIAL_PARMS, parameters        , INPUT );
    ADD_SLOT( double                , rcut              , INPUT );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                , ghost_dist_max    , INPUT_OUTPUT , 0.0 );   
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
 
    ADD_SLOT( EamPotentialScratch   , eam_scratch       , PRIVATE );

  public:

    // Operator execution
    inline void execute () override final
    {
      ///MeamPotential meam( *rcut, *parameters );
      *rcut_max = std::max( *rcut , *rcut_max );
      *ghost_dist_max = std::max( *ghost_dist_max , (*rcut) * 2.0 );
 
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
      
      double phiCut = 0.;
      double rhoCut = 0.;
#     ifdef USTAMP_POTENTIAL_RCUT
      {
        double Rho = 0.;
        double dRho = 0.;
        double Phi = 0.;
        double dPhi = 0.;
        USTAMP_POTENTIAL_EAM_RHO( *parameters, *rcut, Rho, dRho );          
        USTAMP_POTENTIAL_EAM_PHI( *parameters, *rcut, Phi, dPhi );
        phiCut = Phi;
        rhoCut = Rho;
      }
#     endif

      // TODO: emb data should be an in/out slot of this operator, not a member
      //auto cells = grid->cells();
      /*
      eam_scratch->m_offset.resize( n_cells );
      size_t total_particles = 0;
      for(size_t i=0;i<n_cells;i++)
      {
        eam_scratch->m_offset[i] = total_particles;
        total_particles += cells[i].size();
      }
      */
      eam_scratch->m_emb.clear(); // avoid useless copy of scratch values (we don't need previous values)
      eam_scratch->m_emb.resize( grid->number_of_particles() );

      // common bricks for both compute passes
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      ComputePairNullWeightIterator cp_weight {};
      ComputePairOptionalLocks<false> cp_locks;

      // 1st pass parameters : compute per particle EMB term, including ghost particles
      using EmbCPBuf = ComputePairBuffer2<false,false ,NoExtraStorage,DefaultComputePairBufferAppendFunc>;
      ComputePairBufferFactory< EmbCPBuf > emb_buf;
      EmbOp emb_op { *parameters, rhoCut, phiCut, grid->cell_particle_offset_data() /*eam_scratch->m_offset.data()*/ , eam_scratch->m_emb.data() };

      // 2nd pass parameters: compute final force using the emb term, only for non ghost particles (but reading EMB terms from neighbor particles)
      using ForceCPBuf = ComputePairBuffer2<false,false,EamComputeBufferExt,EamCopyParticleEmb>;
      ComputePairBufferFactory< ForceCPBuf , EamCopyParticleEmbInitFunc > force_buf = { grid->cell_particle_offset_data() /*eam_scratch->m_offset.data()*/ , eam_scratch->m_emb.data() };      
      ForceOp force_op { *parameters, rhoCut, phiCut, grid->cell_particle_offset_data()/*eam_scratch->m_offset.data()*/ , eam_scratch->m_emb.data() };

      field_accessor_tuple_from_field_set_t< FieldSet<field::_ep> > emb_fields = {};
      field_accessor_tuple_from_field_set_t<ComputeFields> cp_fields = {};

      // execute the 2 passes
      if( domain->xform_is_identity() )
      {
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , exanb::NullXForm{} , cp_locks );
        compute_cell_particle_pairs( *grid, *rcut, true, optional, emb_buf, emb_op, emb_fields , parallel_execution_context() );
        compute_cell_particle_pairs( *grid, *rcut, false, optional, force_buf, force_op, cp_fields , parallel_execution_context() );
      }
      else
      {
        ldbg << "xform = " << domain->xform() << std::endl;
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , exanb::LinearXForm{ domain->xform() } , cp_locks );
        compute_cell_particle_pairs( *grid, *rcut, true, optional, emb_buf, emb_op, emb_fields , parallel_execution_context() );
        compute_cell_particle_pairs( *grid, *rcut, false, optional, force_buf, force_op, cp_fields , parallel_execution_context() );
      }
    }

  };

  namespace tmplhelper
  {
    template<class GridT> using EamPotentialOperatorName = ::exaStamp::EamPotentialOperatorName<GridT>;
  }

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( EamPotentialStr , make_grid_variant_operator< tmplhelper::EamPotentialOperatorName > );
  }

}

