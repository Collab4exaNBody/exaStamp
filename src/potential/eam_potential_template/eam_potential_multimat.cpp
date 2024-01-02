#pragma xstamp_grid_variant

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

#include "eam_buffer.h"
#include "eam_yaml.h"
#include "potential.h"
#include "eam_force_op_multimat.h"

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

// this allows for parallel compilation of templated operator for each available field set


namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type >
    >
  class EamPotentialOperatorName : public OperatorNode
  {  
    using CellParticles = typename GridT::CellParticles;
    using EamMultiMatParams = EamMultimatParameters<USTAMP_POTENTIAL_PARMS>;
    using EamParamVector = std::vector<EamMultiMatParams>;
    using EamPotParmRO = PRIV_NAMESPACE_NAME::EamMultiMatParamsReadOnly;
    using EamScratch = EamMultimatPotentialScratch< EamPotParmRO >;
    using EmbOp = PRIV_NAMESPACE_NAME::EmbOp;
    using ForceOp = PRIV_NAMESPACE_NAME::ForceOp;

    // compile time constant indicating if grid has virial field
    static inline constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type, field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static inline constexpr ComputeFields force_compute_fields_v{};
    static inline constexpr FieldSet<field::_ep,field::_type> emb_compute_fields_v{};

    // ========= I/O slots =======================
    ADD_SLOT( ParticleSpecies       , species          , INPUT , REQUIRED );
    ADD_SLOT( EamParamVector        , potentials       , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double                , rcut             , INPUT );
    ADD_SLOT( double                , rcut_max         , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( double                , ghost_dist_max   , INPUT_OUTPUT , 0.0 );   
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors  , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( GridT                 , grid             , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain           , INPUT , REQUIRED );
    ADD_SLOT( EamScratch            , eam_scratch      , PRIVATE );

  public:

    // Operator execution
    inline void execute () override final
    {
      //MeamPotential meam( *rcut, *parameters );
      *rcut_max = std::max( *rcut , *rcut_max );
      *ghost_dist_max = std::max( *ghost_dist_max , (*rcut) * 2.0 );
 
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
      
      size_t n_species = species->size();
      size_t n_type_pairs = unique_pair_count( n_species );
      bool compile_eam_parameters = (potentials->size() != n_type_pairs);
      if( ! compile_eam_parameters )
      {
        for(size_t i=0;i<n_species;i++) if( !potentials->at(i).m_type_a.empty() || !potentials->at(i).m_type_b.empty() ) compile_eam_parameters = true;
      }
      if( compile_eam_parameters )
      {
        std::map<size_t , USTAMP_POTENTIAL_PARMS > ordered_parameters;
        std::unordered_map< std::string , size_t > name_to_type;
        for(size_t i=0;i<n_species;i++)
        {
          name_to_type[ species->at(i).m_name ] = i;
        }
        for(const auto& p : *potentials)
        {
          if( name_to_type.find(p.m_type_a) == name_to_type.end() )
          {
            lerr << "Atom type name '"<<p.m_type_a<<"' is invalid in "<< EamPotentialStr << std::endl;
            std::abort();
          }
          if( name_to_type.find(p.m_type_b) == name_to_type.end() )
          {
            lerr << "Atom type name '"<<p.m_type_b<<"' is invalid in "<< EamPotentialStr << std::endl;
            std::abort();
          }
          size_t ta = name_to_type[p.m_type_a];
          size_t tb = name_to_type[p.m_type_b];
          size_t pair_id = unique_pair_id( ta , tb );
          ordered_parameters[pair_id] = p.m_parameters;
        }
	
        potentials->clear();
        potentials->resize( n_type_pairs );
	
        eam_scratch->m_pair_enabled.assign( n_type_pairs , false );
        for(const auto& p : ordered_parameters)
        {
          unsigned int ta=0, tb=0;
          pair_id_to_type_pair(p.first, ta, tb);
          assert( p.first >= 0 && p.first < n_type_pairs );
          eam_scratch->m_pair_enabled[ p.first ] = true;
          auto & pot = potentials->at( p.first );
          pot.m_parameters = p.second;
          pot.m_specy_pair.m_charge_a = species->at(ta).m_charge;
          pot.m_specy_pair.m_charge_b = species->at(tb).m_charge;
        }	
      }

      assert( potentials->size() == n_type_pairs );
      eam_scratch->m_phi_rho_cutoff.resize( n_type_pairs );
      eam_scratch->m_ro_potentials.resize( n_type_pairs );
      for(size_t i=0;i<n_type_pairs;i++)
      {
        double phiCut = 0.;
        double rhoCut = 0.;
#       ifdef USTAMP_POTENTIAL_RCUT
        if( eam_scratch->m_pair_enabled[i] )
        {
          double Rho = 0.;
          double dRho = 0.;
          double Phi = 0.;
          double dPhi = 0.;
          const auto & pot = potentials->at(i);
          USTAMP_POTENTIAL_EAM_RHO_MM( pot.m_parameters, *rcut, Rho, dRho , pot.m_specy_pair );
          USTAMP_POTENTIAL_EAM_PHI_MM( pot.m_parameters, *rcut, Phi, dPhi , pot.m_specy_pair );
          phiCut = Phi;
          rhoCut = Rho;
        }
#       endif
        eam_scratch->m_phi_rho_cutoff[i] = { phiCut , rhoCut };
        const auto & pot = potentials->at(i);
        eam_scratch->m_ro_potentials[i] = { pot.m_parameters , pot.m_specy_pair };
      }

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
      using EmbCPBuf = ComputePairBuffer2<false,false ,EamComputeBufferTypeExt,EamCopyParticleType>;
      ComputePairBufferFactory< EmbCPBuf > emb_buf;
      EmbOp emb_op { eam_scratch->m_ro_potentials.data(), eam_scratch->m_phi_rho_cutoff.data(), grid->cell_particle_offset_data() /*eam_scratch->m_offset.data()*/, eam_scratch->m_emb.data(), eam_scratch->m_pair_enabled.data() };

      // 2nd pass parameters: compute final force using the emb term, only for non ghost particles (but reading EMB terms from neighbor particles)
      using ForceCPBuf = ComputePairBuffer2<false,false,EamComputeBufferEmbTypeExt,EamCopyParticleEmbType>;
      ComputePairBufferFactory< ForceCPBuf , EamCopyParticleEmbInitFunc > force_buf = { grid->cell_particle_offset_data() /*eam_scratch->m_offset.data()*/ , eam_scratch->m_emb.data() };
      ForceOp force_op { eam_scratch->m_ro_potentials.data(), eam_scratch->m_phi_rho_cutoff.data(), grid->cell_particle_offset_data()/*eam_scratch->m_offset.data()*/ , eam_scratch->m_emb.data(), eam_scratch->m_pair_enabled.data() };

      // execute the 2 passes
      if( domain->xform_is_identity() )
      {
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , NullXForm{} , cp_locks );
        compute_cell_particle_pairs( *grid, *rcut, true, optional, emb_buf, emb_op, emb_compute_fields_v, parallel_execution_context() );
        compute_cell_particle_pairs( *grid, *rcut, false, optional, force_buf, force_op, force_compute_fields_v, parallel_execution_context() );
      }
      else
      {
        ldbg << "xform = " << domain->xform() << std::endl;
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , LinearXForm{ domain->xform() } , cp_locks );
        compute_cell_particle_pairs( *grid, *rcut, true, optional, emb_buf, emb_op, emb_compute_fields_v, parallel_execution_context() );
        compute_cell_particle_pairs( *grid, *rcut, false, optional, force_buf, force_op, force_compute_fields_v, parallel_execution_context() );
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

