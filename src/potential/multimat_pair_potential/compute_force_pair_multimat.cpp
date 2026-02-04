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



#include <exanb/core/make_grid_variant_operator.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/cpp_utils.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential_factory/pair_potential_yaml.h>
#include <exanb/core/particle_type_pair.h>
#include <exanb/core/domain.h>
#include <exanb/core/compact_grid_pair_weights.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_pair_multimat.h>

// this allows for parallel compilation of templated operator for each available field set


namespace exaStamp
{
  using namespace exanb;

  // helper function to get the right force operator  
  template<bool has_virial>
  static inline 
  std::shared_ptr< std::conditional_t<has_virial,PairPotentialComputeVirialOperator,PairPotentialComputeOperator> >
  get_potential_force_op( PairPotential& pot, std::integral_constant<bool,has_virial> )
  {
    if constexpr ( has_virial ) return pot.force_virial_op();
    if constexpr ( ! has_virial ) return pot.force_op();
    return nullptr;
  }

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type >
    >
  class ComputeForcePairMultimatChunk : public OperatorNode
  {
    using NeighborPairWeight = std::vector< std::vector< double > >;

    // compile time constant indicating if grid has virial field
    using has_virial_field_t = GridHasField<GridT,field::_virial>;
    static constexpr has_virial_field_t has_virial_field{};

    // attributes processed during computation
    using WorkFieldFieldSet = std::conditional_t< has_virial_field , PairPotentialVirialFieldSet , PairPotentialFieldSet >;  
    static constexpr WorkFieldFieldSet s_compute_field_set {};

    using BaseForceOp = ComputePairOperator<WorkFieldFieldSet>;
    using BaseForceOps = std::vector< std::shared_ptr<BaseForceOp> >;
    using PairForceOps = std::conditional_t< has_virial_field , std::vector< std::shared_ptr<PairPotentialComputeVirialOperator> > , std::vector< std::shared_ptr<PairPotentialComputeOperator> > >;    

    // ========= I/O slots =======================
    ADD_SLOT( Domain             , domain          , INPUT , REQUIRED );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool               , symetric        , INPUT , false );
    ADD_SLOT( bool               , ghost           , INPUT , false );
//    ADD_SLOT( NeighborPairWeight , nbh_weight      , INPUT , OPTIONAL );
    ADD_SLOT( CompactGridPairWeights , compact_nbh_weight      , INPUT , OPTIONAL );

    ADD_SLOT( ParticleSpecies    , species         , INPUT , REQUIRED );
    ADD_SLOT( UserMultiPairPotentials , potentials , INPUT_OUTPUT , REQUIRED );
    
    ADD_SLOT( GridT              , grid            , INPUT_OUTPUT );
    ADD_SLOT( double             , rcut_max        , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( bool               , enable_pair_weights, INPUT, true ); 

  public:

    inline void execute () override final
    {
      assert( chunk_neighbors->m_cell_stream.size() == grid->number_of_cells() );

      /* build up the nul potential we might need to fill in undefined potential pairs */
      YAML::Node zero_pot_node(YAML::NodeType::Map);
      zero_pot_node["potential"] = "zero";
      auto zero_pot = PairPotentialFactory::make_instance( zero_pot_node );
      auto zero_pot_vir = get_potential_force_op( *zero_pot , has_virial_field );
      zero_pot_vir->set_rcut( 0.0 );

      /* assemble and optimize potential functors from user parameters */
      //if(  )
      auto compileInfo = potentials->compile_potentials_for_species(*species,zero_pot_vir,has_virial_field);
      const unsigned int NTypes = compileInfo.NTypes;
      *rcut_max = std::max( *rcut_max , compileInfo.rcut_max );
      auto& compiled_potentials = potentials->compiled_potentials( has_virial_field );

      bool has_weights = *enable_pair_weights && compact_nbh_weight.has_value() ;
      if( has_weights ) { has_weights = has_weights && !compact_nbh_weight->empty(); }

      // only non symetric case handled by now
      if( *symetric )
      {
        lerr << "Symmetric multi-material potential not available, aborting."<<std::endl;
        std::abort();
      }
      
      // optional locking
      ComputePairOptionalLocks<false> cp_locks {};

      // neighbor list traversal
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };

      // this optional indirection table allows for reduction of number of potential functions (less than possible number of pairs)
      const auto& pair_id_map = compiled_potentials.m_pair_id_map;

      if( has_weights )
      {
        ldbg << "USING WEIGHTS" << std::endl;

        assert( compact_nbh_weight.has_value() );
        //ComputePairWeightIterator cp_weight { *nbh_weight };
        CompactPairWeightIterator cp_weight { compact_nbh_weight->m_cell_weights.data() };
        
        auto force_buf = make_compute_pair_buffer< ComputePairBuffer2<true,false> >();
        if( domain->xform_is_identity() )
        {
          NullXForm cp_xform;
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
          compute_pair_multimat( *grid,compiled_potentials.m_rcuts,compiled_potentials.m_force_ops,NTypes,pair_id_map, *ghost, optional, force_buf, s_compute_field_set );
        }
        else
        {
          ldbg << "XFORM = " << domain->xform() << std::endl;
          LinearXForm cp_xform { domain->xform() };
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
          compute_pair_multimat( *grid,compiled_potentials.m_rcuts,compiled_potentials.m_force_ops,NTypes,pair_id_map, *ghost, optional, force_buf, s_compute_field_set );
        }
      }
      else
      {
        ComputePairNullWeightIterator cp_weight {};
        auto force_buf = make_compute_pair_buffer< ComputePairBuffer2<false,false> >();
        if( domain->xform_is_identity() )
        {
          NullXForm cp_xform;
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
          compute_pair_multimat( *grid,compiled_potentials.m_rcuts,compiled_potentials.m_force_ops,NTypes,pair_id_map, *ghost, optional, force_buf, s_compute_field_set );
        }
        else
        {
          ldbg << "XFORM = " << domain->xform() << std::endl;
          LinearXForm cp_xform { domain->xform() };
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
          compute_pair_multimat( *grid,compiled_potentials.m_rcuts,compiled_potentials.m_force_ops,NTypes,pair_id_map, *ghost, optional, force_buf, s_compute_field_set );
        }
      }

    }

  };

  template<class GridT> using ComputeForcePairMultimatChunkTmpl = ComputeForcePairMultimatChunk<GridT>;

   // === register factories ===  
  ONIKA_AUTORUN_INIT(compute_force_pair_multimat)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_force_pair_multimat", make_grid_variant_operator< ComputeForcePairMultimatChunkTmpl > );
  }

}


