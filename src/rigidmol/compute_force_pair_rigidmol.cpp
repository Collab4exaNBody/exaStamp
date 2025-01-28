

#include <exanb/core/make_grid_variant_operator.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/cpp_utils.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential_factory/pair_potential_yaml.h>
#include <exanb/core/particle_type_pair.h>
#include <exanb/core/domain.h>
#include <exanb/core/compact_grid_pair_weights.h>

#include "compute_pair_rigidmol.h"
#include "rigid_molecule_compute.h"

#include <exanb/particle_neighbors/chunk_neighbors.h>


// this allows for parallel compilation of templated operator for each available field set


namespace exaStamp
{
  using namespace exanb;

  // helper function to get the right force operator  
  static inline std::shared_ptr<PairPotentialComputeOperator> get_potential_force_op( PairPotential& pot, std::false_type ) { return pot.force_op(); }
  static inline std::shared_ptr<PairPotentialComputeVirialOperator> get_potential_force_op( PairPotential& pot, std::true_type ) { return pot.force_virial_op(); }

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type, field::_orient, field::_couple >
    >
  class ComputeForcePairRigidMol : public OperatorNode
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
    using PairForceOps = std::conditional_t<has_virial_field, std::vector< std::shared_ptr<PairPotentialComputeVirialOperator> > , std::vector< std::shared_ptr<PairPotentialComputeOperator> > >;    

    // ========= I/O slots =======================
    ADD_SLOT( Domain             , domain          , INPUT , REQUIRED );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool               , ghost           , INPUT , false );
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

      size_t ns = species->size();
      double rigid_mol_max_radius = 0.0;
      for(size_t i=0;i<ns;i++)
      {
        auto & rigidmol = species->at(i);
        size_t na = rigidmol.m_rigid_atom_count;
        for(size_t j=0;j<na;j++)
        {
          rigid_mol_max_radius = std::max( rigid_mol_max_radius , norm( rigidmol.m_rigid_atoms[j].m_pos ) );
        }
      }
      ldbg << "rigidmol max radius = " << rigid_mol_max_radius << std::endl;

      if( grid->number_of_cells() == 0 )
      {
        ldbg << "empty grid, not compiling potentials, only parsing rcuts" << std::endl;
        for(const auto& pot : potentials->m_user_potentials)
        {
          ldbg << "add rcut "<<pot.m_rcut<<std::endl;
          *rcut_max = std::max( *rcut_max , pot.m_rcut + 2 * rigid_mol_max_radius );
        }
        return;
      }

      /* assemble and optimize potential functors from user parameters */
      auto compileInfo = potentials->compile_potentials_for_species(*species,zero_pot_vir,has_virial_field);
      const unsigned int NTypes = compileInfo.NTypes;
      *rcut_max = std::max( *rcut_max , compileInfo.rcut_max + 2 * rigid_mol_max_radius );
      auto& compiled_potentials = potentials->compiled_potentials(has_virial_field);

      bool has_weights = *enable_pair_weights && compact_nbh_weight.has_value() ;
      if( has_weights ) { has_weights = has_weights && !compact_nbh_weight->empty(); }

      // no weighting support by now
      if( has_weights )
      {
        lerr << "Pair weighting not supported for rigid molecules yet" << std::endl;
        std::abort();
      }
      
      // optional locking
      ComputePairOptionalLocks<false> cp_locks {};

      // neighbor list traversal
      //auto cells = grid->cells();
      //using CellT = std::remove_cv_t< std::remove_reference_t< decltype(cells[0]) > >;
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it = { *chunk_neighbors };
      //ChunkParticleNeighborsIterator<CellT,false> nbh_it = { cells , chunk_neighbors->m_data , grid->dimension() , chunk_neighbors->m_chunk_size };

      // this optional indirection table allows for reduction of number of potential functions (less than possible number of pairs)
      const auto& pair_id_map = compiled_potentials.m_pair_id_map;

      ComputePairNullWeightIterator cp_weight {};
      auto force_buf = make_compute_pair_buffer< ComputePairBuffer2<false,false> >();
      LinearXForm cp_xform = { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it , cp_weight , cp_xform, cp_locks );
      RigidMoleculeCompute molfunc {};
      compute_pair_rigidmol( *grid,compiled_potentials.m_rcuts,compiled_potentials.m_force_ops,NTypes,pair_id_map,*species, *ghost, molfunc, optional, force_buf, s_compute_field_set );
    }

  };

  template<class GridT> using ComputeForcePairRigidMolTmpl = ComputeForcePairRigidMol<GridT>;

   // === register factories ===  
  ONIKA_AUTORUN_INIT(compute_force_pair_rigidmol)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_force_pair_rigidmol", make_grid_variant_operator< ComputeForcePairRigidMolTmpl > );
  }

}


