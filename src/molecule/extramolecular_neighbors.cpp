#pragma xstamp_grid_variant

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/amr/amr_grid.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/amr/amr_grid_algorithm.h>
#include <exanb/core/particle_type_pair.h>

#include <onika/cuda/cuda_context.h>
#include <onika/memory/allocator.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <exanb/particle_neighbors/chunk_neighbors_config.h>
#include <exanb/particle_neighbors/chunk_neighbors_scratch.h>
#include <exanb/particle_neighbors/chunk_neighbors_host_write_accessor.h>

#include <exanb/core/domain.h>
#include <exanb/core/xform.h>

#include <exanb/particle_neighbors/chunk_neighbors_execute.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/molecule_compute_param.h>

namespace exaStamp
{
  using namespace exanb;

  template<class CellsT, class FieldIdAccT>
  struct ExtraMolecularNeighborFilter
  {
    CellsT m_cells = {};
    const ChemicalPairPotMap & chemPotMap;
    FieldIdAccT m_field_id = {};
    inline bool operator () (double d2, double rcut2,size_t cell_a,size_t p_a,size_t cell_b,size_t p_b) const
    {
      const uint64_t id_a = m_cells[cell_a][m_field_id][p_a];
      const uint64_t id_b = m_cells[cell_b][m_field_id][p_b];
      onika::oarray_t<uint64_t,2> pair = { id_a , id_b };
      if( pair[0] > pair[1] )
      {
        pair = { id_b , id_a };
      }
      assert( pair[0] != pair[1] );
      const auto it = chemPotMap.find( pair );
      return ( it == chemPotMap.end() ) && ( d2 > 0.0 ) && ( d2 < rcut2 );
    }
  };

  template<typename GridT >
  class ExtraMolecularNeighbors : public OperatorNode
  {
    ADD_SLOT( GridT               , grid            , INPUT );
    ADD_SLOT( AmrGrid             , amr             , INPUT );
    ADD_SLOT( AmrSubCellPairCache , amr_grid_pairs  , INPUT );
    ADD_SLOT( Domain              , domain          , INPUT );
    ADD_SLOT( double              , nbh_dist_lab    , INPUT );
    ADD_SLOT( GridChunkNeighbors  , chunk_neighbors , INPUT_OUTPUT );

    ADD_SLOT( ChemicalPairPotMap , chemical_pair_pot_map , INPUT , DocString{"Map of intramolecular pair potential parameter set and associated weights"} );

    ADD_SLOT( ChunkNeighborsConfig, config , INPUT_OUTPUT, ChunkNeighborsConfig{} );
    ADD_SLOT( ChunkNeighborsScratchStorage, chunk_neighbors_scratch, PRIVATE );

  public:
    inline void execute () override final
    {
      static const ChemicalPairPotMap empty_pair_pot_map = {};

      if( config->chunk_size != 1 )
      {
        lerr << "Warning: chunk_size must be 1 (actual value is "<<config->chunk_size<<"), forcing it to 1." << std::endl;
        config->chunk_size = 1;
      }
      
      const unsigned int cs = 1;
      const unsigned int cs_log2 = 0;
      const bool gpu_enabled = ( global_cuda_ctx() != nullptr ) ? global_cuda_ctx()->has_devices() : false;
      
      LinearXForm xform_filter = {domain->xform()};
      static constexpr std::false_type no_zorder = {};
      
      ldbg << "build extramolecular neighbors using id field" << std::endl;
      auto cells = grid->cells_accessor();
      auto field_id = grid->field_const_accessor( field::id );
      ExtraMolecularNeighborFilter< decltype(cells) , decltype(field_id) > nbh_filter = { cells , chemical_pair_pot_map.has_value() ? (*chemical_pair_pot_map) : empty_pair_pot_map , field_id };
      chunk_neighbors_execute(ldbg,*chunk_neighbors,*grid,*amr,*amr_grid_pairs,*config,*chunk_neighbors_scratch,cs,cs_log2,*nbh_dist_lab, xform_filter, gpu_enabled, no_zorder, nbh_filter );
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(extramolecular_neighbors)
  {
   OperatorNodeFactory::instance()->register_factory("extramolecular_neighbors", make_grid_variant_operator< ExtraMolecularNeighbors > );
  }

}

