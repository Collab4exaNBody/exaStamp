#pragma xstamp_grid_variant

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
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

namespace exaStamp
{
  using namespace exanb;

  template<class CellsT, class FieldIdAccT>
  struct ExtraMolecularNeighborFilter
  {
    CellsT m_cells = {};
    FieldIdAccT m_field_id = {};
    inline bool operator () (double d2, double rcut2,size_t cell_a,size_t p_a,size_t cell_b,size_t p_b) const
    {
      uint64_t idmol_a = molecule_instance_from_id( m_cells[cell_a][m_field_id][p_a] );
      uint64_t idmol_b = molecule_instance_from_id( m_cells[cell_b][m_field_id][p_b] );
      return ( idmol_a != idmol_b ) && ( d2 > 0.0 ) && ( d2 < rcut2 );
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

    ADD_SLOT( ChunkNeighborsConfig, config , INPUT, ChunkNeighborsConfig{} );
    ADD_SLOT( ChunkNeighborsScratchStorage, chunk_neighbors_scratch, PRIVATE );

  public:
    inline void execute () override final
    {
      if( config->chunk_size != 1 )
      {
        fatal_error() << "chunk_size must be 1, otherrwise extramolecular_neighbors is pointless" << std::endl;
      }
      
      const unsigned int cs = 1;
      const unsigned int cs_log2 = 0;
      const bool gpu_enabled = ( global_cuda_ctx() != nullptr ) ? global_cuda_ctx()->has_devices() : false;
      
      LinearXForm xform_filter = {domain->xform()};
      static constexpr std::false_type no_zorder = {};
      
      if( grid->has_allocated_field(field::idmol) )
      {
        ldbg << "build extramolecular neighbors using idmol field" << std::endl;
        auto cells = grid->cells_accessor();
        auto field_idmol = grid->field_const_accessor( field::idmol );
        ExtraMolecularNeighborFilter< decltype(cells) , decltype(field_idmol) > nbh_filter = { cells , field_idmol };
        chunk_neighbors_execute(ldbg,*chunk_neighbors,*grid,*amr,*amr_grid_pairs,*config,*chunk_neighbors_scratch,cs,cs_log2,*nbh_dist_lab, xform_filter, gpu_enabled, no_zorder, nbh_filter );
      }
      else if( grid->has_allocated_field(field::id) )
      {
        ldbg << "build extramolecular neighbors using id field" << std::endl;
        auto cells = grid->cells_accessor();
        auto field_id = grid->field_const_accessor( field::id );
        ExtraMolecularNeighborFilter< decltype(cells) , decltype(field_id) > nbh_filter = { cells , field_id };
        chunk_neighbors_execute(ldbg,*chunk_neighbors,*grid,*amr,*amr_grid_pairs,*config,*chunk_neighbors_scratch,cs,cs_log2,*nbh_dist_lab, xform_filter, gpu_enabled, no_zorder, nbh_filter );
      }
      else
      {
        ldbg << "build extramolecular neighbors without id field => complete neighbors are built" << std::endl;
        chunk_neighbors_execute(ldbg,*chunk_neighbors,*grid,*amr,*amr_grid_pairs,*config,*chunk_neighbors_scratch,cs,cs_log2,*nbh_dist_lab, xform_filter, gpu_enabled, no_zorder );
      }
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("extramolecular_neighbors", make_grid_variant_operator< ExtraMolecularNeighbors > );
  }

}

