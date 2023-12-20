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

  template<class CellsT>
  struct ExtraMolecularNeighborFilter
  {
    CellsT m_cells;
    inline bool operator () (double d2, double rcut2,size_t cell_a,size_t p_a,size_t cell_b,size_t p_b) const
    {
      uint64_t idmol_a = molecule_instance_from_id( m_cells[cell_a][field::idmol][p_a] );
      uint64_t idmol_b = molecule_instance_from_id( m_cells[cell_b][field::idmol][p_b] );
      return ( idmol_a != idmol_b ) && ( d2 < rcut2 );
    }
  };

  template<typename GridT>
  class ExtraMolecularNeighbors : public OperatorNode
  {
    using has_idmol_field_t = typename GridT::CellParticles::template HasField < field::_idmol > ;
    static constexpr bool has_idmol_field = has_idmol_field_t::value;

#ifdef XSTAMP_CUDA_VERSION
    ADD_SLOT( onika::cuda::CudaContext , cuda_ctx , INPUT , OPTIONAL );
#endif

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
      unsigned int cs = config->chunk_size;
      unsigned int cs_log2 = 0;
      while( cs > 1 )
      {
        assert( (cs&1)==0 );
        cs = cs >> 1;
        ++ cs_log2;
      }
      cs = 1<<cs_log2;
      //ldbg << "cs="<<cs<<", log2(cs)="<<cs_log2<<std::endl;
      if( cs != static_cast<size_t>(config->chunk_size) )
      {
        lerr<<"chunk_size is not a power of two"<<std::endl;
        std::abort();
      }

      bool gpu_enabled = false;
#     ifdef XSTAMP_CUDA_VERSION
      if( cuda_ctx.has_value() )
      {
        gpu_enabled = cuda_ctx->has_devices() ;
      }
#     endif
      
      LinearXForm xform_filter = {domain->xform()};
      static constexpr std::false_type no_zorder = {};
      
      if constexpr ( has_idmol_field )
      {
        auto cells = grid->cells();
        ExtraMolecularNeighborFilter< decltype(cells) > nbh_filter = { cells };
        chunk_neighbors_execute(ldbg,*chunk_neighbors,*grid,*amr,*amr_grid_pairs,*config,*chunk_neighbors_scratch,cs,cs_log2,*nbh_dist_lab, xform_filter, gpu_enabled, no_zorder, nbh_filter );
      }
      else
      {
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

