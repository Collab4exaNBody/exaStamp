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
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/particle_species/distance_map.h>
#include <exanb/core/particle_type_id.h>

#include <onika/cuda/cuda_context.h>
#include <onika/memory/allocator.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <exanb/particle_neighbors/chunk_neighbors_config.h>
#include <exanb/particle_neighbors/chunk_neighbors_scratch.h>
#include <exanb/particle_neighbors/chunk_neighbors_host_write_accessor.h>

#include <exanb/core/domain.h>
#include <exanb/core/xform.h>

#include <exanb/particle_neighbors/chunk_neighbors_execute.h>

namespace exaStamp
{
  using namespace exanb;

  template<class CellsT>
  struct PairTypeRCutNeighborFilter
  {
    static constexpr size_t MAX_TYPE_PAIRS = 248;
    CellsT m_cells;
    double m_pair_rcut2[MAX_TYPE_PAIRS];
    inline bool operator () (double d2, double rcut2,size_t cell_a,size_t p_a,size_t cell_b,size_t p_b) const
    {
      const unsigned int type_a = m_cells[cell_a][field::type][p_a];
      const unsigned int type_b = m_cells[cell_b][field::type][p_b];
      const unsigned int pair_id = unique_pair_id(type_a,type_b);
      assert( pair_id < MAX_TYPE_PAIRS );
      return ( d2 < m_pair_rcut2[pair_id] );
    }
  };

  template< typename GridT
          , class = AssertGridHasFields<GridT,field::_type>
          >
  class PairTypeRCutNeighbors : public OperatorNode
  {
#ifdef XNB_CUDA_VERSION
    ADD_SLOT( onika::cuda::CudaContext , cuda_ctx , INPUT , OPTIONAL );
#endif

    ADD_SLOT( GridT               , grid            , INPUT , REQUIRED );
    ADD_SLOT( AmrGrid             , amr             , INPUT , REQUIRED );
    ADD_SLOT( AmrSubCellPairCache , amr_grid_pairs  , INPUT , REQUIRED );
    ADD_SLOT( Domain              , domain          , INPUT , REQUIRED );
    ADD_SLOT( double              , nbh_dist_lab    , INPUT , REQUIRED );
    ADD_SLOT( double              , rcut_inc        , INPUT , REQUIRED );
    ADD_SLOT( GridChunkNeighbors  , chunk_neighbors , INPUT_OUTPUT );

    ADD_SLOT( ParticleSpecies     , species         , INPUT , REQUIRED );
    ADD_SLOT( ParticleTypeMap     , particle_type_map , INPUT , REQUIRED);
    ADD_SLOT( DistanceMap         , pair_distances  , INPUT , REQUIRED , DocString{"Set of pairwise distances. determine neighbor rcut for each atom type pair"} );

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
#     ifdef XNB_CUDA_VERSION
      if( cuda_ctx.has_value() )
      {
        gpu_enabled = cuda_ctx->has_devices() ;
      }
#     endif
      
      LinearXForm xform_filter = {domain->xform()};
      static constexpr std::false_type no_zorder = {};
      
      auto cells = grid->cells();
      PairTypeRCutNeighborFilter< decltype(cells) > nbh_filter = { cells };
      for(size_t i=0;i<nbh_filter.MAX_TYPE_PAIRS;i++) nbh_filter.m_pair_rcut2[i] = 0.0;
      for(const auto& pdist : *pair_distances)
      {
        std::string ts1 = pdist.first.substr(0,pdist.first.find(','));
        std::string ts2 = pdist.first.substr(pdist.first.find(',')+1);
        const unsigned int t1 = particle_type_map->at( ts1 );
        const unsigned int t2 = particle_type_map->at( ts2 );
        const unsigned int pair_id = unique_pair_id(t1,t2);
        assert( pair_id < nbh_filter.MAX_TYPE_PAIRS );
        const double nbh_rcut = pdist.second + (*rcut_inc);
        ldbg << ts1 << " / " << ts2 << " -> "<<pdist.second<<"+"<<(*rcut_inc)<<" = "<<nbh_rcut <<std::endl;
        nbh_filter.m_pair_rcut2[pair_id] = nbh_rcut * nbh_rcut;
      }

      chunk_neighbors_execute(ldbg,*chunk_neighbors,*grid,*amr,*amr_grid_pairs,*config,*chunk_neighbors_scratch,cs,cs_log2,*nbh_dist_lab, xform_filter, gpu_enabled, no_zorder, nbh_filter );
    }

  };

  template<class GridT> using PairTypeRCutNeighborsTmpl = PairTypeRCutNeighbors<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("pair_type_neighbors", make_grid_variant_operator< PairTypeRCutNeighborsTmpl > );
  }

}

