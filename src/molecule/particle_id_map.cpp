#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/particle_id_translation.h>

#include <random>
#include <mpi.h>

#include <vector>
#include <utility>

namespace exaStamp
{

/*
static inline void merge_id_maps( ParticleIdMap& target, ParticleIdMap&& source )
{
  target.merge( std::move(source) );
}

# pragma omp declare reduction(+ : ParticleIdMap : merge_id_maps(omp_out,std::move(omp_in)) ) initializer (omp_priv=ParticleIdMap{})
*/
  template<
      class GridT
    , class = AssertGridHasFields< GridT, field::_id >
    >
  class ParticleIdMapOperator : public OperatorNode
  {
    ADD_SLOT(GridT         , grid   , INPUT );
    ADD_SLOT(ParticleIdMap , id_map , INPUT_OUTPUT );

  public:

    // TODO: remove SLs that not connected to at least one particle in the central area (not ghost)

    inline void execute () override final
    {
      auto cells = grid->cells();     
      IJK dims = grid->dimension();
      
      id_map->clear();

      int nt = omp_get_max_threads();
      std::vector< ParticleIdMap > local_id_maps( nt );
      
#     pragma omp parallel
      {
        size_t tid = omp_get_thread_num();
        GRID_OMP_FOR_BEGIN(dims,i,_,nowait)
        {
          uint64_t const * __restrict__ part_ids = cells[i][field::id];
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            local_id_maps[tid].insert( std::pair<uint64_t,uint64_t>( part_ids[j] , encode_cell_particle(i,j) ) );
          }
        }
        GRID_OMP_FOR_END
      }

/*
      int step = 1;
      while( step < nt )
      {
        int participants = ( nt + step - 1 ) / ( 2 * step );
        assert( participants >= 1 );
#       pragma omp parallel num_threads(participants)
        {
          int tid = omp_get_thread_num() * step * 2;
          assert( tid % (step*2) == 0 );
          assert( tid >= 0 && tid < nt );
          int p = tid + step;
          if( p < nt )
          {
            //local_id_maps[tid].merge( std::move( local_id_maps[p] ) );
            local_id_maps[tid].insert( local_id_maps[p].begin() , local_id_maps[p].end() );
            local_id_maps[p].clear();
          }
        }
        step *= 2;
      }
      *id_map = std::move( local_id_maps[0] );
*/

      for(int i=0;i<nt;i++)
      {
        id_map->insert( local_id_maps[i].begin() , local_id_maps[i].end() );
        local_id_maps[i].clear();
      }

    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return
R"EOF(
Creates a map that tells where is (are) the particle instance(s) with a given unique id.
multiple instances can be found in ghost area because of periodic conditions.
)EOF";
    }

  };

  template<class GridT> using ParticleIdMapOp = ParticleIdMapOperator<GridT>;
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "particle_id_map", make_grid_variant_operator< ParticleIdMapOp > );
  }

}


#if 0
// multithreaded implementation
      int max_nt = omp_get_max_threads();
      ParticleIdMap local_id_map[max_nt];

#     pragma omp parallel
      {
        int nt = omp_get_num_threads();
        int tid = omp_get_thread_num();
        GRID_OMP_FOR_BEGIN(dims,i,_,nowait)
        GRID_OMP_FOR_BEGIN(dims,i,_)
        {
          uint64_t const * __restrict__ part_ids = cells[i][field::id];
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            uint64_t id = part_ids[j];
            local_id_map[tid].insert( std::pair<uint64_t,uint64_t>( id , encode_cell_particle(i,j) ) );
          }
        }
        GRID_OMP_FOR_END
        int round = 1;
        while( round < nt )
        {
#         pragma omp barrier
          if( (tid%round)==0 && (tid+round)<nt )
          {
            local_id_map[tid].merge( std::move(local_id_map[tid+round]) );
          }
          round = round << 1;
        }
      }
      *id_map = std::move( local_id_map[0] );
#endif


