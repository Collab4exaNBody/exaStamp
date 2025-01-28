//_enable_cuda // DO NOT REMOVE THIS LINE

#include <exanb/core/grid.h>
#include <onika/math/basic_types.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/domain.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <onika/cuda/cuda_math.h>

// this allows for parallel compilation of templated operator for each available field set

namespace exanb
{

  struct RemoveZeroDistParticlesFunctor
  {
    const double epsilon2 = 0.0;
    // ComputeBuffer less computation without virial
    template<class CellParticlesT, class NbhDataT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (const Vec3d& dr, double d2, uint64_t& id, CellParticlesT, size_t, size_t, const NbhDataT&) const
    {
      if( d2 <= epsilon2 ) id = onika::cuda::numeric_limits<uint64_t>::max;
    }
  };

  template<> struct ComputePairTraits<RemoveZeroDistParticlesFunctor>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible = true;
    static inline constexpr bool CudaCompatible = true;
  };

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_id >
    >
  class RemoveMergedAtoms : public OperatorNode
  {      
    ADD_SLOT( GridChunkNeighbors , chunk_neighbors , INPUT, REQUIRED, DocString{"neighbor list"} );
    ADD_SLOT( Domain             , domain          , INPUT );
    ADD_SLOT( double             , distance        , 1.e-3 );
    ADD_SLOT( GridT              , grid            , INPUT_OUTPUT );

  public:
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );      
      size_t n_cells = grid->number_of_cells();
      if( n_cells==0 )
      {
        return ;
      }
      
      const double mdist = *distance;
      
      ComputePairOptionalLocks<false> cp_locks {};
      // true means symmetric: we don't want to delete the pair of merged atoms, just one of the two
      exanb::GridChunkNeighborsLightWeightIt<true> nbh_it{ *chunk_neighbors };
      ComputePairBufferFactory< ComputePairBuffer2<> > force_buf = {};
      RemoveZeroDistParticlesFunctor force_op = { mdist*mdist };
      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
      compute_cell_particle_pairs(*grid, mdist, false, optional, force_buf, force_op, FieldSet<field::_id>{}, parallel_execution_context() );

      auto cells = grid->cells();
      size_t number_of_removed_particles = 0;
      
#     pragma omp parallel for schedule(dynamic) reduction(+:number_of_removed_particles)
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        unsigned int n_particles = cells[cell_i].size();
        for(unsigned int p=0;p<n_particles;)
        {
          if( cells[cell_i][field::id][p] == onika::cuda::numeric_limits<uint64_t>::max )
          {
            -- n_particles;
            cells[cell_i].write_tuple( p , cells[cell_i][n_particles] );
            ++ number_of_removed_particles;
          }
          else
          {
            ++ p;
          }
        }
        cells[cell_i].resize( n_particles , grid->cell_allocator() );
      }
      
      grid->rebuild_particle_offsets();
      
      ldbg << "number_of_removed_particles = " << number_of_removed_particles << std::endl;
    }

  };

  template<class GridT> using RemoveMergedAtomsTmpl = RemoveMergedAtoms<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(remove_merged_atoms)
  {  
    OperatorNodeFactory::instance()->register_factory( "remove_merged_atoms" , make_grid_variant_operator<RemoveMergedAtomsTmpl> );
  }

}


