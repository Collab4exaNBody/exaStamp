#pragma once

#include <onika/cuda/cuda.h>
#include <exanb/core/grid.h>
#include <onika/cuda/cuda_context.h>
//#include <onika/cuda/profiling.h>

#include "snapBs.h"
#include "snap_gpu_force_op.h"
#include "snap_gpu_context.h"

namespace SnapExt
{
  using namespace exanb;

  template<class CellsT, class OptionalT, class CPBufFactoryT, class ChunkSizeT, size_t BlockSize, int JMax, class FieldSetT >
  struct CudaSnap
  {

    ONIKA_HOST_DEVICE_FUNC static inline void compute_force(
      CellsT* cells,
      SnapKernelCounters * d_kernel_counters,
      const IJK& dims,
      unsigned int gl,
      const OptionalT& optional,
      const CPBufFactoryT& cpbuf_factory,
      const SnapGpuForceOp<BlockSize,JMax>& force_op,
      double rcut2,
      ChunkSizeT CS,
      FieldSetT
      )
    {
      constexpr onika::BoolConst<true> UseComputeBuffer {};

      const IJK dimsNoGL = { dims.i-2*gl , dims.j-2*gl , dims.k-2*gl };
      const uint64_t ncells_no_gl = dimsNoGL.i * dimsNoGL.j * dimsNoGL.k;
      ONIKA_CU_BLOCK_SHARED unsigned int cell_a_no_gl;
      do
      {
        if( ONIKA_CU_THREAD_IDX == 0 )
        {
          cell_a_no_gl = ONIKA_CU_ATOMIC_ADD( d_kernel_counters->cell_counter , 1u );
        }
        ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
        if( cell_a_no_gl < ncells_no_gl )
        {
          using FieldTuple = field_accessor_tuple_from_field_set_t< FieldSetT >;
          static constexpr size_t nb_fields = onika::tuple_size_const_v<FieldTuple>;
          static constexpr FieldTuple cp_fields = {};
          
          const IJK loc_a_no_gl = grid_index_to_ijk( dimsNoGL, cell_a_no_gl );
          const IJK loc_a = { loc_a_no_gl.i+gl , loc_a_no_gl.j+gl , loc_a_no_gl.k+gl };
          const size_t cell_a = grid_ijk_to_index( dims, loc_a );
          compute_cell_particle_pairs_cell(cells,dims,loc_a,cell_a,rcut2 ,cpbuf_factory, optional, force_op,
                                      cp_fields, CS, std::integral_constant<bool,false>{},
                                      DefaultPositionFields{} , UseComputeBuffer , std::make_index_sequence<nb_fields>{} );
        }
      }
      while( cell_a_no_gl < ncells_no_gl );
    }
  };

  template<class CellsT, class OptionalT, class CPBufFactoryT, class ChunkSizeT, size_t BlockSize, int JMax, class FieldSetT>
  ONIKA_DEVICE_KERNEL_FUNC void cuda_snap_force_kernel(
    CellsT* cells,
    SnapKernelCounters * d_kernel_counters,
    IJK dims,
    unsigned int gl,
    const OptionalT optional,
    const CPBufFactoryT cpbuf_factory,
    const SnapGpuForceOp<BlockSize,JMax> force_op,
    double rcut2,
    ChunkSizeT CS,
    FieldSetT)
  {
    CudaSnap<CellsT,OptionalT,CPBufFactoryT,ChunkSizeT,BlockSize,JMax,FieldSetT>::compute_force( cells, d_kernel_counters, dims, gl, optional, cpbuf_factory, force_op, rcut2, CS , FieldSetT{} );
  }

  template<class SnapGutContextT, class GridT, class OptionalArgsT, class FieldSetT>
  inline void cuda_snap_force( SnapGutContextT& ctx, snapBs& snap_bs, GridT& grid
                             , double rcut, double factor, double coef, double rfac0, double rmin0, bool bzflag
                             , bool ghost, const OptionalArgsT& optional
                             , FieldSetT )
  {
    //static constexpr onika::IntConst<1> const_1{};
    static constexpr onika::IntConst<4> const_4{};
    static constexpr onika::IntConst<8> const_8{};
    
    auto* cells = grid.cells();
    const IJK dims = grid.dimension();
    const unsigned int gl = ghost ? grid.ghost_layers() : 0;
    //const size_t N = grid.number_of_cells();
    const double rcut2 = rcut*rcut;

    // execution configuration
    static constexpr int BlockSize = CUDA_BLOCK_SIZE;
    static constexpr int JMax = SnapGutContextT::JMax;

    // reset cell counter to 0 before each run
    ctx.reset_cell_counters();

    auto cpbuf_factory = make_compute_pair_buffer<SnapComputeBuffer>();
    SnapConstantPointers constants = { ctx.n_fblocks, ctx.d_bs_fblock };
    SnapGpuForceOp<BlockSize,JMax> force_op { rcut, factor, coef, rfac0, rmin0 , ( bzflag ? snap_bs.en_zero_val() : 0.0 ) , constants , { ctx.d_block_scratch } , ctx.d_kernel_counters };
    const unsigned int cs = optional.nbh.m_chunk_size;
    onikaStream_t custr = ctx.custream; // ctx.m_threadStream[0];
/*
    static bool first_time=true;
    if(first_time)
    {
      first_time=false;
      using CellsT = std::remove_reference_t< decltype(cells[0]) >;
      using OptionalT = OptionalArgsT;
      using CPBufFactoryT = std::remove_reference_t< decltype(cpbuf_factory) >;
      using ChunkSizeT = onika::IntConst<8>;
      cudaOccupancyMaxActiveBlocksPerMultiprocessor( &mo , CudaSnap<CellsT,OptionalT,CPBufFactoryT,ChunkSizeT,BlockSize,JMax>::compute_force , BlockSize , 0 );
      std::cout << "max overlap detected = "<<mo<<"\n";
    }
*/
    std::cerr << "Deactivated GPU impelmentation, see file " << __FILE__ <<" at line " << __LINE__ << std::endl;
    std::abort();
    /*
    // need to be fixed match new function API
    switch( cs )
    {
      //case 1 : ONIKA_CU_LAUNCH_KERNEL( ctx.n_cu_blocks, BlockSize, 0, custr , cuda_snap_force_kernel, cells, ctx.d_kernel_counters, dims, gl, optional, cpbuf_factory, force_op, rcut2, const_1 , FieldSetT{} ); break;
      case 4 : ONIKA_CU_LAUNCH_KERNEL( ctx.n_cu_blocks, BlockSize, 0, custr , cuda_snap_force_kernel, cells, ctx.d_kernel_counters, dims, gl, optional, cpbuf_factory, force_op, rcut2, const_4 , FieldSetT{} ); break;
      case 8 : ONIKA_CU_LAUNCH_KERNEL( ctx.n_cu_blocks, BlockSize, 0, custr , cuda_snap_force_kernel, cells, ctx.d_kernel_counters, dims, gl, optional, cpbuf_factory, force_op, rcut2, const_8 , FieldSetT{} ); break;
      default: ONIKA_CU_LAUNCH_KERNEL( ctx.n_cu_blocks, BlockSize, 0, custr , cuda_snap_force_kernel, cells, ctx.d_kernel_counters, dims, gl, optional, cpbuf_factory, force_op, rcut2, cs      , FieldSetT{} ); break;
    }
    */
  }

}

