#pragma once

#include <onika/integral_constant.h>
#include <onika/cuda/stl_adaptors.h>
#include <onika/force_assert.h>
#include <onika/parallel/parallel_execution_context.h>
#include <onika/cuda/cuda_context.h>
#include <onika/cuda/cuda_error.h>

#include "snap_ext.h"
#include "snap_block_scratch.h"

namespace SnapExt
{

  template<int _BlockSize, int _JMax>
  struct SnapGPUContext
  {
    static inline constexpr int BlockSize = _BlockSize;
    static inline constexpr int JMax = _JMax;
    
    BSFullBlockWorkItem* d_bs_fblock = nullptr;
    SnapBlockScratch<BlockSize,JMax>* d_block_scratch = nullptr;
    SnapKernelCounters * d_kernel_counters = nullptr;
    onikaStream_t custream = 0;
    unsigned int n_fblocks = 0;
    unsigned int n_cu_blocks = 0;
    snapBs* bound_bs_ctx = nullptr;

    std::vector< BSFullBlockWorkItem > fb_data;

    inline void initialize( onika::cuda::CudaContext& ctx , snapBs& _bs  )
    {
      if( bound_bs_ctx != &_bs || d_bs_fblock==nullptr || d_block_scratch==nullptr || d_kernel_counters==nullptr )
      {
        finalize();
        bound_bs_ctx = &_bs;
        
        custream = ctx.m_threadStream[0];
        
        n_cu_blocks = ctx.m_devices[0].m_deviceProp.multiProcessorCount * 4; //onika::parallel::ParallelExecutionContext::gpu_sm_mult();
        
        //std::vector< SnapExt::BSFullBlockWorkItem > fb_data;
        bound_bs_ctx->generate_lbs_compute_blocks( fb_data );
        n_fblocks = fb_data.size();
        
        size_t scratch_size = sizeof(SnapBlockScratch<BlockSize,JMax>) * n_cu_blocks;
        size_t fb_data_size = sizeof(BSFullBlockWorkItem)  * fb_data.size();
        size_t kernel_counters_size = sizeof(SnapKernelCounters);
	
        // std::cout << "SnapGpuContext @"<<(void*)this<<"\n";
        // std::cout << "Scratch: blocks="<<n_cu_blocks<<", bytes ="<<scratch_size<<std::endl;
        // std::cout << "LBS: BlockSize="<<BlockSize<<" , JMax="<<JMax<<" , n_fblocks="<<n_fblocks<<std::endl;

        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MALLOC(&d_block_scratch   , scratch_size ) );        
        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MALLOC(&d_bs_fblock       , fb_data_size  ) );            
        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MALLOC(&d_kernel_counters , kernel_counters_size ) );

        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMSET(d_block_scratch  , 0, scratch_size ) );
        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMSET(d_bs_fblock      , 0, fb_data_size ) );
        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMSET(d_kernel_counters, 0, kernel_counters_size ) );
        
        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMCPY_KIND(d_bs_fblock  , fb_data.data()  , fb_data_size , onikaMemcpyHostToDevice) );
      }
    }
    
    inline void reset_cell_counters()
    {
      ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMSET( d_kernel_counters, 0, sizeof(SnapKernelCounters) ) );
      //ONIKA_CU_CHECK_ERRORS( cudaMemsetAsync( d_kernel_counters, 0, sizeof(SnapExt::SnapKernelCounters) , custream ) );
    }

    inline void synchronize()
    {
      ONIKA_CU_CHECK_ERRORS( ONIKA_CU_STREAM_SYNCHRONIZE(custream) );
      
#     ifdef SNAP_PROFILING_ENABLE
      {
        SnapKernelCounters kernel_counters;
        ONIKA_CU_CHECK_ERRORS( cudaMemcpy(&kernel_counters, d_kernel_counters, sizeof(SnapKernelCounters), cudaMemcpyDeviceToHost) );
        const auto& stats = kernel_counters.stats;
        double scaling = 0.001 / stats.n_neighbors;
        std::cout <<"n_atoms="<<stats.n_atoms<<", n_nbh="<<stats.n_neighbors<<", reduce="<< stats.reduction_time*scaling<<", flush="<<stats.flush_time*scaling<<", LBS="<<stats.compute_lbs_time*scaling
                  <<", CMM="<<stats.compute_cmm_time*scaling<<", CMM_CACHE="<<stats.cache_cmm_time*scaling <<"\n";
      }
#     endif

    }

    inline void finalize()
    {
      // std::cout<<"Free SnapGPUContext\n";
      if( d_bs_fblock       != nullptr ) { ONIKA_CU_CHECK_ERRORS( ONIKA_CU_FREE(d_bs_fblock      ) ); d_bs_fblock      =nullptr; }
      if( d_block_scratch   != nullptr ) { ONIKA_CU_CHECK_ERRORS( ONIKA_CU_FREE(d_block_scratch  ) ); d_block_scratch  =nullptr; }
      if( d_kernel_counters != nullptr ) { ONIKA_CU_CHECK_ERRORS( ONIKA_CU_FREE(d_kernel_counters) ); d_kernel_counters=nullptr; }
      bound_bs_ctx = nullptr;
      n_fblocks = 0;
      n_cu_blocks = 0;
      custream = 0;
    }

    inline ~SnapGPUContext()
    {
      finalize();
    }

  };

}

