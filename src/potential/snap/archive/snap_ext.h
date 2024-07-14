#pragma once

#include "snap_3Dtypes.h"
#include <cstdint>
#include <onika/cuda/cuda.h>
#include "snap_constants.h"
//#define SNAP_PROFILING_ENABLE 1

#ifdef SNAP_PROFILING_ENABLE
#define SNAP_PROFILING_CODE(...) __VA_ARGS__
#else
#define SNAP_PROFILING_CODE(...) (void)0
#endif


namespace SnapExt
{
//  static inline constexpr unsigned int BLOCKS_PERS_SM = 6; // should be 6 // now controled through onika::task::ParallelTaskConfig::gpu_sm_mult()
  static inline constexpr unsigned int SNAP_OPT_MAX_NEIGHBORS = 320;
  static inline constexpr unsigned int CUDA_BLOCK_SIZE = 32; // should be 64 or more
  static inline constexpr unsigned int SNAP_CMM_COMPUTE_BUFFER_SIZE = CUDA_BLOCK_SIZE; //SNAP_OPT_MAX_NEIGHBORS;

  struct SnapTimerStats
  {
#   ifdef SNAP_PROFILING_ENABLE
    unsigned long long reduction_time = 0;
    unsigned long long flush_time = 0;
    unsigned long long compute_lbs_time = 0;
    unsigned long long compute_cmm_time = 0;
    unsigned long long cache_cmm_time = 0;
    unsigned long long n_neighbors = 0;
    unsigned long long n_atoms = 0;
#   endif
    ONIKA_HOST_DEVICE_FUNC inline SnapTimerStats& operator += (const SnapTimerStats& rhs)
    {
#   ifdef SNAP_PROFILING_ENABLE
      reduction_time   += rhs.reduction_time;
      flush_time       += rhs.flush_time;
      compute_lbs_time += rhs.compute_lbs_time;
      cache_cmm_time   += rhs.cache_cmm_time;
      n_neighbors      += rhs.n_neighbors;
      n_atoms          += rhs.n_atoms;
#   endif
      return *this;
    }
    ONIKA_HOST_DEVICE_FUNC inline void gather_all (SnapTimerStats& rhs)
    {
#   ifdef SNAP_PROFILING_ENABLE
      rhs.block_reduce();
      if( ONIKA_CU_THREAD_IDX == 0 )
      {
        ONIKA_CU_ATOMIC_ADD( reduction_time   , rhs.reduction_time );
        ONIKA_CU_ATOMIC_ADD( flush_time       , rhs.flush_time );
        ONIKA_CU_ATOMIC_ADD( compute_lbs_time , rhs.compute_lbs_time );
        ONIKA_CU_ATOMIC_ADD( compute_cmm_time , rhs.compute_cmm_time );
        ONIKA_CU_ATOMIC_ADD( cache_cmm_time   , rhs.cache_cmm_time );
        ONIKA_CU_ATOMIC_ADD( n_neighbors , rhs.n_neighbors );
        ONIKA_CU_ATOMIC_ADD( n_atoms     , rhs.n_atoms );
      }
#   endif
    }
    ONIKA_HOST_DEVICE_FUNC inline void block_reduce ()
    {
#   ifdef SNAP_PROFILING_ENABLE
      reduction_time   = onika::cuda::block_reduce_add( reduction_time   , onika::IntConst<CUDA_BLOCK_SIZE>{} );
      flush_time       = onika::cuda::block_reduce_add( flush_time       , onika::IntConst<CUDA_BLOCK_SIZE>{} );
      compute_lbs_time = onika::cuda::block_reduce_add( compute_lbs_time , onika::IntConst<CUDA_BLOCK_SIZE>{} );
      compute_cmm_time = onika::cuda::block_reduce_add( compute_cmm_time , onika::IntConst<CUDA_BLOCK_SIZE>{} );
      cache_cmm_time   = onika::cuda::block_reduce_add( cache_cmm_time   , onika::IntConst<CUDA_BLOCK_SIZE>{} );
      n_neighbors = onika::cuda::block_reduce_add( n_neighbors , onika::IntConst<CUDA_BLOCK_SIZE>{} );
      n_atoms     = onika::cuda::block_reduce_add( n_atoms , onika::IntConst<CUDA_BLOCK_SIZE>{} );
#   endif
    }
  };

  struct alignas(64) SnapKernelCounters
  {
    unsigned int cell_counter = 0;
    alignas(64) SnapTimerStats stats;
  };
  
  struct alignas(64) CMMData
  {
    Complexd cmm;
    Complexd dcmm_x;
    Complexd dcmm_y;
    Complexd dcmm_z;
  };

  template<class ComplexArray, class Complex3DArray>
  struct SnapScratchBuffers
  {
    ComplexArray cmm;
    ComplexArray dcmm_x;
    ComplexArray dcmm_y;
    ComplexArray dcmm_z;
    ComplexArray gsh;
    Complex3DArray dgsh;
    
    ONIKA_HOST_DEVICE_FUNC inline ComplexArray get_cmm() const { return cmm; }
    ONIKA_HOST_DEVICE_FUNC inline ComplexArray get_dcmm_x() const { return dcmm_x; }
    ONIKA_HOST_DEVICE_FUNC inline ComplexArray get_dcmm_y() const { return dcmm_y; }
    ONIKA_HOST_DEVICE_FUNC inline ComplexArray get_dcmm_z() const { return dcmm_z; }
//    ONIKA_HOST_DEVICE_FUNC inline ComplexArray get_gsh() const { return gsh; }
//    ONIKA_HOST_DEVICE_FUNC inline Complex3DArray get_dgsh() const { return dgsh; }
  };

  using SnapBsIndexT = uint16_t;

  struct SubIterT
  {
    SnapBsIndexT idx_a;
    SnapBsIndexT idx_b;
  };


  struct BSFullBlockWorkItem{ double cg_prod; uint16_t bs; uint16_t idx_a; uint16_t idx_b; uint16_t pad; };

  template<int BlockSize> struct LBSBlockConstants;

  struct SnapConstantPointers
  {
    const unsigned int n_bs_fblock = 0;
    const BSFullBlockWorkItem* __restrict__ bs_fblock = nullptr;
  };

  struct SinCosTheta
  {
    double sin_theta;
    double cos_theta;
  };

}

