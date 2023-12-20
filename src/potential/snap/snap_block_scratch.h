#pragma once

#include "snap_3Dtypes.h"
#include "snap_constants.h"
#include <cstdint>
#include <onika/cuda/cuda.h>


namespace SnapExt
{

  struct Complex2DArrayAccessor
  {
    double & r;
    double & i;
    ONIKA_HOST_DEVICE_FUNC inline operator Complexd () const { return {r,i}; }
    ONIKA_HOST_DEVICE_FUNC inline Complex2DArrayAccessor& operator = (const Complexd& c) { r=c.r; i=c.i; return *this; }
    ONIKA_HOST_DEVICE_FUNC inline Complex2DArrayAccessor& operator += (const Complexd& c) { r+=c.r; i+=c.i; return *this; }
  };

  template<unsigned int ArraySize>
  struct alignas(64) Complex2DArray
  {
    alignas(64) double r[ArraySize];
    alignas(64) double i[ArraySize];
    ONIKA_HOST_DEVICE_FUNC inline Complexd operator [] (unsigned int idx) const { return { r[idx] , i[idx] }; }
    ONIKA_HOST_DEVICE_FUNC inline Complex2DArrayAccessor operator [] (unsigned int idx) { return { r[idx] , i[idx] }; }
  };

  template<size_t BlockSize, int JMax>
  struct alignas(64) SnapBlockScratch
  {
    static inline constexpr size_t ArraySize = SNAP_CMM_COMPUTE_BUFFER_SIZE * SnapConstants<JMax>::compact_gsh_size;
    alignas(64) Complexd dcmm_x[ ArraySize ];
    alignas(64) Complexd dcmm_y[ ArraySize ];
    alignas(64) Complexd dcmm_z[ ArraySize ];    
    /*
    Complex2DArray<ArraySize> dcmm_x;
    Complex2DArray<ArraySize> dcmm_y;
    Complex2DArray<ArraySize> dcmm_z;
    */
  };

  template<size_t BlockSize, int JMax>
  struct SnapBlockScratchAccessor
  {
    static inline constexpr size_t ArraySize = SnapBlockScratch<BlockSize,JMax>::ArraySize;    
    SnapBlockScratch<BlockSize,JMax>* __restrict__ blocks;
    ONIKA_HOST_DEVICE_FUNC inline Complexd* get_dcmm_x() const { return blocks[ONIKA_CU_BLOCK_IDX].dcmm_x; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd* get_dcmm_y() const { return blocks[ONIKA_CU_BLOCK_IDX].dcmm_y; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd* get_dcmm_z() const { return blocks[ONIKA_CU_BLOCK_IDX].dcmm_z; }
/*
    ONIKA_HOST_DEVICE_FUNC inline Complex2DArray<ArraySize>& get_dcmm_x() const { return blocks[ONIKA_CU_BLOCK_IDX].dcmm_x; }
    ONIKA_HOST_DEVICE_FUNC inline Complex2DArray<ArraySize>& get_dcmm_y() const { return blocks[ONIKA_CU_BLOCK_IDX].dcmm_y; }
    ONIKA_HOST_DEVICE_FUNC inline Complex2DArray<ArraySize>& get_dcmm_z() const { return blocks[ONIKA_CU_BLOCK_IDX].dcmm_z; }
*/
  };

}

