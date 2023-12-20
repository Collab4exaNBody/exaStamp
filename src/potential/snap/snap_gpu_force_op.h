#pragma once

#include <exanb/compute/compute_pair_buffer.h>
#include <exanb/compute/compute_pair_traits.h>

#include <onika/integral_constant.h>
#include <onika/cuda/cuda_error.h>
#include <onika/force_assert.h>

#include "snap_block_scratch.h"
#include "snap_ext.h"
#include "snap_constants.h"
#include "snap_opt.h"

namespace SnapExt
{

  using namespace exanb;

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) SnapExtraStorage
  {
    // temporary arrays
    alignas(onika::memory::DEFAULT_ALIGNMENT) double3d force[SNAP_OPT_MAX_NEIGHBORS];
  };

  using SnapComputeBuffer = ComputePairBuffer2<false,false,SnapExtraStorage,DefaultComputePairBufferAppendFunc,SNAP_OPT_MAX_NEIGHBORS>;

  template<size_t BlockSize,int _JMax>
  struct SnapGpuForceOp
  {
    //using exanb::IJK;

    static inline constexpr int JMax = _JMax;
  
    const double rcut = 0.0;
    const double factor = 0.0;
    const double coefs_0 = 0.0;
    const double rfac0 = 0.0;
    const double rmin0 = 0.0;
    const double energy0 = 0.0;
    
    SnapConstantPointers constants;
    SnapBlockScratchAccessor<BlockSize,JMax> scratch;
    SnapKernelCounters * d_kernel_counters = nullptr;

    // without virial computation
    template<class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (
      size_t n,
      SnapComputeBuffer& tab,
      double& ep,
      double& fx,
      double& fy,
      double& fz,
      CellParticlesT* cells
      ) const
    {
      // for(unsigned int i=0;i<n;++i) { tab.d2[i] = sqrt( tab.d2[i] ); } // ont used in snap_force_jmax3
      SnapTimerStats stats;
      double E = snap_force_opt( tab.drx, tab.dry, tab.drz, nullptr/*tab.d2*/, n, rcut, rfac0, rmin0, factor, coefs_0, tab.ext.force, stats, constants, scratch, onika::IntConst<JMax>{} );
      d_kernel_counters->stats.gather_all(stats);

      E -= energy0;

      if( n > 0 )
      {
        double _fx = 0.;
        double _fy = 0.;
        double _fz = 0.;

        for(unsigned int i=0;i<n;++i)
        {
          const double3d F = tab.ext.force[i]; // snap_bs.force_val(i);
          _fx -= F.x;
          _fy -= F.y;
          _fz -= F.z;

          size_t cell_b=0, p_b=0;
          tab.nbh.get(i, cell_b, p_b);
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fx][p_b] , F.x );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fy][p_b] , F.y );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fz][p_b] , F.z );
        }

        ep += E ; //* 0.5;
        ONIKA_CU_BLOCK_ATOMIC_ADD( fx , _fx );
        ONIKA_CU_BLOCK_ATOMIC_ADD( fy , _fy );
        ONIKA_CU_BLOCK_ATOMIC_ADD( fz , _fz );
      }
    }

    // without virial computation
    template<class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (
      size_t n,
      SnapComputeBuffer& tab,
      double& ep,
      double& fx,
      double& fy,
      double& fz,
      Mat3d& virial,
      CellParticlesT* cells
      ) const
    {
      // for(unsigned int i=0;i<n;++i) { tab.d2[i] = sqrt( tab.d2[i] ); }  // ont used in snap_force_jmax3
      SnapTimerStats stats;
      double E = snap_force_opt( tab.drx, tab.dry, tab.drz, nullptr /*tab.d2*/, n, rcut, rfac0, rmin0, factor, coefs_0, tab.ext.force, stats, constants, scratch, onika::IntConst<JMax>{} );
      d_kernel_counters->stats.gather_all(stats);

      E -= energy0;

      if( n > 0 )
      {
        Mat3d _vir;
        double _fx = 0.;
        double _fy = 0.;
        double _fz = 0.;

        for(unsigned int i=0;i<n;++i)
        {
          const double3d F = tab.ext.force[i]; // snap_bs.force_val(i);
          _fx -= F.x;
          _fy -= F.y;
          _fz -= F.z;
          auto v_contrib = tensor( Vec3d{F.x,F.y,F.z}, Vec3d{tab.drx[i],tab.dry[i],tab.drz[i]} ) ; //* -0.5;
          _vir += v_contrib;

          size_t cell_b=0, p_b=0;
          tab.nbh.get(i, cell_b, p_b);
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fx][p_b] , F.x );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fy][p_b] , F.y );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fz][p_b] , F.z );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m11 , - v_contrib.m11 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m12 , - v_contrib.m12 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m13 , - v_contrib.m13 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m21 , - v_contrib.m21 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m22 , - v_contrib.m22 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m23 , - v_contrib.m23 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m31 , - v_contrib.m31 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m32 , - v_contrib.m32 );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::virial][p_b].m33 , - v_contrib.m33 );
        }

        ep += E ; //* 0.5;
        ONIKA_CU_BLOCK_ATOMIC_ADD( fx , _fx );
        ONIKA_CU_BLOCK_ATOMIC_ADD( fy , _fy );
        ONIKA_CU_BLOCK_ATOMIC_ADD( fz , _fz );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m11 , _vir.m11 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m12 , _vir.m12 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m13 , _vir.m13 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m21 , _vir.m21 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m22 , _vir.m22 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m23 , _vir.m23 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m31 , _vir.m31 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m32 , _vir.m32 );
        ONIKA_CU_BLOCK_ATOMIC_ADD( virial.m33 , _vir.m33 );
      }
    }

  };

}

namespace exanb
{
  // partial specialization to require synchronous functor call
  template<size_t BlockSize, int JMax>
  struct ComputePairTraits< SnapExt::SnapGpuForceOp<BlockSize,JMax> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = true;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
  };
}

