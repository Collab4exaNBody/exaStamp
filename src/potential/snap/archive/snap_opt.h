#pragma once


//#include <iostream>
//#include <fstream>
//#include <iomanip> // setw.
#include <cmath> // sqrt.
//#include <strings.h>
//#include <algorithm> 
//#include <set> 
#include <cassert>
// #define SNAP_VERBOSE_DBG 1

#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_reduction.h>

#include "snap_ext.h"
#include "snap_3Dtypes.h"
#include "snap_gsh_opt.h"
 
#include "snap_constants.h"

namespace SnapExt
{


//#define SNAP_SHARED_DCMM_CACHE 1
//#define SNAP_MULTI_CMM_BUFFER 1

ONIKA_HOST_DEVICE_FUNC
inline void inplace_reduce_add_4d(double& a, double &b, double &c, double &d )
{
  static constexpr unsigned int QBLOCK_SIZE = SnapExt::CUDA_BLOCK_SIZE / 4;

  ONIKA_CU_BLOCK_SHARED double sdata[ SnapExt::CUDA_BLOCK_SIZE ];

  const unsigned int group_tid = ONIKA_CU_THREAD_IDX % QBLOCK_SIZE;
  const unsigned int group     = ONIKA_CU_THREAD_IDX / QBLOCK_SIZE;

//  ONIKA_FORCE_ASSERT(group<4);
//  ONIKA_FORCE_ASSERT(group_tid<QBLOCK_SIZE);

  if( group == 0 )
  {
    sdata[ 0*QBLOCK_SIZE + group_tid ] = a;
    sdata[ 1*QBLOCK_SIZE + group_tid ] = b;
    sdata[ 2*QBLOCK_SIZE + group_tid ] = c;
    sdata[ 3*QBLOCK_SIZE + group_tid ] = d;
  }
  ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
  for(unsigned int step=1;step<4;step++)
  {
    if( group == step )
    {
      sdata[ 0*QBLOCK_SIZE + group_tid ] += a;
      sdata[ 1*QBLOCK_SIZE + group_tid ] += b;
      sdata[ 2*QBLOCK_SIZE + group_tid ] += c;
      sdata[ 3*QBLOCK_SIZE + group_tid ] += d;
    }
    ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
  }

  for (unsigned int s=QBLOCK_SIZE/2; s>0; s>>=1)
  {
    if (group_tid < s)
    {
      sdata[ group*QBLOCK_SIZE + group_tid ] += sdata[ group*QBLOCK_SIZE + group_tid + s];
    }
    ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
  }
  a = sdata[ 0*QBLOCK_SIZE + 0 ];
  b = sdata[ 1*QBLOCK_SIZE + 0 ];
  c = sdata[ 2*QBLOCK_SIZE + 0 ];
  d = sdata[ 3*QBLOCK_SIZE + 0 ];
}

ONIKA_HOST_DEVICE_FUNC
inline void my_cu_memcpy(void* _dst, const void* _src, unsigned long _n)
{
  uint32_t* dst = (uint32_t*)_dst;
  const uint32_t* src = (const uint32_t*)_src;
  unsigned int n = _n / sizeof(uint32_t);
  ONIKA_CU_BLOCK_SIMD_FOR(unsigned int,i,0,n)
  {
    dst[i] = src[i];
  }
}

template<class CMMArray, class DCMMXArray, class DCMMYArray, class DCMMZArray, class SnapConstantPointersT>
ONIKA_HOST_DEVICE_FUNC
static inline void compute_all_lbs2(
    CMMArray cmm, DCMMXArray dcmm_x, DCMMYArray dcmm_y, DCMMZArray dcmm_z
  , const SnapConstantPointersT& cptrs
  , double& out_energy_contrib
  , double3d& out_force_contrib
  )
{
  double energy_contrib = 0.0;
  double3d force_contrib = { 0., 0., 0. };

  int last_bs = -1;
  Complexd lbs = { 0. , 0. };
  Complexd dlbs_x = { 0. , 0. };
  Complexd dlbs_y = { 0. , 0. };
  Complexd dlbs_z = { 0. , 0. };
  
  ONIKA_CU_BLOCK_SIMD_FOR(int,widx,0,cptrs.n_bs_fblock)
  //for(int widx=0;widx<cptrs.n_bs_fblock;widx++)
  {
    const auto wi = cptrs.bs_fblock[widx]; //.work_item[ONIKA_CU_THREAD_IDX];
    int bs = wi.bs; // if( bs >= 32768 ) bs -= 32768;
    if( bs != last_bs )
    {
      if( last_bs != -1 )
      {
        Complexd main_cmm = cmm[ last_bs ];
        Complexd main_dcmm_x = dcmm_x[ last_bs ];
        Complexd main_dcmm_y = dcmm_y[ last_bs ];
        Complexd main_dcmm_z = dcmm_z[ last_bs ];

        Complexd conj_cmm = conj( main_cmm );
        Complexd conj_dcmm_x = conj( main_dcmm_x );
        Complexd conj_dcmm_y = conj( main_dcmm_y );
        Complexd conj_dcmm_z = conj( main_dcmm_z );

        // add contribution to force and energy
        force_contrib.x += ( conj_dcmm_x * lbs + conj_cmm * dlbs_x ).real();
        force_contrib.y += ( conj_dcmm_y * lbs + conj_cmm * dlbs_y ).real();
        force_contrib.z += ( conj_dcmm_z * lbs + conj_cmm * dlbs_z ).real();
        energy_contrib  += (conj_cmm * lbs ).real();
      }
      
      last_bs = bs;
      lbs = Complexd{ 0. , 0. };
      dlbs_x = Complexd{ 0. , 0. };
      dlbs_y = Complexd{ 0. , 0. };
      dlbs_z = Complexd{ 0. , 0. };
    }

    auto cmm_a_cmm = cmm[ wi.idx_a ];
    auto cmm_a_dcmm_x = dcmm_x[ wi.idx_a ];
    auto cmm_a_dcmm_y = dcmm_y[ wi.idx_a ];
    auto cmm_a_dcmm_z = dcmm_z[ wi.idx_a ];
    
    auto cmm_b_cmm = cmm[ wi.idx_b ];
    auto cmm_b_dcmm_x = dcmm_x[ wi.idx_b ];
    auto cmm_b_dcmm_y = dcmm_y[ wi.idx_b ];
    auto cmm_b_dcmm_z = dcmm_z[ wi.idx_b ];
    
    // compute
    const double cg_prod = wi.cg_prod; //lbs_group.cg_prod[widx];
    lbs += cg_prod * cmm_a_cmm * cmm_b_cmm ;
    dlbs_x += cg_prod * ( cmm_a_dcmm_x * cmm_b_cmm + cmm_a_cmm * cmm_b_dcmm_x );
    dlbs_y += cg_prod * ( cmm_a_dcmm_y * cmm_b_cmm + cmm_a_cmm * cmm_b_dcmm_y );
    dlbs_z += cg_prod * ( cmm_a_dcmm_z * cmm_b_cmm + cmm_a_cmm * cmm_b_dcmm_z );
  }

  if( last_bs != -1 )
  {
    Complexd main_cmm = cmm[ last_bs ];
    Complexd main_dcmm_x = dcmm_x[ last_bs ];
    Complexd main_dcmm_y = dcmm_y[ last_bs ];
    Complexd main_dcmm_z = dcmm_z[ last_bs ];

    Complexd conj_cmm = conj( main_cmm );
    Complexd conj_dcmm_x = conj( main_dcmm_x );
    Complexd conj_dcmm_y = conj( main_dcmm_y );
    Complexd conj_dcmm_z = conj( main_dcmm_z );

    force_contrib.x += ( conj_dcmm_x * lbs + conj_cmm * dlbs_x ).real();
    force_contrib.y += ( conj_dcmm_y * lbs + conj_cmm * dlbs_y ).real();
    force_contrib.z += ( conj_dcmm_z * lbs + conj_cmm * dlbs_z ).real();
    energy_contrib += (conj_cmm * lbs ).real();
  }

  out_force_contrib -= force_contrib;
  out_energy_contrib += energy_contrib;
}

template<class SnapScratchBuffersT, int JMax>
ONIKA_HOST_DEVICE_FUNC
static inline double snap_force_opt(
                 double const * __restrict__ m_rx
               , double const * __restrict__ m_ry
               , double const * __restrict__ m_rz
               , const double * /* __restrict__ m_radius */
               , unsigned int m_N_atom
               , double rcut
               , double rfac0
               , double rmin0
               , double m_factor
               , double m_coefs_0
               //, const snapBsIdxCache* m_idx_cache
               , double3d * __restrict__ m_force
               , SnapExt::SnapTimerStats& stats
               , const SnapExt::SnapConstantPointers& cptrs
               , const SnapScratchBuffersT& scratch
               , onika::IntConst<JMax> JMaxCst
               )
{
  using namespace SnapExt;

  //static constexpr double rmin0 = 0.;
  
  using onika::cuda::block_reduce_add;
  using onika::cuda::UnitializedPlaceHolder;

  /******** constants depending on jmax ************/
  static constexpr auto compact_gsh_size = SnapExt::SnapConstants<JMax>::compact_gsh_size;
  static constexpr auto gsh_start = SnapExt::SnapConstants<JMax>::gsh_start;
  static constexpr auto m_two_jmax = SnapExt::SnapConstants<JMax>::m_two_jmax;
  static constexpr auto n_cmm_ones = SnapExt::SnapConstants<JMax>::n_cmm_ones;
  static constexpr auto cmm_init_ones = SnapExt::SnapConstants<JMax>::cmm_init_ones;
  //static constexpr auto n_bs = SnapExt::SnapConstants<JMax>::n_bs;
  /************************************************/
  
  static constexpr onika::IntConst<SnapExt::CUDA_BLOCK_SIZE> BlockSizeCst = {};

# ifdef SNAP_MULTI_CMM_BUFFER
  ONIKA_CU_BLOCK_SHARED UnitializedPlaceHolder<Complexd> _cmm[ compact_gsh_size ]; Complexd* __restrict__ buf_cmm = (Complexd*) _cmm;
  ONIKA_CU_BLOCK_SHARED UnitializedPlaceHolder<Complexd> _next_cmm[ compact_gsh_size ]; Complexd* __restrict__ buf_next_cmm = (Complexd*) _next_cmm;
# else
  ONIKA_CU_BLOCK_SHARED UnitializedPlaceHolder<Complexd> _cmm[ compact_gsh_size ]; Complexd * /*__restrict__*/  cmm = (Complexd*) _cmm;
# endif

  auto * /*__restrict__*/ dcmm_x = scratch.get_dcmm_x();
  auto * /*__restrict__*/ dcmm_y = scratch.get_dcmm_y();
  auto * /*__restrict__*/ dcmm_z = scratch.get_dcmm_z();

/*
  if( m_N_atom > SNAP_OPT_MAX_NEIGHBORS )
  {
    printf("N=%d MAX=%d\n",int(m_N_atom),int(SNAP_OPT_MAX_NEIGHBORS) );
  }
  ONIKA_FORCE_ASSERT( m_N_atom < SNAP_OPT_MAX_NEIGHBORS );
*/

  // ONIKA_CU_BLOCK_SYNC();

# ifdef SNAP_MULTI_CMM_BUFFER
  ONIKA_CU_BLOCK_SHARED int coop_N_atom;
  ONIKA_CU_BLOCK_SHARED int coop_NbhBufferSize;
  ONIKA_CU_BLOCK_SHARED int coop_BufferStart;
  ONIKA_CU_BLOCK_SHARED int coop_SkipSize;
  ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,compact_gsh_size) buf_next_cmm[k] = 0.;
  ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
  ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,n_cmm_ones) buf_next_cmm[cmm_init_ones.idx[k]] = 1.;
  if( ONIKA_CU_THREAD_IDX == 0 ) { coop_BufferStart=0; coop_N_atom=0; coop_NbhBufferSize=0; }
  ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
# else
# define coop_SkipSize 0
# define coop_NbhBufferSize coop_N_atom
# endif

  double m_energy = 0.;

  for(int coop_round=0;coop_round<SnapExt::CUDA_BLOCK_SIZE;coop_round++)
  {
#   ifndef SNAP_MULTI_CMM_BUFFER
    ONIKA_CU_BLOCK_SHARED int coop_N_atom;
#   else
    // swap cmm buffers
    { auto tmp_cmm = buf_next_cmm; buf_next_cmm = buf_cmm; buf_cmm = tmp_cmm; }
#   endif      

    // ========= void snapBs::compute_cmm( double rcut ) ================
    {
      ONIKA_CU_BLOCK_SHARED double coop_rcut;
      ONIKA_CU_BLOCK_SHARED double coop_factor;
      ONIKA_CU_BLOCK_SHARED double coop_rx[SNAP_CMM_COMPUTE_BUFFER_SIZE];
      ONIKA_CU_BLOCK_SHARED double coop_ry[SNAP_CMM_COMPUTE_BUFFER_SIZE];
      ONIKA_CU_BLOCK_SHARED double coop_rz[SNAP_CMM_COMPUTE_BUFFER_SIZE];
    
      //ONIKA_CU_BLOCK_SHARED double coop_radius[SNAP_OPT_MAX_NEIGHBORS];
      if( ONIKA_CU_THREAD_IDX==coop_round )
      {
#       ifdef SNAP_MULTI_CMM_BUFFER
        coop_SkipSize = coop_NbhBufferSize - coop_N_atom;
        coop_BufferStart = ( coop_BufferStart + coop_N_atom ) % SNAP_CMM_COMPUTE_BUFFER_SIZE;
        coop_NbhBufferSize = m_N_atom;
#       endif
        coop_N_atom = m_N_atom;
        coop_rcut = rcut;
        coop_factor = m_factor;
        for(int i=coop_SkipSize;i<m_N_atom;i++)
        {
          coop_rx[i] = m_rx[i];
          coop_ry[i] = m_ry[i];
          coop_rz[i] = m_rz[i];
        }
      }
      ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();

#     ifdef SNAP_MULTI_CMM_BUFFER
      if( ONIKA_CU_THREAD_IDX==(coop_round+1) && coop_SkipSize<coop_N_atom )
      {
        for(int i=0; coop_NbhBufferSize<SNAP_CMM_COMPUTE_BUFFER_SIZE && i<m_N_atom ; ++i, ++coop_NbhBufferSize)
        {
          coop_rx[coop_NbhBufferSize] = m_rx[i];
          coop_ry[coop_NbhBufferSize] = m_ry[i];
          coop_rz[coop_NbhBufferSize] = m_rz[i];
        }
      }
      ONIKA_CU_BLOCK_SYNC();
      if( coop_NbhBufferSize > coop_N_atom ) { ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,compact_gsh_size) buf_next_cmm[k] = 0.; }
      if( coop_SkipSize==0 && coop_N_atom>0 ) { ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,compact_gsh_size) buf_cmm[k] = 0.; }
      ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
      if( coop_NbhBufferSize > coop_N_atom ) { ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,n_cmm_ones) buf_next_cmm[cmm_init_ones.idx[k]] = 1.; }
      if( coop_SkipSize==0 && coop_N_atom>0 ) { ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,n_cmm_ones) buf_cmm[cmm_init_ones.idx[k]] = 1.; }  
      ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
#     else
      ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,compact_gsh_size) cmm[k] = 0.;
      ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
      ONIKA_CU_BLOCK_SIMD_FOR(int,k,0,n_cmm_ones) cmm[cmm_init_ones.idx[k]] = 1.;
      ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
#     endif

      SNAP_PROFILING_CODE( auto T0 = ONIKA_CU_CLOCK() );

      ONIKA_CU_BLOCK_SIMD_FOR_UNGUARDED(int,i,coop_SkipSize,coop_NbhBufferSize)
      {
        Complexd gsh[ compact_gsh_size + gsh_start ];
        double fcut = 0.0;
 
#       ifdef SNAP_MULTI_CMM_BUFFER
        const unsigned int cmm_idx = (coop_BufferStart+i) % SNAP_CMM_COMPUTE_BUFFER_SIZE;
        Complexd * __restrict__ cmm = ( i < coop_N_atom ) ? buf_cmm : buf_next_cmm;
#       else
#       define cmm_idx i
#       endif
 
        if( i < coop_NbhBufferSize )
        {
          complex3d dgsh[ compact_gsh_size + gsh_start ];
          double3d r_vec = { coop_rx[i] , coop_ry[i] , coop_rz[i] };
          const double radius_i = sqrt(r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z); //coop_radius[i]; // m_radius[i];
          snap_compute_gsh_ext(r_vec,radius_i,coop_rcut,rfac0,rmin0,onika::IntConst<m_two_jmax>{},gsh, dgsh );
          fcut = 0.5*(cos(M_PI*(radius_i-rmin0)/(coop_rcut-rmin0))+1.)*coop_factor;
          const double dfcut=-0.5*M_PI/(coop_rcut-rmin0)*sin(M_PI*(radius_i-rmin0)/(coop_rcut-rmin0))*coop_factor;
          r_vec = r_vec * dfcut / radius_i;
          for(int cmm_k=0;cmm_k<compact_gsh_size;cmm_k++)
          {
            int gsh_k = cmm_k + gsh_start;
            complex3d dc = r_vec * gsh[gsh_k] + fcut * dgsh[gsh_k];
            dcmm_x[ cmm_idx*compact_gsh_size + cmm_k ] = dc.x;
            dcmm_y[ cmm_idx*compact_gsh_size + cmm_k ] = dc.y;
            dcmm_z[ cmm_idx*compact_gsh_size + cmm_k ] = dc.z;
          }
        }

        if constexpr ( SnapExt::CUDA_BLOCK_SIZE <= compact_gsh_size )
        {
          //ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
          for(int cmm_k=0, rot_cmm_k=ONIKA_CU_THREAD_IDX; cmm_k<compact_gsh_size ; ++cmm_k)
          {
            int gsh_k = rot_cmm_k + gsh_start;
            Complexd cmm_contrib = { 0. , 0. };
            if( i < coop_NbhBufferSize ) cmm_contrib = gsh[gsh_k] * fcut;
            cmm[rot_cmm_k] += cmm_contrib;
            ++rot_cmm_k; if( rot_cmm_k==compact_gsh_size ) rot_cmm_k=0;
            ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
          }
        }
        if constexpr ( SnapExt::CUDA_BLOCK_SIZE > compact_gsh_size )
        {
          for(int cmm_k=0;cmm_k<compact_gsh_size;cmm_k++)
          {
            int gsh_k = cmm_k + gsh_start;
            Complexd cmm_contrib = { 0. , 0. };
            if( i < coop_NbhBufferSize ) cmm_contrib = gsh[gsh_k] * fcut;
            cmm_contrib = block_reduce_add( cmm_contrib , BlockSizeCst );
            if( ONIKA_CU_THREAD_IDX==0 ) cmm[cmm_k] += cmm_contrib;
          }
        }
        
      }
      ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();

      SNAP_PROFILING_CODE( stats.compute_cmm_time += ONIKA_CU_CLOCK()-T0 );

    } // scoped life for temporary shared variables


    // ===================== void snapBs::compute_bs() ==========================

    for (int i=0; i < coop_N_atom ; i++) // N_atom can be bound to a maximum constant value
    {
      SNAP_PROFILING_CODE( auto T0 = ONIKA_CU_CLOCK() );

#     ifdef SNAP_MULTI_CMM_BUFFER
      const unsigned int cmm_idx = (coop_BufferStart+i) % SNAP_CMM_COMPUTE_BUFFER_SIZE;
      const Complexd * __restrict__ cmm_ro = ( i < coop_N_atom ) ? buf_cmm : buf_next_cmm;
#     else
      const Complexd * /*__restrict__*/ cmm_ro = cmm;
#     endif

      //const int ai = (coop_base+i) % SNAP_OPT_MAX_NEIGHBORS;
      double3d f_val = { 0., 0., 0. };
      double energy_contrib = 0.0;

      // scoped shared variable to limit shared memory usage
      {
  #     ifdef SNAP_SHARED_DCMM_CACHE
  
        ONIKA_CU_BLOCK_SHARED UnitializedPlaceHolder<Complexd> _s_dcmm_x[compact_gsh_size];
        ONIKA_CU_BLOCK_SHARED UnitializedPlaceHolder<Complexd> _s_dcmm_y[compact_gsh_size];
        ONIKA_CU_BLOCK_SHARED UnitializedPlaceHolder<Complexd> _s_dcmm_z[compact_gsh_size];

        ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();

        my_cu_memcpy( _s_dcmm_x , dcmm_x+cmm_idx*compact_gsh_size , compact_gsh_size*sizeof(Complexd) );
        my_cu_memcpy( _s_dcmm_y , dcmm_y+cmm_idx*compact_gsh_size , compact_gsh_size*sizeof(Complexd) );
        my_cu_memcpy( _s_dcmm_z , dcmm_z+cmm_idx*compact_gsh_size , compact_gsh_size*sizeof(Complexd) );
        const Complexd * __restrict__ dcmm_x_ro = (const Complexd*) _s_dcmm_x;
        const Complexd * __restrict__ dcmm_y_ro = (const Complexd*) _s_dcmm_y;
        const Complexd * __restrict__ dcmm_z_ro = (const Complexd*) _s_dcmm_z;

  #     else
  
        const Complexd* /*__restrict__*/ dcmm_x_ro = dcmm_x + cmm_idx*compact_gsh_size;
        const Complexd* /*__restrict__*/ dcmm_y_ro = dcmm_y + cmm_idx*compact_gsh_size;
        const Complexd* /*__restrict__*/ dcmm_z_ro = dcmm_z + cmm_idx*compact_gsh_size;
        
  #     endif

        ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();

        SNAP_PROFILING_CODE( stats.cache_cmm_time += ONIKA_CU_CLOCK()-T0 );

        SNAP_PROFILING_CODE( T0=ONIKA_CU_CLOCK() );

        compute_all_lbs2(cmm_ro, dcmm_x_ro, dcmm_y_ro, dcmm_z_ro, cptrs, energy_contrib, f_val);        
      }
      
      SNAP_PROFILING_CODE( auto T1=ONIKA_CU_CLOCK(); stats.compute_lbs_time += T1-T0; T0=T1 );

      f_val = block_reduce_add( f_val , BlockSizeCst );
      energy_contrib = block_reduce_add( energy_contrib , BlockSizeCst );

      SNAP_PROFILING_CODE( T1=ONIKA_CU_CLOCK(); stats.reduction_time += T1-T0 );
       
      // we re-use shared rx to store output force
      if(ONIKA_CU_THREAD_IDX==coop_round)
      {
        m_energy += energy_contrib;
        m_force[i] = f_val;
      }

      ONIKA_CU_BLOCK_FENCE(); ONIKA_CU_BLOCK_SYNC();
    
    } // end of loop on neighbor atoms

  } // end of cooperation rounds loop

# ifdef SNAP_PROFILING_ENABLE
  stats.n_neighbors = m_N_atom;
  stats.n_atoms = ( m_N_atom > 0 ) ? 1 : 0;
# endif

  //if( ONIKA_CU_THREAD_IDX == 0 && ONIKA_CU_BLOCK_IDX==3 ) { printf("base=%d\n",coop_base); }

  m_energy *= m_factor / m_N_atom;
  m_energy += m_coefs_0;
  return m_energy;
}

} // end of namespace SnapExt

