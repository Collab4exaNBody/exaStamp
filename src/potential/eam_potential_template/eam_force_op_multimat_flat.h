#pragma once

#include <cmath>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <onika/parallel/parallel_for.h>

#include "eam_force_op_multimat.h"
#include "potential.h"

# define EamPotentialFlatName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_flat_force)
# define EamPotentialFlatStr USTAMP_STR(EamPotentialFlatName)

namespace exaStamp
{
  using namespace exanb;

  namespace PRIV_NAMESPACE_NAME
  {

    template<bool NewtonSym, class XFormT>
    struct FlatSymRhoOp
    {
      using NeighborOffset = uint64_t;
      using ParticleIndex = uint32_t;
      using NeighborCount = uint16_t;
    
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> m_eam_parms = {};
      double rcut2 = 0.0;

      const NeighborOffset * __restrict__ m_neighbor_offset = nullptr;
      const ParticleIndex * __restrict__ m_neighbor_list = nullptr;
      const NeighborCount * __restrict__ m_half_count = nullptr;

      const uint64_t * __restrict__ m_ghost_flag = nullptr;
      const uint8_t * __restrict__ m_pair_enabled = nullptr;

      double * __restrict__ m_rho = nullptr;
      const uint8_t * __restrict__ m_types = nullptr;
      const double * __restrict__ m_rx = nullptr;
      const double * __restrict__ m_ry = nullptr;
      const double * __restrict__ m_rz = nullptr;
      
      XFormT m_xform;
      spin_mutex * __restrict__ m_locks = nullptr;

      ONIKA_HOST_DEVICE_FUNC inline void operator () (size_t i) const
      {
        static constexpr bool CPAA = NewtonSym &&   onika::cuda::gpu_device_execution_t::value;
        static constexpr bool LOCK = NewtonSym && ! onika::cuda::gpu_device_execution_t::value;

        if ( m_ghost_flag != nullptr && ( m_ghost_flag[i/64] & ( 1ull << (i%64) ) ) != 0 ) return;

        const int type_a = m_types[i];
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

        const Vec3d ri = { m_rx[i] , m_ry[i] , m_rz[i] };
        const auto j_start = m_neighbor_offset[i];
        const unsigned int nj = NewtonSym ? m_half_count[i] : ( m_neighbor_offset[i+1] - j_start );
        const auto * __restrict__ neighbors = m_neighbor_list + j_start;
        
        double rho_i = 0.0;
        for(unsigned int jj=0;jj<nj;jj++)
        {
          const auto j = neighbors[jj];
          const Vec3d rj = { m_rx[j] , m_ry[j] , m_rz[j] };
          const Vec3d dr = m_xform.transformCoord( rj - ri );
          const double r2 = norm2( dr );
          if( r2 > 0.0 && r2 <= rcut2 )
          {
            const int type_b = m_types[j];  
            if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
            {
              const double r = sqrt( r2 );
              rho_i += USTAMP_POTENTIAL_EAM_RHO_NODERIV( m_eam_parms, r, type_b, type_a );
              if constexpr ( NewtonSym )
              {
                const double rho_j = USTAMP_POTENTIAL_EAM_RHO_NODERIV( m_eam_parms, r, type_a, type_b );
                if constexpr ( LOCK ) m_locks[j].lock();
                if constexpr ( CPAA ) { ONIKA_CU_BLOCK_ATOMIC_ADD( m_rho[j] , rho_j ); }
                if constexpr (!CPAA ) { m_rho[j] += rho_j; }
                if constexpr ( LOCK ) m_locks[j].unlock();
              }
            }
          }
        }
        if constexpr ( LOCK ) m_locks[i].lock();
        if constexpr ( CPAA ) { ONIKA_CU_BLOCK_ATOMIC_ADD( m_rho[i] , rho_i ); }
        if constexpr (!CPAA ) { m_rho[i] += rho_i; }
        if constexpr ( LOCK ) m_locks[i].unlock();
      }
    };
    
    template<bool EnergyFlag=true>
    struct FlatRho2EmbOp
    {      
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> m_eam_parms = {};

      const uint64_t * __restrict__ m_ghost_flag = nullptr;
      const uint8_t * __restrict__ m_pair_enabled = nullptr;

      const uint8_t * __restrict__ m_types = nullptr;
      double * __restrict__ m_rho_emb = nullptr;
      double * __restrict__ m_ep = nullptr;
      
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () ( size_t i ) const
      {
        if ( m_ghost_flag != nullptr && ( m_ghost_flag[i/64] & ( 1ull << (i%64) ) ) != 0 ) return;
        const unsigned int type_a = m_types[i];
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_a)] )
        {
          double emb = 0.;
          double dEmb = 0.;
          USTAMP_POTENTIAL_EAM_EMB( m_eam_parms, m_rho_emb[i], emb, dEmb , type_a );
          m_rho_emb[i] = dEmb;
          m_ep[i] += emb;
        }
      }
    };

    template<>
    struct FlatRho2EmbOp<false>
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> m_eam_parms = {};

      const uint64_t * __restrict__ m_ghost_flag = nullptr;
      const uint8_t * __restrict__ m_pair_enabled = nullptr;

      const uint8_t * __restrict__ m_types = nullptr;
      double * __restrict__ m_rho_emb = nullptr;
      
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () ( size_t i ) const
      {
        if ( m_ghost_flag != nullptr && ( m_ghost_flag[i/64] & ( 1ull << (i%64) ) ) != 0 ) return;
        const unsigned int type_a = m_types[i];
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_a)] )
        {
          m_rho_emb[i] = USTAMP_POTENTIAL_EAM_EMB_DERIV( m_eam_parms , m_rho_emb[i] , type_a );
        }
      }
    };
 
 
    template<bool NewtonSym, class XFormT>
    struct FlatSymForceOp
    {
      using NeighborOffset = uint64_t;
      using ParticleIndex = uint32_t;
      using NeighborCount = uint16_t;
    
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> m_eam_parms = {};
      double rcut2 = 0.0;

      const NeighborOffset * __restrict__ m_neighbor_offset = nullptr;
      const ParticleIndex * __restrict__ m_neighbor_list = nullptr;
      const NeighborCount * __restrict__ m_half_count = nullptr;

      const uint64_t * __restrict__ m_ghost_flag = nullptr;
      const uint8_t * __restrict__ m_pair_enabled = nullptr;

      const double * __restrict__ m_dEmb = nullptr;
      const uint8_t * __restrict__ m_types = nullptr;
      const double * __restrict__ m_rx = nullptr;
      const double * __restrict__ m_ry = nullptr;
      const double * __restrict__ m_rz = nullptr;

      double * __restrict__ m_fx = nullptr;
      double * __restrict__ m_fy = nullptr;
      double * __restrict__ m_fz = nullptr;
      double * __restrict__ m_ep = nullptr;

      XFormT m_xform;
      spin_mutex * __restrict__ m_locks = nullptr;

      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () ( size_t i ) const
      {
        static constexpr bool CPAA = NewtonSym &&   onika::cuda::gpu_device_execution_t::value;
        static constexpr bool LOCK = NewtonSym && ! onika::cuda::gpu_device_execution_t::value;

        if ( m_ghost_flag != nullptr && ( m_ghost_flag[i/64] & ( 1ull << (i%64) ) ) != 0 ) return;

        const int type_a = m_types[i];
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

        const Vec3d ri = { m_rx[i] , m_ry[i] , m_rz[i] };
        const auto j_start = m_neighbor_offset[i];
        const unsigned int nj = NewtonSym ? m_half_count[i] : ( m_neighbor_offset[i+1] - j_start );
        const auto * __restrict__ neighbors = m_neighbor_list + j_start;
        
        Vec3d force_i = { 0.0 , 0.0 , 0.0 };
        double ep_i = 0.0;
        for(unsigned int jj=0;jj<nj;jj++)
        {
          const auto j = neighbors[jj];
          const Vec3d rj = { m_rx[j] , m_ry[j] , m_rz[j] };
          Vec3d dr_fe = m_xform.transformCoord( rj - ri );
          const double r2 = norm2( dr_fe );
          if( r2 > 0.0 && r2 <= rcut2 )
          {
            const int type_b = m_types[j];  
            if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
            {
              const double r = sqrt( r2 );
              double phi = 0.;
              USTAMP_POTENTIAL_EAM_MM_FORCE( m_eam_parms, dr_fe, phi, r, m_dEmb[i], m_dEmb[j], type_a, type_b );
              force_i += dr_fe;
              ep_i += .5 * phi;          
              if constexpr ( NewtonSym )
              {
                if constexpr ( LOCK ) m_locks[j].lock();
                if constexpr ( CPAA )
                {
                  ONIKA_CU_BLOCK_ATOMIC_ADD( m_fx[j] , - dr_fe.x );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( m_fy[j] , - dr_fe.y );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( m_fz[j] , - dr_fe.z );
                  if( m_ep != nullptr ) ONIKA_CU_BLOCK_ATOMIC_ADD( m_ep[j] , .5 * phi );
                }
                if constexpr (!CPAA )
                {
                  m_fx[j] -= dr_fe.x;
                  m_fy[j] -= dr_fe.y;
                  m_fz[j] -= dr_fe.z;
                  if( m_ep != nullptr ) m_ep[j] += .5 * phi;
                }
                if constexpr ( LOCK ) m_locks[j].unlock();
              }
            }
          }
        }
        if constexpr ( LOCK ) m_locks[i].lock();
        if constexpr ( CPAA )
        {
          ONIKA_CU_BLOCK_ATOMIC_ADD( m_fx[i] , force_i.x );
          ONIKA_CU_BLOCK_ATOMIC_ADD( m_fy[i] , force_i.y );
          ONIKA_CU_BLOCK_ATOMIC_ADD( m_fz[i] , force_i.z );
          if( m_ep != nullptr ) ONIKA_CU_BLOCK_ATOMIC_ADD( m_ep[i] , ep_i );
        }
        if constexpr (!CPAA )
        {
          m_fx[i] += force_i.x;
          m_fy[i] += force_i.y;
          m_fz[i] += force_i.z;
          if( m_ep != nullptr ) m_ep[i] += ep_i;
        }
        if constexpr ( LOCK ) m_locks[i].unlock();
      }
      
    };


   
  }
}


namespace onika
{
  namespace parallel
  {  
    template<bool NewtonSym, class XFormT>
    struct ParallelForFunctorTraits< exaStamp::PRIV_NAMESPACE_NAME::FlatSymRhoOp<NewtonSym,XFormT> >
    {      
      static inline constexpr bool CudaCompatible = true;
    };
    
    template<bool EnergyFlag>
    struct ParallelForFunctorTraits< exaStamp::PRIV_NAMESPACE_NAME::FlatRho2EmbOp<EnergyFlag> >
    {      
      static inline constexpr bool CudaCompatible = true;
    };

    template<bool NewtonSym, class XFormT>
    struct ParallelForFunctorTraits< exaStamp::PRIV_NAMESPACE_NAME::FlatSymForceOp<NewtonSym,XFormT> >
    {      
      static inline constexpr bool CudaCompatible = true;
    };

  }
}


