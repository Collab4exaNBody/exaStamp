#pragma once

#include <cmath>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <onika/parallel/parallel_for.h>

#include <exaStamp/potential/eam/eam_buffer.h>
#include "potential.h"

#ifndef PRIV_NAMESPACE_NAME
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)
#endif

# define EamPotentialOperatorName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force)
# define EamPotentialFlatName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_flat_force)
# define EamParameterInitName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_init)
# define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)
# define EamPotentialFlatStr USTAMP_STR(EamPotentialFlatName)
# define EamParameterInitStr USTAMP_STR(EamParameterInitName)

#ifdef USTAMP_POTENTIAL_ENABLE_CUDA
#define USTAMP_POTENTIAL_CUDA_COMPATIBLE true
#else
#define USTAMP_POTENTIAL_CUDA_COMPATIBLE false
#endif

#ifndef USTAMP_POTENTIAL_EAM_MM_FORCE
#define USTAMP_POTENTIAL_EAM_MM_FORCE(p,dr_fe,phi,r,fpi,cells[cell_b][m_dEmb_field][p_b],type_a,type_b) \
{ \
  double Rho = 0.; \
  double phip = 0.; \
  double rhoip = 0.; \
  double rhojp = 0.; \
  USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhoip , type_a, type_b ); Rho=0.; \
  USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhojp , type_b, type_a ); Rho=0.; \
  USTAMP_POTENTIAL_EAM_PHI( p, r, phi, phip  , type_a, type_b ); \
  const double recip = 1.0/r; \
  const double fpj = cells[cell_b][m_dEmb_field][p_b]; \
  const double psip = fpi * rhojp + fpj * rhoip + phip; \
  const double fpair = psip*recip; \
  dr_fe *= fpair; \
}
#endif

#ifndef USTAMP_POTENTIAL_EAM_EMB_DERIV
# define USTAMP_POTENTIAL_EAM_EMB_DERIV USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_emb_deriv)
# ifdef USTAMP_POTENTIAL_ENABLE_CUDA
  ONIKA_HOST_DEVICE_FUNC
# endif
  static inline double USTAMP_POTENTIAL_EAM_EMB_DERIV ( const USTAMP_POTENTIAL_PARMS &p , double rho , int type_a )
  {
    double emb = 0.;
    double dEmb = 0.;
    USTAMP_POTENTIAL_EAM_EMB( p, rho, emb, dEmb , type_a );
    return dEmb;
  }
# endif

#ifndef USTAMP_POTENTIAL_EAM_RHO_NODERIV
# define USTAMP_POTENTIAL_EAM_RHO_NODERIV USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_rho_noderiv)
# ifdef USTAMP_POTENTIAL_ENABLE_CUDA
  ONIKA_HOST_DEVICE_FUNC
# endif
  static inline double USTAMP_POTENTIAL_EAM_RHO_NODERIV( const USTAMP_POTENTIAL_PARMS &p , double r, int type_a, int type_b )
  {
    double rho=0.0, drho=0.0;
    USTAMP_POTENTIAL_EAM_RHO( p, r, rho, drho, type_a, type_b );
    return rho;
  }
#endif

namespace exaStamp
{
  using namespace exanb;

  namespace PRIV_NAMESPACE_NAME
  {

    using EamMultiMatParams = EamMultimatParameters< USTAMP_POTENTIAL_PARMS >;
    using EamMultiMatParamsReadOnly = EamMultimatParametersRO< onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> >;

    struct RhoOpExtStorage
    {
      double rho = 0.0;
    };

    template<bool NewtonSym, bool Ghost, class ParticleLocksT>
    struct FlatSymRhoOp
    {
      using NeighborOffset = uint64_t;
      using ParticleIndex = uint32_t;
      using NeighborCount = uint16_t;
    
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> m_eam_parms = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      double rcut2 = 0.0;

      const NeighborOffset * __restrict__ m_neighbor_offset = nullptr;
      const ParticleIndex * __restrict__ m_neighbor_list = nullptr;
      const NeighborCount * __restrict__ m_half_count = nullptr;
      const uint64_t * __restrict__ m_ghost_flag = nullptr;

      double * __restrict__ m_rho = nullptr;
      const uint8_t * __restrict__ m_types = nullptr;
      const double * __restrict__ m_rx = nullptr;
      const double * __restrict__ m_ry = nullptr;
      const double * __restrict__ m_rz = nullptr;
      
      ParticleLocksT & m_locks;

      ONIKA_HOST_DEVICE_FUNC inline void operator () (size_t i) const
      {
        static constexpr bool CPAA = NewtonSym &&   onika::cuda::gpu_device_execution_t::value;
        static constexpr bool LOCK = NewtonSym && ! onika::cuda::gpu_device_execution_t::value;

        if constexpr ( ! Ghost ) if( ( m_ghost_flag[i/64] & ( 1ull << (i%64) ) ) != 0 ) return;

        const int type_a = m_types[i];
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

        const Vec3d ri = { m_rx[i] , m_ry[i] , m_rz[i] };
        const auto j_start = m_neighbor_offset[i];
        const unsigned int nj = NewtonSym ? m_half_count[i] : ( m_neighbor_offset[i+1] - j_start );
        const auto * __restrict__ neighbors = m_neighbor_list + j_start;
        
        double rho = 0.0;
        for(unsigned int jj=0;jj<nj;jj++)
        {
          const auto j = neighbors[jj];
          const Vec3d rj = { m_rx[j] , m_ry[j] , m_rz[j] };
          const Vec3d dr = rj - ri;
          const double r2 = norm2( dr );
          if( r2 < rcut2 )
          {
            const int type_b = m_types[j];  
            if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
            {
              const double r = sqrt( r2 );
              rho += USTAMP_POTENTIAL_EAM_RHO_NODERIV( m_eam_parms, r, type_b, type_a );
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
        if constexpr ( CPAA ) { ONIKA_CU_BLOCK_ATOMIC_ADD( m_rho[i] , rho ); }
        if constexpr (!CPAA ) { m_rho[i] += rho; }
        if constexpr ( LOCK ) m_locks[i].unlock();
      }
    };

    template<bool NewtonSym, class CPLocksT>
    struct SymRhoOp
    {      
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      const size_t * __restrict__ m_cell_particle_offset = nullptr;
      double * __restrict__ m_rho_data = nullptr;
      CPLocksT & cp_locks;
      
      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t , size_t , exanb::ComputePairParticleContextStart ) const
      {
        ctx.ext.rho = 0.0;
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
      {
        static constexpr bool CPAA = NewtonSym &&   onika::cuda::gpu_device_execution_t::value;
        static constexpr bool LOCK = NewtonSym && ! onika::cuda::gpu_device_execution_t::value;
        if constexpr ( LOCK ) cp_locks[cell_a][p_a].lock();
        if constexpr ( CPAA ) { ONIKA_CU_BLOCK_ATOMIC_ADD( m_rho_data[m_cell_particle_offset[cell_a]+p_a] , ctx.ext.rho ); }
        if constexpr (!CPAA ) { m_rho_data[m_cell_particle_offset[cell_a]+p_a] += ctx.ext.rho; }
        if constexpr ( LOCK ) cp_locks[cell_a][p_a].unlock();
      }
      
      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator() (
        ComputeBufferT& ctx, Vec3d dr,double d2, int type_a,
        CellParticlesT cells,size_t cell_b,size_t p_b , double /*scale */ ) const
      {
        static constexpr bool CPAA = NewtonSym &&   onika::cuda::gpu_device_execution_t::value;
        static constexpr bool LOCK = NewtonSym && ! onika::cuda::gpu_device_execution_t::value;
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;
        const double r = sqrt( d2 );
        //assert( cell_b < m_nb_cells );
        //assert( p_b < ( m_cell_offsets[cell_b+1] - m_cell_offsets[cell_b] ) );
        const int type_b = cells[cell_b][field::type][p_b];  
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
        {
          ctx.ext.rho += USTAMP_POTENTIAL_EAM_RHO_NODERIV( p, r, type_b, type_a );
          if constexpr ( NewtonSym )
          {
            const double rholoc = USTAMP_POTENTIAL_EAM_RHO_NODERIV( p, r, type_a, type_b );
            if constexpr ( LOCK ) cp_locks[cell_b][p_b].lock();
            if constexpr ( CPAA ) { ONIKA_CU_BLOCK_ATOMIC_ADD( m_rho_data[m_cell_particle_offset[cell_b]+p_b] , rholoc ); }
            if constexpr (!CPAA ) { m_rho_data[m_cell_particle_offset[cell_b]+p_b] += rholoc; }
            if constexpr ( LOCK ) cp_locks[cell_b][p_b].unlock();
          }
        }
      }
    };

    struct Rho2EmbOp
    {      
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () ( double& ep, int type_a, double& rho_dEmb ) const
      {
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_a)] )
        {
          double rho = rho_dEmb;
          double emb = 0.;
          double dEmb = 0.;
          USTAMP_POTENTIAL_EAM_EMB( p, rho, emb, dEmb , type_a );
          rho_dEmb = dEmb;
          ep += emb;
        }
      }
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () ( int type_a, double& rho_dEmb ) const
      {
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_a)] )
        {
          rho_dEmb = USTAMP_POTENTIAL_EAM_EMB_DERIV( p , rho_dEmb , type_a );
        }
      }
    };

    struct ForceOpExtStorage
    {
      double fpi = 0.0;
      Vec3d f = { 0. , 0. , 0. };
    };

    struct ForceEnergyOpExtStorage
    {
      double fpi = 0.0;
      Vec3d f = { 0. , 0. , 0. };
      double ep = 0.0;
    };

    template<bool NewtonSym, class CPLocksT>
    struct SymForceOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      const size_t * __restrict__ m_cell_particle_offset = nullptr;
      const double * __restrict__ m_dEmb_data = nullptr;
      CPLocksT& cp_locks;

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart ) const
      {
        ctx.ext.fpi = m_dEmb_data[m_cell_particle_offset[cell_a]+p_a];
        ctx.ext.f = Vec3d { 0. , 0. , 0. };
        if constexpr ( std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> )
        {
          ctx.ext.ep = 0.0;
        }
        //ctx.vir = Mat3d{};
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (const ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
      {
        static constexpr bool CPAA = NewtonSym &&   onika::cuda::gpu_device_execution_t::value;
        static constexpr bool LOCK = NewtonSym && ! onika::cuda::gpu_device_execution_t::value;
        if constexpr ( LOCK ) cp_locks[cell_a][p_a].lock();
        if constexpr ( CPAA )
        {
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::fx][p_a] , ctx.ext.f.x );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::fy][p_a] , ctx.ext.f.y );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::fz][p_a] , ctx.ext.f.z );
          if constexpr ( std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> )
          {
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::ep][p_a] , ctx.ext.ep );
          }
          //cells[cell_a][field::virial][p_a] += ctx.vir;
        }
        if constexpr ( !CPAA )
        {
          cells[cell_a][field::fx][p_a] += ctx.ext.f.x ;
          cells[cell_a][field::fy][p_a] += ctx.ext.f.y ;
          cells[cell_a][field::fz][p_a] += ctx.ext.f.z ;
          if constexpr ( std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> )
          {
            cells[cell_a][field::ep][p_a] += ctx.ext.ep ;
          }
          //cells[cell_a][field::virial][p_a] += ctx.vir;
        }
        if constexpr ( LOCK ) cp_locks[cell_a][p_a].unlock();
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
        ComputeBufferT& ctx,
        Vec3d dr_fe,double d2, int type_a, 
        CellParticlesT cells,size_t cell_b, size_t p_b
        , double /*scale*/) const
      {
        static constexpr bool CPAA = NewtonSym &&   onika::cuda::gpu_device_execution_t::value;
        static constexpr bool LOCK = NewtonSym && ! onika::cuda::gpu_device_execution_t::value;
        static constexpr bool eflag = std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> ;
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

        assert( p_b < cells[cell_b].size() );
        const int type_b = cells[cell_b][field::type][p_b];

        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
        {
          const double r = sqrt( d2 );
          assert( r > 0.0 );
          double phi = 0.;
          USTAMP_POTENTIAL_EAM_MM_FORCE(p,dr_fe,phi,r,ctx.ext.fpi,m_dEmb_data[m_cell_particle_offset[cell_b]+p_b],type_a,type_b);
          ctx.ext.f += dr_fe;
          if constexpr ( eflag ) ctx.ext.ep += .5 * phi;          
          if constexpr ( NewtonSym )
          {
            if constexpr ( LOCK ) cp_locks[cell_b][p_b].lock();
            if constexpr ( CPAA )
            {            
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fx][p_b] , - dr_fe.x );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fy][p_b] , - dr_fe.y );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fz][p_b] , - dr_fe.z );
              if constexpr ( eflag ) ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::ep][p_b] , .5 * phi );
            }
            if constexpr (!CPAA )
            {
              cells[cell_b][field::fx][p_b] -= dr_fe.x;
              cells[cell_b][field::fy][p_b] -= dr_fe.y;
              cells[cell_b][field::fz][p_b] -= dr_fe.z;
              if constexpr ( eflag ) cells[cell_b][field::ep][p_b] += .5 * phi;
            }
            if constexpr ( LOCK ) cp_locks[cell_b][p_b].unlock();
          }
        }
      }
      
    };

  }

}

namespace exanb
{

  template<bool NewtonSym, class CPLocksT> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::SymRhoOp<NewtonSym,CPLocksT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = USTAMP_POTENTIAL_CUDA_COMPATIBLE;
    static inline constexpr bool HasParticleContextStart = true;    
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
//    static inline constexpr bool RequiresNbhOptionalData = false; // interaction weighting
  };

  template<bool NewtonSym, class CPLocksT> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::SymForceOp<NewtonSym,CPLocksT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = USTAMP_POTENTIAL_CUDA_COMPATIBLE;
    static inline constexpr bool HasParticleContextStart = true;    
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
//    static inline constexpr bool RequiresNbhOptionalData = false; // interaction weighting
  };

  template<> struct ComputeCellParticlesTraits< exaStamp::PRIV_NAMESPACE_NAME::Rho2EmbOp >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = USTAMP_POTENTIAL_CUDA_COMPATIBLE;
  };
    
}

namespace onika
{
  namespace parallel
  {  
    template<bool NewtonSym, bool Ghost, class ParticleLocksT>
    struct ParallelForFunctorTraits< exaStamp::PRIV_NAMESPACE_NAME::FlatSymRhoOp<NewtonSym,Ghost,ParticleLocksT> >
    {      
      static inline constexpr bool CudaCompatible = true;
    };
  }
}


