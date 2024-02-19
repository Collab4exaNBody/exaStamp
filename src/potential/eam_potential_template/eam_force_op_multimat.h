#pragma once

#include <cmath>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>

#include <exaStamp/potential/eam/eam_buffer.h>
#include "potential.h"

#ifndef PRIV_NAMESPACE_NAME
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)
#endif

# define EamPotentialOperatorName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force)
# define EamParameterInitName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_init)
# define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)
# define EamParameterInitStr USTAMP_STR(EamParameterInitName)

#ifdef USTAMP_POTENTIAL_ENABLE_CUDA
#define USTAMP_POTENTIAL_CUDA_COMPATIBLE true
#else
#define USTAMP_POTENTIAL_CUDA_COMPATIBLE false
#endif

namespace exaStamp
{
  using namespace exanb;

  namespace PRIV_NAMESPACE_NAME
  {

    /*
    notes:
      - emb/rho same array
      - Warning about compute_rho&co functions sets or adds quantity to rho/drho (resp. phi/dphi, emd/dEmb)
      - try buffer less version for force op
    */

    using EamMultiMatParams = EamMultimatParameters< USTAMP_POTENTIAL_PARMS >;
    using EamMultiMatParamsReadOnly = EamMultimatParametersRO< onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> >;

    struct RhoOpExtStorage
    {
      double rho = 0.0;
    };

    template<class RhoFieldAccT, bool NewtonSym, class CPLocksT, bool CPAA>
    struct SymRhoOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      RhoFieldAccT m_rho_field = {};
      CPLocksT & cp_locks;
      
      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t , size_t , exanb::ComputePairParticleContextStart ) const
      {
        ctx.ext.rho = 0.0;
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
      {
        cp_locks[cell_a][p_a].lock();
        if constexpr ( CPAA ) { ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_rho_field][p_a] , ctx.ext.rho ); }
        else { cells[cell_a][m_rho_field][p_a] += ctx.ext.rho; }
        cp_locks[cell_a][p_a].unlock();
      }
      
      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, Vec3d dr,double d2, int type_a, double& rho, CellParticlesT cells,size_t cell_b, size_t p_b, double /*scale*/ ) const
      {
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;
        const double r = sqrt( d2 );
        assert( cell_b < m_nb_cells );
        assert( p_b < ( m_cell_offsets[cell_b+1] - m_cell_offsets[cell_b] ) );
        const int type_b = cells[cell_b][field::type][p_b];  
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
        {
          double rholoc = 0.;
          double drholoc = 0.;
          USTAMP_POTENTIAL_EAM_RHO( p, r, rholoc, drholoc, type_b, type_a );
          ctx.ext.rho += rholoc;
          if constexpr ( NewtonSym )
          {
            rholoc = 0.;
            drholoc = 0.;
            USTAMP_POTENTIAL_EAM_RHO( p, r, rholoc, drholoc, type_a, type_b );
            cp_locks[cell_b][p_b].lock();
            if constexpr ( CPAA ) { ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_rho_field][p_b] , rholoc ); }
            else { cells[cell_b][m_rho_field][p_b] += rholoc; }
            cp_locks[cell_b][p_b].unlock();
          }
        }
      }
    };

    struct Rho2EmbOp
    {      
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      ONIKA_HOST_DEVICE_FUNC inline void operator () ( double& ep, int type_a, double rho, double& dEmb ) const
      {
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_a)] )
        {
          double emb = 0.;
          /*double*/ dEmb = 0.;
          USTAMP_POTENTIAL_EAM_EMB( p, rho, emb, dEmb , type_a );
          ep += emb;
        }
      }
      ONIKA_HOST_DEVICE_FUNC inline void operator () ( int type_a, double rho, double& dEmb ) const
      {
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_a)] )
        {
          double emb = 0.;
          /*double*/ dEmb = 0.;
          USTAMP_POTENTIAL_EAM_EMB( p, rho, emb, dEmb , type_a );
        }
      }
    };

    struct ForceOpExtStorage
    {
      Vec3d f = { 0. , 0. , 0. };
    };

    struct ForceEnergyOpExtStorage
    {
      Vec3d f = { 0. , 0. , 0. };
      double ep = 0.0;
    };

    template<class EmbFieldAccT, bool NewtonSym, class CPLocksT, bool CPAA>
    struct SymForceOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      EmbFieldAccT m_dEmb_field = {};
      CPLocksT& cp_locks;

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart ) const
      {
        ctx.ext.f = Vec3d { 0. , 0. , 0. };
        if constexpr ( std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> )
        {
          ctx.ext.ep = 0.0;
        }
        //ctx.vir = Mat3d{};
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (const ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
      {
        cp_locks[cell_a][p_a].lock();
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
        else
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
        cp_locks[cell_a][p_a].unlock();
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        ComputeBufferT& ctx,
        Vec3d dr,double d2,
        double& fx, double& fy, double& fz, int type_a,
        double fpi, 
        CellParticlesT cells,size_t cell_b, size_t p_b,
        double /*scale*/) const
      {
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;
        
        assert( p_b < cells[cell_b].size() );
        const int type_b = cells[cell_b][field::type][p_b];
        
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
        {
          const double r = sqrt( d2 );
          assert( r > 0.0 );

#         ifdef USTAMP_POTENTIAL_EAM_MM_FORCE
          const Vec3d fe = USTAMP_POTENTIAL_EAM_MM_FORCE(p,dr,r,fpi,cells[cell_b][m_dEmb_field][p_b],type_a,type_b);
#         else
          double Rho = 0.;
          double phi = 0.;
          double phip = 0.;
          double rhoip = 0.;
          double rhojp = 0.;
          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhoip , type_a, type_b ); Rho=0.;
          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhojp , type_b, type_a ); Rho=0.;
          USTAMP_POTENTIAL_EAM_PHI( p, r, phi, phip  , type_a, type_b ); phi=0.;
          const double recip = 1.0/r;
	        const double fpj = cells[cell_b][m_dEmb_field][p_b];
	        const double psip = fpi * rhojp + fpj * rhoip + phip;
	        const double fpair = psip*recip;
	        const Vec3d fe = { dr.x * fpair , dr.y * fpair , dr.z * fpair };
#         endif

          ctx.ext.f += fe;
          
          if constexpr ( NewtonSym )
          {
            cp_locks[cell_b][p_b].lock();
            if constexpr ( CPAA )
            {            
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fx][p_b] , -fe.x );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fy][p_b] , -fe.y );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fz][p_b] , -fe.z );
            }
            else
            {
              cells[cell_b][field::fx][p_b] -= fe.x;
              cells[cell_b][field::fy][p_b] -= fe.y;
              cells[cell_b][field::fz][p_b] -= fe.z;
            }
            cp_locks[cell_b][p_b].unlock();
          }
        }
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        ComputeBufferT& ctx,
        Vec3d dr,double d2,
        double& fx, double& fy, double& fz, double& ep, int type_a,
        double fpi, 
        CellParticlesT cells,size_t cell_b, size_t p_b,
        double /*scale*/) const
      {
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;
        
        assert( p_b < cells[cell_b].size() );
        const int type_b = cells[cell_b][field::type][p_b];
        
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
        {
          const double r = sqrt( d2 );
          assert( r > 0.0 );

          double Rho = 0.;
          double phi = 0.;
          double phip = 0.;
          double rhoip = 0.;
          double rhojp = 0.;	  

          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhoip , type_a, type_b ); Rho=0.0;
          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhojp , type_b, type_a ); Rho=0.0;
          USTAMP_POTENTIAL_EAM_PHI( p, r, phi, phip  , type_a, type_b );

          const double recip = 1.0/r;
	        const double fpj = cells[cell_b][m_dEmb_field][p_b]; //tab.nbh_pt[i][field::dEmb];
	        const double psip = fpi * rhojp + fpj * rhoip + phip;
	        const double fpair = psip*recip;
	        const Vec3d fe = { dr.x * fpair , dr.y * fpair , dr.z * fpair };

          ctx.ext.f += fe;
          ctx.ext.ep += .5 * phi;
          //const Mat3d vir = tensor( Vec3d{fe_x,fe_y,fe_z}, Vec3d{dr.x,dr.y,dr.z} ) * -0.5;
          //virial += vir;
          
          if constexpr ( NewtonSym )
          {
            cp_locks[cell_b][p_b].lock();
            if constexpr ( CPAA )
            {            
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fx][p_b] , -fe.x );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fy][p_b] , -fe.y );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fz][p_b] , -fe.z );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::ep][p_b] , .5 * phi );
              //cells[cell_b][field::virial][p_b] += vir;
            }
            else
            {
              cells[cell_b][field::fx][p_b] -= fe.x;
              cells[cell_b][field::fy][p_b] -= fe.y;
              cells[cell_b][field::fz][p_b] -= fe.z;
              cells[cell_b][field::ep][p_b] += .5 * phi;
              //cells[cell_b][field::virial][p_b] += vir;
            }
            cp_locks[cell_b][p_b].unlock();
          }
        }
      }
      
    };

  }

}

namespace exanb
{

  template<class FieldAccT, bool NewtonSym, class CPLocksT, bool CPAA> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::SymRhoOp<FieldAccT,NewtonSym,CPLocksT,CPAA> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = USTAMP_POTENTIAL_CUDA_COMPATIBLE;
    static inline constexpr bool HasParticleContextStart = true;    
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };

  template<class FieldAccT, bool NewtonSym, class CPLocksT, bool CPAA> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::SymForceOp<FieldAccT,NewtonSym,CPLocksT,CPAA> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = USTAMP_POTENTIAL_CUDA_COMPATIBLE;
    static inline constexpr bool HasParticleContextStart = true;    
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };

  template<> struct ComputeCellParticlesTraits< exaStamp::PRIV_NAMESPACE_NAME::Rho2EmbOp >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = USTAMP_POTENTIAL_CUDA_COMPATIBLE;
  };

}




