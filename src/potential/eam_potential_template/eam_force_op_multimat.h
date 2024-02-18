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
# define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)

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

    template<class RhoFieldAccT, bool NewtonSym>
    struct SymRhoOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      RhoFieldAccT m_rho_field = {};
      //exanb::ComputePairOptionalLocks<NewtonSym> & cp_locks;
      
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () ( Vec3d dr,double d2, int type_a, double& rho, CellParticlesT cells,size_t cell_b, size_t p_b, double /*scale*/ ) const
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
          rho += rholoc;
          if constexpr ( NewtonSym )
          {
            rholoc = 0.;
            drholoc = 0.;
            USTAMP_POTENTIAL_EAM_RHO( p, r, rholoc, drholoc, type_a, type_b );
            cells[cell_b][m_rho_field][p_b] += rholoc;
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
    };

    template<class EmbFieldAccT, bool NewtonSym>
    struct SymForceOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      EmbFieldAccT m_dEmb_field = {};

      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        Vec3d dr,double d2,
        double& ep, double& fx, double& fy, double& fz, int type_a /*, Mat3d& vir*/,
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
          const double fe_x = dr.x * fpair;
          const double fe_y = dr.y * fpair;
          const double fe_z = dr.z * fpair;

          fx += fe_x;
          fy += fe_y;
          fz += fe_z;

          ep += .5 * phi;
          //const Mat3d vir = tensor( Vec3d{fe_x,fe_y,fe_z}, Vec3d{dr.x,dr.y,dr.z} ) * -0.5;
          //virial += vir;
          if constexpr ( NewtonSym )
          {
            cells[cell_b][field::fx][p_b] -= fe_x;
            cells[cell_b][field::fy][p_b] -= fe_y;
            cells[cell_b][field::fz][p_b] -= fe_z;
            cells[cell_b][field::ep][p_b] += .5 * phi;
            //cells[cell_b][field::virial][p_b] += vir;
          }
        }
      }
      
    };



    struct EmbOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;

      template<class ComputePairBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,

        // particles fields to compute                
        double& ep,
        int type_a,
        double& dEmb,

        // data and locks accessors for neighbors (not used)
        CellParticlesT
        ) const
      {
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

        double rho = 0.;

#       pragma omp simd reduction(+:rho)
        for(size_t i=0;i<n;i++)
        {
          const double r = sqrt( tab.d2[i] );
          const int type_b = tab.nbh_pt[i][field::type];
          [[maybe_unused]] const bool pair_inversed = type_a > type_b;
          const int pair_id = unique_pair_id( type_a , type_b );
	  
          double rholoc = 0.;
          double drholoc = 0.;
          if( m_pair_enabled!=nullptr || m_pair_enabled[pair_id] )
          {
            USTAMP_POTENTIAL_EAM_RHO( p, r, rholoc, drholoc, type_b, type_a );
            rho += rholoc;
          }
        }

        double emb = 0.;
        /*double*/ dEmb = 0.;
        USTAMP_POTENTIAL_EAM_EMB( p, rho, emb, dEmb , type_a );
        //m_particle_emb[ m_cell_emb_offset[tab.cell] + tab.part ] = dEmb;	
        //	      if (rho > rhomax) emb += dEmb * (rho-rhomax);
        ep += emb;
      }

    };

    struct ForceOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      
      template<class ComputePairBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        int type_a,
        double dEmb,
        CellParticlesT
        ) const
      {
        FakeMat3d virial;
        (*this) ( n,tab,_ep,_fx,_fy,_fz, type_a, virial, dEmb, nullptr );
      }

      template<class ComputePairBufferT, class Mat3dT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        int type_a,
        Mat3dT& virial,
        double fpi, /* aka dEmb for central atom*/
        CellParticlesT
        ) const
      {
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

      	// derivative of Embedding function at atom i
        // const double fpi = m_particle_emb[ m_cell_emb_offset[tab.cell] + tab.part ];
        double ep=0., fx=0., fy=0., fz=0.;

        Mat3dT vir; // default constructor defines all elements to 0
        /* assert( vir.m11==0 && vir.m12==0 && vir.m13==0 &&
                vir.m21==0 && vir.m22==0 && vir.m23==0 &&
                vir.m31==0 && vir.m32==0 && vir.m33==0 ); */

#       pragma omp simd
        for(size_t i=0;i<n;i++)
        {
          tab.d2[i] = sqrt( tab.d2[i] );
        }

#       pragma omp simd reduction(+:fx,fy,fz,ep,vir)
        for(size_t i=0;i<n;i++)
        {
          const double r = tab.d2[i];
          const int type_b = tab.nbh_pt[i][field::type];
          //[[maybe_unused]] const bool pair_inversed = ( type_a > type_b );
          //const int pair_id = unique_pair_id(type_a,type_b);

	        // std::cout << "\t type_b        = " << type_b << std::endl;
	        // std::cout << "\t pair_inversed = " << pair_inversed << std::endl;
	        // std::cout << "\t pair_id       = " << pair_id << std::endl;
	        // std::cout << "\t rho n_b       = " << p[pair_id].m_parameters.n_a << std::endl;
	        // std::cout << "\t rho n_a       = " << p[pair_id].m_parameters.n_a << std::endl;	  
	  
          double Rho = 0.;
          double phi = 0.;
          double phip = 0.;
          double rhoip = 0.;
          double rhojp = 0.;	  
	        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
	        {
	          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhoip , type_a, type_b );
	          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, rhojp , type_b, type_a );
	          USTAMP_POTENTIAL_EAM_PHI( p, r, phi, phip  , type_a, type_b );
	        }
          double recip = 1.0/r;
	        double fpj = tab.nbh_pt[i][field::dEmb];
	        double psip = fpi * rhojp + fpj * rhoip + phip;
	        double fpair = psip*recip;
	  
          const double drx = tab.drx[i];
          const double dry = tab.dry[i];
          const double drz = tab.drz[i];
          const double fe_x = drx * fpair;
          const double fe_y = dry * fpair;
          const double fe_z = drz * fpair;

          fx  += fe_x;
          fy  += fe_y;
          fz  += fe_z;

          ep  += .5 * phi;
          vir += tensor( Vec3d{fe_x,fe_y,fe_z}, Vec3d{drx,dry,drz} ) * -0.5;
        }

        _ep += ep;
        _fx += fx;
        _fy += fy;
        _fz += fz;        
        virial += vir;
      }
    };

  }

}

namespace exanb
{
  template<> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::EmbOp>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif
  };

  template<> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::ForceOp>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif
  };

  template<class FieldAccT, bool NewtonSym> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::SymRhoOp<FieldAccT,NewtonSym> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible = true;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif
  };

  template<class FieldAccT, bool NewtonSym> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::SymForceOp<FieldAccT,NewtonSym> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible = true;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif
  };

  template<> struct ComputeCellParticlesTraits< exaStamp::PRIV_NAMESPACE_NAME::Rho2EmbOp >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

}




