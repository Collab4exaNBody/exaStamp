#pragma once

#include <cmath>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_pair_traits.h>

#include "potential.h"

#ifndef PRIV_NAMESPACE_NAME
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)
#endif

namespace exaStamp
{
  using namespace exanb;

  namespace PRIV_NAMESPACE_NAME
  {

    struct EmbOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p;
      const double rhoCut;
      const double phiCut;
      // const size_t* m_cell_emb_offset = nullptr;
      // double* m_particle_emb = nullptr;

      template<class ComputePairBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,

        // particles fields to compute                
        double& ep,
        double& dEmb,

        // data and locks accessors for neighbors (not used)
        CellParticlesT
        ) const
      {
        double particle_rho = 0.;

#       ifndef XSTAMP_EAM_CHECK_PAIR_DISTANCE
#       pragma omp simd reduction(+:particle_rho)
#       endif
        for(size_t i=0;i<n;i++)
        {
          double r = sqrt( tab.d2[i] );
#         ifdef XSTAMP_EAM_CHECK_PAIR_DISTANCE
          if( r < XSTAMP_EAM_CHECK_PAIR_DISTANCE )
          {
            printf("EAM particle distance limit reached : r=%.6e , below limit %.6e\n",r,XSTAMP_EAM_CHECK_PAIR_DISTANCE);
            ONIKA_CU_ABORT();
          }
#         endif
          double Rho = 0.;
          double dRho = 0.;
          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, dRho );
          Rho -= rhoCut;      
          particle_rho += Rho;
        }

        double Emb = 0.;
        /*double*/ dEmb = 0.;
        USTAMP_POTENTIAL_EAM_EMB( p, particle_rho, Emb, dEmb );
        ep += Emb;
        // m_particle_emb[ m_cell_emb_offset[tab.cell] + tab.part ] = dEmb;
      }

    };

    struct ForceOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p;
      const double rhoCut;
      const double phiCut;
      //const size_t* m_cell_emb_offset = nullptr;
      //const double* m_particle_emb = nullptr;

      template<class ComputePairBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        double dEmb,
        CellParticlesT 
        ) const
      {
        FakeMat3d virial;
        (*this) ( n,tab,_ep,_fx,_fy,_fz, virial, dEmb, nullptr );
      }

      template<class ComputePairBufferT, class Mat3dT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        Mat3dT& virial,
        double dEmb,
        CellParticlesT
        ) const
      {
        //const double emb = m_particle_emb[ m_cell_emb_offset[tab.cell] + tab.part ];
        double ep=0., fx=0., fy=0., fz=0.;

        Mat3dT vir; // default constructor defines all elements to 0
        /* assert( vir.m11==0 && vir.m12==0 && vir.m13==0 &&
                vir.m21==0 && vir.m22==0 && vir.m23==0 &&
                vir.m31==0 && vir.m32==0 && vir.m33==0 ); */

#       ifndef XSTAMP_EAM_CHECK_PAIR_DISTANCE
#       pragma omp simd reduction(+:fx,fy,fz,ep,vir)
#       endif
        for(size_t i=0;i<n;i++)
        {
          double r = sqrt( tab.d2[i] );
#         ifdef XSTAMP_EAM_CHECK_PAIR_DISTANCE
          if( r < XSTAMP_EAM_CHECK_PAIR_DISTANCE )
          {
            printf("EAM particle distance limit reached : r=%.6e , below limit %.6e\n",r,XSTAMP_EAM_CHECK_PAIR_DISTANCE);
            ONIKA_CU_ABORT();
          }
#         endif
          double Rho = 0.;
          double dRho = 0.;
          double Phi = 0.;
          double dPhi = 0.;
          USTAMP_POTENTIAL_EAM_RHO( p, r, Rho, dRho );          
          USTAMP_POTENTIAL_EAM_PHI( p, r, Phi, dPhi );          
          Rho -= rhoCut;   
          Phi -= phiCut;
          double de = ( dRho * ( dEmb + tab.nbh_pt[i][field::dEmb] ) + dPhi ) / r;

          const double drx = tab.drx[i];
          const double dry = tab.dry[i];
          const double drz = tab.drz[i];
          const double fe_x = de * drx;
          const double fe_y = de * dry;
          const double fe_z = de * drz;

          fx  += fe_x;
          fy  += fe_y;
          fz  += fe_z;
          ep  += .5 * Phi;
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

}



