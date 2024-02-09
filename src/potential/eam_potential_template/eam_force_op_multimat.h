#pragma once

#include <cmath>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>

#include <exaStamp/potential/eam/eam_buffer.h>
#include "potential.h"

#ifndef PRIV_NAMESPACE_NAME
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)
#endif

# define EamPotentialOperatorName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force)
# define EamPotentialEmbOnlyName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_emb)
# define EamPotentialForceOnlyName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force_reuse_emb)

# define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)
# define EamPotentialEmbOnlyStr USTAMP_STR(EamPotentialEmbOnlyName)
# define EamPotentialForceOnlyStr USTAMP_STR(EamPotentialForceOnlyName)

namespace exaStamp
{
  using namespace exanb;

  namespace PRIV_NAMESPACE_NAME
  {

    using EamMultiMatParams = EamMultimatParameters< USTAMP_POTENTIAL_PARMS >;
    using EamMultiMatParamsReadOnly = EamMultimatParametersRO< onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> >;

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
        if( ! m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

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
          if( m_pair_enabled[pair_id] )
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
        double fpi,
        CellParticlesT
        ) const
      {
        if( ! m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

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
          [[maybe_unused]] const bool pair_inversed = ( type_a > type_b );
          const int pair_id = unique_pair_id( type_a , type_b );

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
	        if( m_pair_enabled[pair_id] )
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



