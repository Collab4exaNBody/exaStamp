#pragma once

#include <cmath>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>

#include "eam_buffer.h"
#include "potential.h"

#ifndef PRIV_NAMESPACE_NAME
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)
#endif

# define EamPotentialOperatorName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force_multimat)
# define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)

namespace exaStamp
{
  using namespace exanb;

# ifndef USTAMP_POTENTIAL_EAM_PHI_MM
# define EamPotentialPhiMMName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_phi_mm)
# define USTAMP_POTENTIAL_EAM_PHI_MM EamPotentialPhiMMName
  inline void EamPotentialPhiMMName(const USTAMP_POTENTIAL_PARMS & p, double r, double& phiValue, double& dphi, const EAMSpecyPairInfo& /*unused*/ , bool = false /*unused*/ )
  {
    USTAMP_POTENTIAL_EAM_PHI( p, r, phiValue, dphi );
  }
# endif

# ifndef USTAMP_POTENTIAL_EAM_RHO_MM
# define EamPotentialRhoMMName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_rho_mm)
# define USTAMP_POTENTIAL_EAM_RHO_MM EamPotentialRhoMMName
  inline void EamPotentialRhoMMName(const USTAMP_POTENTIAL_PARMS & p, double r, double& rhoValue, double& drho, const EAMSpecyPairInfo& /*unused*/ , bool = false /*unused*/ )
  {
    USTAMP_POTENTIAL_EAM_RHO( p, r, rhoValue, drho );
  }
# endif

  namespace PRIV_NAMESPACE_NAME
  {

    using EamMultiMatParams = EamMultimatParameters<USTAMP_POTENTIAL_PARMS>;

    struct EmbOp
    {
      const EamMultiMatParams* p = nullptr;
      const PhiRhoCutoff* phi_rho_cutoff = nullptr;
      const size_t* m_cell_emb_offset = nullptr;
      double* m_particle_emb = nullptr;
      const uint8_t* m_pair_enabled = nullptr;

      template<class ComputePairBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,

        // particles fields to compute                
        double& ep,
        int type_a,

        // data and locks accessors for neighbors (not used)
        CellParticlesT*
        ) const
      {
        if( ! m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

        double rho = 0.;
	
#       pragma omp simd reduction(+:rho)
        for(size_t i=0;i<n;i++)
        {
          double r = sqrt( tab.d2[i] );
          int type_b = tab.ext.type[i];
          bool pair_inversed = type_a > type_b;
          int pair_id = unique_pair_id( type_a , type_b );
	  
          double rholoc = 0.;
          double drholoc = 0.;
	  if( m_pair_enabled[pair_id] )
	    {
	      USTAMP_POTENTIAL_EAM_RHO_MM( p[unique_pair_id( type_b, type_b )].m_parameters, r, rholoc, drholoc, p[pair_id].m_specy_pair, pair_inversed );
	      //rholoc -= phi_rho_cutoff[unique_pair_id( type_a, type_a )].m_rho_cutoff;
	      rho += rholoc;
	    }
        }

        double emb = 0.;
        double demb = 0.;
	double rhomax = 100.;
        USTAMP_POTENTIAL_EAM_EMB( p[unique_pair_id(type_a,type_a)].m_parameters, rho, emb, demb );
        m_particle_emb[ m_cell_emb_offset[tab.cell] + tab.part ] = demb;	
	if (rho > rhomax) emb += demb * (rho-rhomax);
        ep += emb;
      }

    };


    struct ForceOp
    {
      const EamMultiMatParams* p;
      const PhiRhoCutoff* phi_rho_cutoff = nullptr;
      const size_t* m_cell_emb_offset = nullptr;
      double* m_particle_emb = nullptr;
      const uint8_t* m_pair_enabled = nullptr;

      template<class ComputePairBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBufferT& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        int type_a,
        CellParticlesT* _unused
        ) const
      {
        FakeMat3d virial;
        (*this) ( n,tab,_ep,_fx,_fy,_fz, type_a, virial, _unused );
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
        CellParticlesT*
        ) const
      {
        if( ! m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

	// derivative of Embedding function at atom i
        const double fpi = m_particle_emb[ m_cell_emb_offset[tab.cell] + tab.part ];
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
          double r = tab.d2[i];
          int type_b = tab.ext.type[i];
          bool pair_inversed = ( type_a > type_b );
          int pair_id = unique_pair_id( type_a , type_b );

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
	  if( m_pair_enabled[pair_id] ) {
	    // rhoip = derivative of (density at neighbor atom j due to central atom i)
	    USTAMP_POTENTIAL_EAM_RHO_MM( p[unique_pair_id( type_a, type_a )].m_parameters, r, Rho, rhoip , p[pair_id].m_specy_pair , pair_inversed );
	    // rhojp = derivative of (density at central atom i due to neighbor atom j)
	    USTAMP_POTENTIAL_EAM_RHO_MM( p[unique_pair_id( type_b, type_b )].m_parameters, r, Rho, rhojp , p[pair_id].m_specy_pair , pair_inversed );
	    // phi = pair potential energy
	    // phip = phi'
	    USTAMP_POTENTIAL_EAM_PHI_MM( p[pair_id].m_parameters, r, phi, phip , p[pair_id].m_specy_pair , pair_inversed );
	  }
	  double fpj = tab.ext.emb[i];
	  double psip = fpi * rhojp + fpj * rhoip + phip;	  
	  double fpair = psip/r;
	  
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



