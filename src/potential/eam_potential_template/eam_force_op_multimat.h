/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

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

#define EamPotentialOperatorName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_force)
#define EamPotentialInitName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_init)

#define EamPotentialStr USTAMP_STR(EamPotentialOperatorName)
#define EamParameterInitStr USTAMP_STR(EamPotentialInitName)

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

#ifndef USTAMP_POTENTIAL_EAM_MM_INIT_TYPES
#define USTAMP_POTENTIAL_EAM_MM_INIT_TYPES(p,nt,pe) /**/
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
      //unsigned int cell_a = 0;
      //unsigned int p_a = 0;
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
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a , size_t p_a, exanb::ComputePairParticleContextStart ) const
      {
        ctx.ext.rho = 0.0;
        //ctx.ext.cell_a = cell_a;
        //ctx.ext.p_a = p_a;
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
      {
        static constexpr bool CPAA = NewtonSym &&   gpu_device_execution();
        static constexpr bool LOCK = NewtonSym && ! gpu_device_execution();
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
        static constexpr bool CPAA = NewtonSym &&   gpu_device_execution();
        static constexpr bool LOCK = NewtonSym && ! gpu_device_execution();
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;
        //assert( cell_b < m_nb_cells );
        //assert( p_b < ( m_cell_offsets[cell_b+1] - m_cell_offsets[cell_b] ) );
        const int type_b = cells[cell_b][field::type][p_b];  
        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
        {
          const double r = sqrt( d2 );
          ctx.ext.rho += USTAMP_POTENTIAL_EAM_RHO_NODERIV( p, r, type_b, type_a );

          //size_t p_a_idx = m_cell_particle_offset[ctx.ext.cell_a] + ctx.ext.p_a;
          //size_t p_b_idx = m_cell_particle_offset[cell_b] + p_b;
          //printf("i=%05d , j=%05d , rij=%g,%g,%g\n",int(p_a_idx),int(p_b_idx),cells[cell_a][posfields.e0][p_a],cells[cell_a][posfields.e1][p_a],cells[cell_a][posfields.e2][p_a]);

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
      Mat3d vir = Mat3d{ 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. };
    };

    template<bool NewtonSym, class CPLocksT, class VirFieldT>
    struct SymForceOp
    {
      const onika::cuda::ro_shallow_copy_t<USTAMP_POTENTIAL_PARMS> p = {};
      const uint8_t * __restrict__ m_pair_enabled = nullptr;
      const size_t * __restrict__ m_cell_particle_offset = nullptr;
      const double * __restrict__ m_dEmb_data = nullptr;
      CPLocksT& cp_locks;
      VirFieldT m_vir_field = {};
      
      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart ) const
      {
        ctx.ext.fpi = m_dEmb_data[m_cell_particle_offset[cell_a]+p_a];
        ctx.ext.f = Vec3d { 0. , 0. , 0. };
        if constexpr ( std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> )
        {
          ctx.ext.ep = 0.0;
          ctx.ext.vir = Mat3d{ 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. };
        }
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (const ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
      {
        static constexpr bool CPAA = NewtonSym &&   gpu_device_execution();
        static constexpr bool LOCK = NewtonSym && ! gpu_device_execution();
        if constexpr ( LOCK ) cp_locks[cell_a][p_a].lock();
        if constexpr ( CPAA )
        {
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::fx][p_a] , ctx.ext.f.x );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::fy][p_a] , ctx.ext.f.y );
          ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::fz][p_a] , ctx.ext.f.z );
          if constexpr ( std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> )
          {
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][field::ep][p_a] , ctx.ext.ep );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m11 , ctx.ext.vir.m11 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m12 , ctx.ext.vir.m12 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m13 , ctx.ext.vir.m13 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m21 , ctx.ext.vir.m21 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m22 , ctx.ext.vir.m22 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m23 , ctx.ext.vir.m23 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m31 , ctx.ext.vir.m31 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m32 , ctx.ext.vir.m32 );
            ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_a][m_vir_field][p_a].m33 , ctx.ext.vir.m33 );
          }
        }
        if constexpr ( !CPAA )
        {
          cells[cell_a][field::fx][p_a] += ctx.ext.f.x ;
          cells[cell_a][field::fy][p_a] += ctx.ext.f.y ;
          cells[cell_a][field::fz][p_a] += ctx.ext.f.z ;
          if constexpr ( std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> )
          {
            cells[cell_a][field::ep][p_a] += ctx.ext.ep ;
            cells[cell_a][m_vir_field][p_a] += ctx.ext.vir;
          }
        }
        if constexpr ( LOCK ) cp_locks[cell_a][p_a].unlock();
      }

      template<class ComputeBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (
        ComputeBufferT& ctx,
        const Vec3d& dr,double d2, int type_a, const Mat3d& /*virial*/ ,
        CellParticlesT cells,size_t cell_b, size_t p_b
        , double /*scale*/) const
      {
        static constexpr bool CPAA = NewtonSym &&   gpu_device_execution();
        static constexpr bool LOCK = NewtonSym && ! gpu_device_execution();
        static constexpr bool eflag = std::is_same_v<decltype(ctx.ext),ForceEnergyOpExtStorage> ;
        if( m_pair_enabled!=nullptr && !m_pair_enabled[unique_pair_id(type_a,type_a)] ) return;

        assert( p_b < cells[cell_b].size() );
        const int type_b = cells[cell_b][field::type][p_b];

        if( m_pair_enabled==nullptr || m_pair_enabled[unique_pair_id(type_a,type_b)] )
        {
          const double r = sqrt( d2 );
          assert( r > 0.0 );
          double phi = 0.;
          // dr vaut rj - ri
          Vec3d dr_fe = dr; // en entrée dr_fe vaut rj-ri, en sortie il vaudra (rj-ri) * fpair;
          // TODO: modif USTAMP_POTENTIAL_EAM_MM_FORCE, aka eam_alloy_mm_force pour renvoyer la force seule, et faire la mult à la main en sortie
          USTAMP_POTENTIAL_EAM_MM_FORCE(p,dr_fe,phi,r,ctx.ext.fpi,m_dEmb_data[m_cell_particle_offset[cell_b]+p_b],type_a,type_b);
          ctx.ext.f += dr_fe;
          Mat3d virlocal = tensor( dr_fe, dr ) * -0.5;
          if constexpr ( eflag )
            {
              ctx.ext.ep += .5 * phi;
              ctx.ext.vir += virlocal;
            }
          if constexpr ( NewtonSym )
          {
            if constexpr ( LOCK ) cp_locks[cell_b][p_b].lock();
            if constexpr ( CPAA )
            {            
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fx][p_b] , - dr_fe.x );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fy][p_b] , - dr_fe.y );
              ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::fz][p_b] , - dr_fe.z );
              if constexpr ( eflag )
                {
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][field::ep][p_b] , .5 * phi );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m11 , virlocal.m11 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m12 , virlocal.m12 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m13 , virlocal.m13 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m21 , virlocal.m21 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m22 , virlocal.m22 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m23 , virlocal.m23 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m31 , virlocal.m31 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m32 , virlocal.m32 );
                  ONIKA_CU_BLOCK_ATOMIC_ADD( cells[cell_b][m_vir_field][p_b].m33 , virlocal.m33 );
                }
            }
            if constexpr (!CPAA )
            {
              cells[cell_b][field::fx][p_b] -= dr_fe.x;
              cells[cell_b][field::fy][p_b] -= dr_fe.y;
              cells[cell_b][field::fz][p_b] -= dr_fe.z;
              if constexpr ( eflag )
                {
                  cells[cell_b][field::ep][p_b] += .5 * phi;
                  cells[cell_b][m_vir_field][p_b] += virlocal;
                }
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

  template<bool NewtonSym, class CPLocksT, class VirFieldT> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::SymForceOp<NewtonSym,CPLocksT,VirFieldT> >
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


