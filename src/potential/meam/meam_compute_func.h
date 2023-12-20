#pragma once

#include "meam_parameters.h"
#include "meam.h"
#include <exanb/compute/compute_pair_buffer.h>

//#include <exanb/debug/print_particle.h>

namespace exaStamp
{
  using onika::memory::DEFAULT_ALIGNMENT ;

  // MEAM Extra storage needed along with compute buffer
  template<size_t NMaxParticles = MeamParameters::MAX_PARTICLE_NEIGHBORS >
  struct alignas(DEFAULT_ALIGNMENT) MeamExtraStorageT
  {
    // temporary arrays
    alignas(DEFAULT_ALIGNMENT) double S[NMaxParticles];
//    alignas(DEFAULT_ALIGNMENT) double dfx[NMaxParticles];
//    alignas(DEFAULT_ALIGNMENT) double dfy[NMaxParticles];
//    alignas(DEFAULT_ALIGNMENT) double dfz[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double dfSx[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double dfSy[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double dfSz[NMaxParticles];
//    alignas(DEFAULT_ALIGNMENT) double den[NMaxParticles];

    alignas(DEFAULT_ALIGNMENT) double drho0x[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho0y[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho0z[NMaxParticles];
    
    alignas(DEFAULT_ALIGNMENT) double drho1x[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho1y[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho1z[NMaxParticles];

    alignas(DEFAULT_ALIGNMENT) double drho2x[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho2y[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho2z[NMaxParticles];

    alignas(DEFAULT_ALIGNMENT) double drho3x[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho3y[NMaxParticles];
    alignas(DEFAULT_ALIGNMENT) double drho3z[NMaxParticles];
  };
  using MeamExtraStorage = MeamExtraStorageT<>;

  // MEAM Compute functor
  template<class CellParticles, class ComputeBufferT, bool CPAA = true> // CPAA = Central particle's field requires CU atomic add intrinsic
  struct MeamForceComputeFunctor
  {
    // core computations and parameters
    MeamPotential m_pot;

    ONIKA_HOST_DEVICE_FUNC inline void operator ()
      (
      size_t n,
      ComputeBufferT& tab,
      double& _ep,
      double& _fx,
      double& _fy,
      double& _fz,
      Mat3d& , // UNUSED virial
      CellParticles* cells
      ) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock particle_lock;
      (*this) ( n,tab,_ep,_fx,_fy,_fz,cells,locks,particle_lock );
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator ()
      (
      size_t n,
      ComputeBufferT& tab,
      double& _ep,
      double& _fx,
      double& _fy,
      double& _fz,
      CellParticles* cells
      ) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock particle_lock;
      (*this) ( n,tab,_ep,_fx,_fy,_fz,cells,locks,particle_lock );
    }

    template<class GridCellParticleLocksT , class ParticleLockT >
    ONIKA_HOST_DEVICE_FUNC inline void operator ()
      (
      size_t n,
      ComputeBufferT& tab,
      double& _ep,
      double& _fx,
      double& _fy,
      double& _fz,
      Mat3d& , // UNUSED virial
      CellParticles* cells,
      GridCellParticleLocksT locks,
      ParticleLockT & particle_lock
      ) const
    {
      (*this) ( n,tab,_ep,_fx,_fy,_fz,cells,locks,particle_lock );
    }

    template<class GridCellParticleLocksT , class ParticleLockT >
    ONIKA_HOST_DEVICE_FUNC inline void operator ()
      (
      size_t n,
      ComputeBufferT& tab,
      double& _ep,
      double& _fx,
      double& _fy,
      double& _fz,
      CellParticles* cells,
      GridCellParticleLocksT locks,
      ParticleLockT & particle_lock
      ) const
    {
      using onika::cuda::max;

      assert( /* n >= 0 && */ n <= MeamParameters::MAX_PARTICLE_NEIGHBORS );
      
      // compute screening term
      m_pot.screeningFunction( tab.drx, tab.dry, tab.drz, tab.ext.S, n );


      // remove screened neighbors  
      for(size_t i=0; i<n; )
      {
        if( tab.ext.S[i] == 0. )
        {
          tab.copy( n-1 , i );
          tab.ext.S[i] = tab.ext.S[n-1];
          -- n;
        }
        else { ++ i; }
      }

      double Fe = 0.;
      double Fx = 0.;          
      double Fy = 0.;          
      double Fz = 0.;

      {
        // use those arrays as scratch arrays for phi computation internals.
        // it's legit because they're only needed and used afterward, in derivedRhoS
        double * tmp_dfx = tab.ext.drho0x;
        double * tmp_dfy = tab.ext.drho0y;
        double * tmp_dfz = tab.ext.drho0z;
        double * tmp_den = tab.ext.drho1x;

        // compute phi
        m_pot.phi( tmp_dfx, tmp_dfy, tmp_dfz, tab.ext.dfSx, tab.ext.dfSy, tab.ext.dfSz, tmp_den, tab.drx, tab.dry, tab.drz, tab.ext.S, n );

        // intermediate step, symetric contributions are added to the local and neighbor atoms
        // ---- writeForceMEAM begin ----
#       pragma omp simd reduction(+:Fx,Fy,Fz,Fe)
        for(size_t i=0;i<n;i++)
        {
          Fx += tmp_dfx[i];
          Fy += tmp_dfy[i];
          Fz += tmp_dfz[i];
          Fe += tmp_den[i];
          // TODO: virial
        }
      }

#     ifndef XSTAMP_MEAM_MERGE_FORCE_UPDATES
      for(size_t i=0;i<n;i++)
      {
        size_t cell_j=0, p_j=0;
        tab.nbh.get( i, cell_j, p_j );
        auto & cell = cells[cell_j];
        locks[cell_j][p_j].lock();
        ONIKA_CU_BLOCK_ATOMIC_ADD( cell[field::fx][p_j] , -tab.ext.dfSx[i] );
        ONIKA_CU_BLOCK_ATOMIC_ADD( cell[field::fy][p_j] , -tab.ext.dfSy[i] );
        ONIKA_CU_BLOCK_ATOMIC_ADD( cell[field::fz][p_j] , -tab.ext.dfSz[i] );
        locks[cell_j][p_j].unlock();
        tab.ext.dfSx[i] = 0.;
        tab.ext.dfSy[i] = 0.;
        tab.ext.dfSz[i] = 0.;
      }

      particle_lock.lock();
      if constexpr ( CPAA )
      {
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fx, Fx ); 
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fy, Fy );
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fz, Fz );
      }
      if constexpr ( !CPAA )
      {
        _fx += Fx;
        _fy += Fy;
        _fz += Fz;
      }
      particle_lock.unlock();
      _ep += Fe;
      Fx=0.; Fy=0.; Fz=0.; Fe=0.;
#     endif

      // -- Electronic term --
  	  // compute rho0 (scalar)
  	  double rho0 = m_pot.compute_rho_0( tab.ext.S, tab.drx, tab.dry, tab.drz, n );
    	//to avoid issue in case rho0=0 => rho=0 and log(rho) * rho != 0
    	rho0 = max(rho0,1.e-30);

    	//compute first part of rho0,rho 2, rho 3 and rho1 (pour l'instant)
      Vec3d rhod0 = { 0., 0., 0. };
      Vec3d rhod1 = { 0., 0., 0. };
      Vec3d rhod2 = { 0., 0., 0. };
      Vec3d rhod3 = { 0., 0., 0. };
      double rho1 = 0.;
      double rho2 = 0.;
      double rho3 = 0.;
    	m_pot.derivedRhoS(tab.ext.S,
    	                  tab.drx, tab.dry, tab.drz,
                              tab.ext.drho0x, tab.ext.drho0y, tab.ext.drho0z, rhod0,
                        rho1, tab.ext.drho1x, tab.ext.drho1y, tab.ext.drho1z, rhod1,
                        rho2, tab.ext.drho2x, tab.ext.drho2y, tab.ext.drho2z, rhod2,
                        rho3, tab.ext.drho3x, tab.ext.drho3y, tab.ext.drho3z, rhod3,
                        n );

      // ==== final step (writeRho_0123_MEAM in ExaStamp v1)  =========
      // Apply the potential to get the rho terms
      double rho = 0.0;
      m_pot.rho(rho, rho0,rho1,rho2,rho3);

      // compute F(rho/Z)
      double fEmbed = 0.0, df = 0.0;
      m_pot.fEmbed(rho, fEmbed, df);

      double tmpRho0 = 0.0, tmpRho1 = 0.0, tmpRho2 = 0.0, tmpRho3 = 0.0, tmpG = 0.0, tmpdG = 0.0;
      m_pot.getCoeffDerivateRho(rho0,rho1,rho2,rho3, tmpRho0, tmpRho1, tmpRho2, tmpRho3, tmpG, tmpdG );
      const double Coeff = ( df /*/ m_mass */ ) / m_pot.p.Z; // Z=12

      // update force if not in a ghost cell
      //if (ghostLayer[cell]==0)
      //{
      // Update derivate of rho 0,1,2,3 on the particle
      Fx += Coeff*(tmpG*rhod0.x+rho0*tmpdG*(tmpRho0*rhod0.x+tmpRho1*rhod1.x+tmpRho2*rhod2.x+tmpRho3*rhod3.x)) ;
      Fy += Coeff*(tmpG*rhod0.y+rho0*tmpdG*(tmpRho0*rhod0.y+tmpRho1*rhod1.y+tmpRho2*rhod2.y+tmpRho3*rhod3.y)) ;
      Fz += Coeff*(tmpG*rhod0.z+rho0*tmpdG*(tmpRho0*rhod0.z+tmpRho1*rhod1.z+tmpRho2*rhod2.z+tmpRho3*rhod3.z)) ;
      Fe += fEmbed;

/*      if( ids!=nullptr && ids[tab.part]==332350 )
      {
        lout <<ids[tab.part]<<" : N="<<n<<" Fx="<<Fx<< " Fy="<<Fy<< " Fz="<<Fz<<" Ep="<<Fe<< std::endl;
      }*/

      particle_lock.lock();
      if constexpr ( CPAA )
      {
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fx, Fx ); 
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fy, Fy );
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fz, Fz );
      }
      if constexpr ( !CPAA )
      {
        _fx += Fx;
        _fy += Fy;
        _fz += Fz;
      }
      particle_lock.unlock();
      _ep += Fe; // not concurently wirtten by neighbor atoms (local computation)
      //}

      // add contribution on neighbors
      for(size_t i=0;i<n;i++)
      {
        size_t cell_j=0, p_j=0;
        tab.nbh.get( i, cell_j, p_j );
        auto & cell = cells[cell_j];
        const double nbh_fx_contrib = tab.ext.dfSx[i] +
          Coeff*(rho0*tmpdG*(tmpRho0*tab.ext.drho0x[i]+tmpRho1*tab.ext.drho1x[i]+tmpRho2*tab.ext.drho2x[i]+tmpRho3*tab.ext.drho3x[i])+tmpG*tab.ext.drho0x[i]) ;
        const double nbh_fy_contrib = tab.ext.dfSy[i] +
          Coeff*(rho0*tmpdG*(tmpRho0*tab.ext.drho0y[i]+tmpRho1*tab.ext.drho1y[i]+tmpRho2*tab.ext.drho2y[i]+tmpRho3*tab.ext.drho3y[i])+tmpG*tab.ext.drho0y[i]) ;
        const double nbh_fz_contrib = tab.ext.dfSz[i] +
          Coeff*(rho0*tmpdG*(tmpRho0*tab.ext.drho0z[i]+tmpRho1*tab.ext.drho1z[i]+tmpRho2*tab.ext.drho2z[i]+tmpRho3*tab.ext.drho3z[i])+tmpG*tab.ext.drho0z[i]) ;
        locks[cell_j][p_j].lock();
        ONIKA_CU_BLOCK_ATOMIC_ADD( cell[field::fx][p_j] , -nbh_fx_contrib );
        ONIKA_CU_BLOCK_ATOMIC_ADD( cell[field::fy][p_j] , -nbh_fy_contrib );
        ONIKA_CU_BLOCK_ATOMIC_ADD( cell[field::fz][p_j] , -nbh_fz_contrib );
        locks[cell_j][p_j].unlock();
      }
    }

    ONIKA_HOST_DEVICE_FUNC inline double rcut() const { return m_pot.p.Rcut; }
  };

}

namespace exanb
{

  // specialize functor traits to allow Cuda execution space
  template<class CellParticles, class ComputeBufferT>
  struct ComputePairTraits< exaStamp::MeamForceComputeFunctor<CellParticles,ComputeBufferT> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = true;
  };

}

