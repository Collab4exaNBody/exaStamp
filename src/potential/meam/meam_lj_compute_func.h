#pragma once

#include "meam_compute_func.h"
#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>
#include <exanb/core/particle_type_pair.h>

#ifndef XSTAMP_MEAM_MULTIMAT_MAX_TYPES
#define XSTAMP_MEAM_MULTIMAT_MAX_TYPES 4
#endif

namespace exaStamp
{
  struct MeamLJParms
  {
    LennardJonesParms lj;
    double rcut2 = 0.0;
    double ecut = 0.0;
  };

  struct UserMeamLJPair
  {
    int pair_id = -1;
    double rcut;
    MeamLJParms lj;
    std::string type_a;
    std::string type_b;
  };

  using UserMeamLJParameters = std::vector<UserMeamLJPair>;

  struct MeamLJMultiMatParms
  {
    static constexpr unsigned int MAX_TYPES = XSTAMP_MEAM_MULTIMAT_MAX_TYPES;
    static constexpr unsigned int MAX_TYPE_PAIRS = exanb::UnorderedPairCount<MAX_TYPES>::count;
    int m_meam_type = -1;
    int m_meam_pair_id = -1;
    double m_meam_rcut2 = 0.0;
    MeamLJParms m_lj_parms[MAX_TYPE_PAIRS];
  };

  // MEAM Extra storage needed along with compute buffer
  template<size_t NMaxParticles = MeamParameters::MAX_PARTICLE_NEIGHBORS >
  struct alignas(DEFAULT_ALIGNMENT) MeamLJExtraStorageT
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
    
    const MeamLJMultiMatParms * m_meam_lj_multimat = nullptr;
    size_t m_meam_lj_nb_lj_pairs = 0;
  };
  using MeamLJExtraStorage = MeamLJExtraStorageT<>;

  struct MeamLJNeighborFilter
  {
    template<typename ComputeBufferT, typename FieldArraysT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, double weight) const noexcept
    {      
      const MeamLJMultiMatParms* mmp = tab.ext.m_meam_lj_multimat;
      const unsigned int type_a = tab.ta;
      const unsigned int type_b = cells[cell_b][field::type][p_b];
      const unsigned int pair_id = unique_pair_id(type_a,type_b);
      
      if( int(pair_id) != mmp->m_meam_pair_id )
      {
        if( d2 <= mmp->m_lj_parms[pair_id].rcut2 )
        {
          const double r = sqrt(d2);
          //const double weight = tab.nbh_data.get(i);
          double e=0.0, de=0.0;
          lj_compute_energy( mmp->m_lj_parms[pair_id].lj , PairPotentialMinimalParameters{} , r , e , de );
          e -= mmp->m_lj_parms[pair_id].ecut;
          de /= r;
          e *= weight;
          de *= weight;        
          const Vec3d fe = dr * de;
          tab.ext.S[0] += fe.x;
          tab.ext.S[1] += fe.y;
          tab.ext.S[2] += fe.z;
          tab.ext.S[3] += .5 * e;
          ++ tab.ext.m_meam_lj_nb_lj_pairs;
        }
      }
      else if ( d2 <= mmp->m_meam_rcut2 )
      {
#       ifdef XSTAMP_MEAM_ENFORCE_OVERFLOW_CHECK
        const auto overflow_check = std::true_type {};
#       else
        const auto overflow_check = std::false_type {};
#       endif
        exanb::DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, weight, overflow_check );
      }
    }
  };

  template<class CellParticles, class ComputeBufferT>
  struct MeamLJForceComputeFunctor
  {
#   ifdef XSTAMP_MEAM_MERGE_FORCE_UPDATES
    using BaseMeamFunctor = MeamForceComputeFunctor<CellParticles,ComputeBufferT,false>;
#   else
    using BaseMeamFunctor = MeamForceComputeFunctor<CellParticles,ComputeBufferT,true>;
#   endif
    BaseMeamFunctor m_nofilter_func;
    MeamLJMultiMatParms m_meam_lj_multimat;

    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, CellParticles* cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart) const
    {
      tab.ext.m_meam_lj_multimat = & m_meam_lj_multimat;
      tab.ext.m_meam_lj_nb_lj_pairs = 0;
      tab.ta = cells[cell_a][field::type][p_a];
      tab.ext.S[0] = 0.0;
      tab.ext.S[1] = 0.0;
      tab.ext.S[2] = 0.0;
      tab.ext.S[3] = 0.0;
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () (int n, ComputeBufferT& tab, double& _ep, double& _fx, double& _fy, double& _fz, Mat3d&, CellParticles* cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock particle_lock;
      (*this) ( n,tab,_ep,_fx,_fy,_fz,cells,locks,particle_lock );
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator () (int n, ComputeBufferT& tab,double& _ep,double& _fx,double& _fy,double& _fz,CellParticles* cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock particle_lock;
      (*this) ( n,tab,_ep,_fx,_fy,_fz,cells,locks,particle_lock );
    }

    template<class GridCellParticleLocksT , class ParticleLockT >
    ONIKA_HOST_DEVICE_FUNC inline void operator () (int n, ComputeBufferT& tab,double& _ep,double& _fx,double& _fy,double& _fz,Mat3d&,CellParticles* cells,GridCellParticleLocksT locks,ParticleLockT & particle_lock) const
    {
      (*this) ( n,tab,_ep,_fx,_fy,_fz,cells,locks,particle_lock );
    }

    template<class GridCellParticleLocksT , class ParticleLockT >
    ONIKA_HOST_DEVICE_FUNC inline void operator ()
      (
      int n,
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
      // get central particle's type
      //const unsigned int type_a = tab.ta; // static_cast<unsigned int>( tab.ext.S[0] );
      assert( n==0 || tab.ta==m_meam_lj_multimat.m_meam_type );

      // get contributions from LJ pairs stored in S array

      // this one doesn't need to be guarded
#     ifdef XSTAMP_MEAM_MERGE_FORCE_UPDATES
      double fx = 0.0;
      double fy = 0.0;
      double fz = 0.0;
      double ep = 0.0; 
      if( tab.ext.m_meam_lj_nb_lj_pairs > 0 )
      {
        fx = tab.ext.S[0];
        fy = tab.ext.S[1];
        fz = tab.ext.S[2];
        ep = tab.ext.S[3];
      }
#     else
      if( tab.ext.m_meam_lj_nb_lj_pairs > 0 )
      {
        particle_lock.lock();
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fx , tab.ext.S[0] );
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fy , tab.ext.S[1] );
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fz , tab.ext.S[2] );
        particle_lock.unlock();
        _ep += tab.ext.S[3];
      }
#     endif

      if( n > 0 )
      {
#       ifdef XSTAMP_MEAM_MERGE_FORCE_UPDATES
        FakeParticleLock no_particle_lock;
        m_nofilter_func( n, tab, ep, fx, fy, fz, cells, locks, no_particle_lock );
#       else
        m_nofilter_func( n, tab, _ep, _fx, _fy, _fz, cells, locks, particle_lock );
#       endif
      }
    
#     ifdef XSTAMP_MEAM_MERGE_FORCE_UPDATES
      if( n > 0 || tab.ext.m_meam_lj_nb_lj_pairs > 0 )
      {
        particle_lock.lock();
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fx , fx );
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fy , fy );
        ONIKA_CU_BLOCK_ATOMIC_ADD( _fz , fz );
        particle_lock.unlock();
        _ep += ep;
      }
#     endif

      tab.ext.m_meam_lj_nb_lj_pairs = 0;
    }
    
    ONIKA_HOST_DEVICE_FUNC inline double rcut() const { return m_nofilter_func.m_pot.p.Rcut; }
  };

}

namespace exanb
{

  template<class CellParticles, class ComputeBufferT>
  struct ComputePairTraits< exaStamp::MeamLJForceComputeFunctor<CellParticles,ComputeBufferT> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = true;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<class CellParticles, class ComputeBufferT>
  struct ComputePairParticleContextTraits< exaStamp::MeamLJForceComputeFunctor<CellParticles,ComputeBufferT> >
  {
    static inline constexpr bool HasParticleContextStart = true;    
    static inline constexpr bool HasParticleContextStop = false;    
  };

}

#include <yaml-cpp/yaml.h>

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::UserMeamLJPair>
  {
    static bool decode(const Node& node, exaStamp::UserMeamLJPair& v)
    {
      v = exaStamp::UserMeamLJPair{};
      if( !node.IsMap() ) { return false; }
      v.type_a = node["type_a"].as<std::string>();
      v.type_b = node["type_b"].as<std::string>();
      v.lj.lj.epsilon = node["epsilon"].as<exanb::Quantity>().convert();
      v.lj.lj.sigma = node["sigma"].as<exanb::Quantity>().convert();
      v.rcut = node["rcut"].as<exanb::Quantity>().convert();
      double e=0.0, de=0.0;
      if( v.rcut > 0.0 ) exaStamp::lj_compute_energy( v.lj.lj , exaStamp::PairPotentialMinimalParameters{} , v.rcut , e , de );    
      v.lj.ecut = e;
      v.lj.rcut2 = v.rcut * v.rcut;
      return true;
    }
  };
}


