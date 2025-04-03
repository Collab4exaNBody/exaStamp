#pragma once

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>

namespace exaStamp
{

  using namespace exanb;
  
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PaceComputeBuffer
  {
    alignas(onika::memory::DEFAULT_ALIGNMENT) int type[exanb::MAX_PARTICLE_NEIGHBORS];
  };

  struct CopyParticleType
  {
    template<class ComputeBufferT, class FieldArraysT>
    /* ONIKA_HOST_DEVICE_FUNC */ inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, FieldArraysT cells, size_t cell_b, size_t p_b, double weight) const noexcept
    {
      assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
      tab.ext.type[tab.count] = cells[cell_b][field::type][p_b];
      exanb::DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, weight );
    }
  };

  // Force operator
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) PaceForceOp 
  {
    PaceThreadContext* m_thread_ctx = nullptr;
    const size_t n_thread_ctx = 0;
    const bool conv_energy_units = true;
    static constexpr double conv_energy_factor = EXASTAMP_CONST_QUANTITY( 1. * eV );

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells
      ) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3d& virial,
      CellParticlesT cells
      ) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {
      FakeMat3d virial;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class Mat3dT,class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (
      int jnum ,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3dT& virial ,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {
      static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;
      //static constexpr double scale[1][1] = { {1.0} };

      Mat3dT _vir;
      double _en = 0.;
      double _fx = 0.;
      double _fy = 0.;
      double _fz = 0.;
      
      size_t tid = omp_get_thread_num();
      assert( tid < n_thread_ctx );
      PaceThreadContext & pace_ctx = m_thread_ctx[tid];
      auto aceimplptr = pace_ctx.aceimpl;
      aceimplptr->ace->resize_neighbours_cache(jnum);

      double **x = new double*[jnum + 1];
      int *typevec = new int[jnum+1];
      int *jlist = new int[jnum];
      for (int jj = 0; jj < jnum + 1; jj++) {
        x[jj] = new double[3];
      }

      typevec[0] = type+1;
      for (int jj = 0; jj < jnum; jj++) {
        typevec[jj+1] = buf.ext.type[jj]+1;
        jlist[jj] = jj+1;
      }

      x[0][0] = 0.;
      x[0][1] = 0.;
      x[0][2] = 0.;
      
      for (int jj = 1; jj < jnum+1; jj++) {
        x[jj][0] = buf.drx[jj-1];
        x[jj][1] = buf.dry[jj-1];
        x[jj][2] = buf.drz[jj-1];
      }

      // for (int jj = 0; jj < jnum; ++jj) {
      //   int j = jlist[jj];
      //   int type_j = typevec[j];
      //   std::cout << " jj = " << jj << ", j = " << j << ", type_j = " << type_j << std::endl;
      // }

      const double xtmp = x[0][0];
      const double ytmp = x[0][2];
      const double ztmp = x[0][1];
      // std::cout << "xtmp = " << xtmp << "," << ytmp << "," << ztmp << std::endl;
      // std::cout << "xlas = " << x[jnum][0] << "," << x[jnum][1] << "," << x[jnum][2] << std::endl;
      
      // std::cout << "coucou" << std::endl;
      //      aceimplptr->ace->compute_atom(0, x, typevec, jnum, jlist);
      
      // Clean up allocated memory
      for (int jj = 0; jj < jnum + 1; ++jj) {
        delete[] x[jj];
      }
      delete[] x;
      delete[] typevec;
      delete[] jlist;
      
      // double fij[3];
      // for (int jj = 0; jj < jnum; ++jj) {
      //   fij[0]=0.;
      //   fij[1]=0.;
      //   fij[2]=0.;
      //   fij[0] = aceimplptr->ace->neighbours_forces(jj,0);
      //   fij[1] = aceimplptr->ace->neighbours_forces(jj,1);
      //   fij[2] = aceimplptr->ace->neighbours_forces(jj,2);
      //   fij[0] *= conv_energy_factor;
      //   fij[1] *= conv_energy_factor;
      //   fij[2] *= conv_energy_factor;
        
      //   Mat3d v_contrib = tensor( Vec3d{ fij[0] , fij[1] , fij[2] }, Vec3d{ buf.drx[jj], buf.dry[jj], buf.drz[jj] } );
      //   if constexpr ( compute_virial ) { _vir += v_contrib * -1.0; }
        
      //   _fx += fij[0];
      //   _fy += fij[1];
      //   _fz += fij[2];

      //   //        _en += aceimplptr->ace->e_atom * conv_energy_factor;
      //   size_t cell_b=0, p_b=0;
      //   buf.nbh.get(jj, cell_b, p_b);
      //   auto& lock_b = locks[cell_b][p_b];
      //   lock_b.lock();
      //   cells[cell_b][field::fx][p_b] -= fij[0];
      //   cells[cell_b][field::fy][p_b] -= fij[1];
      //   cells[cell_b][field::fz][p_b] -= fij[2];
      //   lock_b.unlock();
      // }
      
    }
  };

}


