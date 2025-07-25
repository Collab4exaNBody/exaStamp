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

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>

namespace exaStamp
{

  using VecPaceThreadContext = std::vector<PaceThreadContext>;
  
  bool hasExtension(const std::string& filename, const std::string& extension) {
    if (filename.length() >= extension.length()) {
      return std::equal(extension.rbegin(), extension.rend(), filename.rbegin());
    }
    return false;
  }

  static char const *const elements_pace[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
  static constexpr int elements_num_pace = sizeof(elements_pace) / sizeof(const char *);
  
  static int AtomicNumberByName_pace(const char *elname)
  {
    for (int i = 1; i < elements_num_pace; i++)
      if (strcmp(elname, elements_pace[i]) == 0) return i;
    return -1;
  }
  
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
    //    VecPaceThreadContext& m_thread_ctx;
    std::vector< std::shared_ptr<ACEImpl>> & m_pace_interfaces;
    const bool conv_energy_units = true;
    static constexpr double conv_energy_factor = EXASTAMP_CONST_QUANTITY( 1. * eV );
    const bool eflag;
    
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

      Mat3dT _vir;
      double _en = 0.;
      double _fx = 0.;
      double _fy = 0.;
      double _fz = 0.;
      
      size_t tid = omp_get_thread_num();

      // ACEImpl* acecalc = m_thread_ctx[tid].aceimpl;
      // acecalc->ace->resize_neighbours_cache(jnum);
      ACEImpl& acecalc = *m_pace_interfaces[tid];
      acecalc.ace->resize_neighbours_cache(jnum);

      // compute_atom needs lammps-like arrays
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
      
      for (int jj = 0; jj < jnum; jj++) {
        x[jj+1][0] = buf.drx[jj];
        x[jj+1][1] = buf.dry[jj];
        x[jj+1][2] = buf.drz[jj];
      }
      
      //      acecalc->ace->compute_atom(0, x, typevec, jnum, jlist);
      acecalc.ace->compute_atom(0, x, typevec, jnum, jlist);      
      
      for (int jj = 0; jj < jnum + 1; ++jj) {
        delete[] x[jj];
      }
      delete[] x;
      delete[] typevec;
      delete[] jlist;
      
      double fij[3];
      fij[0]=0.;
      fij[1]=0.;
      fij[2]=0.;      
      for (int jj = 0; jj < jnum; ++jj) {
        // fij[0]=acecalc->ace->neighbours_forces(jj,0) * conv_energy_factor;
        // fij[1]=acecalc->ace->neighbours_forces(jj,1) * conv_energy_factor;
        // fij[2]=acecalc->ace->neighbours_forces(jj,2) * conv_energy_factor;
        fij[0]=acecalc.ace->neighbours_forces(jj,0) * conv_energy_factor;
        fij[1]=acecalc.ace->neighbours_forces(jj,1) * conv_energy_factor;
        fij[2]=acecalc.ace->neighbours_forces(jj,2) * conv_energy_factor;

        Mat3d v_contrib = tensor( Vec3d{ fij[0] , fij[1] , fij[2] }, Vec3d{ buf.drx[jj], buf.dry[jj], buf.drz[jj] } );
        if constexpr ( compute_virial ) { _vir += v_contrib * -1.0; }
        
        _fx += fij[0];
        _fy += fij[1];
        _fz += fij[2];
        
        size_t cell_b=0, p_b=0;
        buf.nbh.get(jj, cell_b, p_b);
        auto& lock_b = locks[cell_b][p_b];
        lock_b.lock();
        cells[cell_b][field::fx][p_b] -= fij[0];
        cells[cell_b][field::fy][p_b] -= fij[1];
        cells[cell_b][field::fz][p_b] -= fij[2];
        lock_b.unlock();
      }

      lock_a.lock();
      if (eflag) {      
        //        en += acecalc->ace->e_atom * conv_energy_factor;
        en += acecalc.ace->e_atom * conv_energy_factor;
      }
      fx += _fx;
      fy += _fy;
      fz += _fz;
      if constexpr ( compute_virial ) { virial += _vir; }
      lock_a.unlock();
      
    }
  };

}


