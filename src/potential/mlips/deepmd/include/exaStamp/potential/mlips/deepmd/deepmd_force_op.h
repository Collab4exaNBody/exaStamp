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
#include <exanb/core/concurent_add_contributions.h>
#include <onika/cuda/cuda.h>

using namespace exanb;
using namespace exaStamp;

template<class TypeFieldT>
struct alignas(onika::memory::DEFAULT_ALIGNMENT) DPMDForceOp
{

  DPMDThreadContext* m_thread_ctx = nullptr;

  TypeFieldT m_type_field = {};
  
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
   size_t n,
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

    static constexpr bool CPAA = true &&   gpu_device_execution();
    static constexpr bool LOCK = true && ! gpu_device_execution();
    
    Mat3dT _vir;
    double _en = 0.;
    double _fx = 0.;
    double _fy = 0.;
    double _fz = 0.;
    
    static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;
    size_t tid = omp_get_thread_num();
    DPMDThreadContext & dpmd_ctx = m_thread_ctx[tid];
    auto dpmd_ptr = dpmd_ctx.dpmd_model;

    // ================================================================
    // The buffer contains ALL neighbors within 2*rcut of atom 0.
    // n = total number of neighbors in the buffer.
    //
    // DeePMD array layout:
    //   index 0         : central atom (placed at origin)
    //   index 1 .. n    : neighbors from buffer (all treated as ghost)
    // ================================================================
    const int nn = static_cast<int>(n);
    const int np = nn + 1;  // total atoms in cluster

    std::vector<double> dcoord( 3 * np, 0.0 );
    std::vector<int>    atype( np );
    std::vector<double> cell = {};  // no PBC for local cluster

    // Central atom at origin
    atype[0] = type;
    for (int i = 0; i < nn; ++i) {
      dcoord[3 * (i + 1) + 0] = buf.drx[i];
      dcoord[3 * (i + 1) + 1] = buf.dry[i];
      dcoord[3 * (i + 1) + 2] = buf.drz[i];
      //      atype[i + 1] = buf.ext.type[i];
      size_t cell_b = 0, p_b = 0;
      buf.nbh.get(i, cell_b, p_b);
      atype[i + 1] = cells[cell_b][m_type_field][p_b];
    }

    // ================================================================
    // Build FULL SYMMETRIC neighbor list.
    //
    // Atom 0 (local):
    //   lists all atoms j+1 where dist(0, j) < rcut
    //
    // Ghost atom j+1:
    //   lists atom 0   if dist(0, j) < rcut
    //   lists atom k+1 if dist(j, k) < rcut
    //
    // The 2*rcut gather sphere ensures that all atoms needed to
    // complete the environment of atom 0's rcut-neighbors are present.
    // ================================================================
    std::vector<std::vector<int>> nlist_vec(np);

    // Atom 0's neighbors: atoms within rcut of center
    for (int i = 0; i < nn; ++i) nlist_vec[0].push_back(i + 1);

    // Convert to LAMMPS-style InputNlist
    std::vector<int>  ilist(np), numneigh(np);
    std::vector<int*> firstneigh(np);
    deepmd::InputNlist nlist(np, &ilist[0], &numneigh[0], &firstneigh[0]);
    deepmd::convert_nlist(nlist, nlist_vec);

    // ================================================================
    // Call DeePMD compute.
    //
    // nloc   = 1   (only atom 0 is local)
    // nghost = nn  (atoms 1..np-1 are ghost)
    //
    // Result:
    //   e    = E_0       (atomic energy of central atom ONLY)
    //   f[k] = -dE_0/dr_k  for all atoms k in the cluster
    // ================================================================
    double e;
    std::vector<double> f, v;
    const int nghost = nn;
    const int ago    = 0;
    dpmd_ptr->compute(e, f, v, dcoord, atype, cell, nghost, nlist, ago);

    // ================================================================
    // Unit conversion
    //   DeePMD: energy in eV, forces in eV/Angstrom
    // ================================================================
    const double conv_e_factor = ONIKA_QUANTITY( 1. * eV ).convert();
    const double conv_f_factor = ONIKA_QUANTITY( 1. * eV/ang ).convert();

    // ================================================================
    // Energy: E_0 (atomic energy of the central atom only)
    // ================================================================
    _en = e * conv_e_factor;

    // ================================================================
    // Force on central atom (index 0): -dE_0/dr_0
    // ================================================================
    _fx = f[0] * conv_f_factor;
    _fy = f[1] * conv_f_factor;
    _fz = f[2] * conv_f_factor;

    // ================================================================
    // Distribute partial forces to neighbors WITHIN rcut.
    //
    // f[3*(j+1) .. 3*(j+1)+2] = -dE_0 / dr_j
    //
    // Only neighbors within rcut get a non-negligible contribution:
    // E_0's descriptor smoothly vanishes at rcut, so atoms beyond
    // rcut have (essentially) zero force from E_0.
    // Those outer-shell atoms (rcut < r < 2*rcut) exist only to
    // provide complete environments for the inner neighbors.
    //
    // When every atom in the system has been processed as a center,
    // the total force on any atom j is:
    //   F_j = sum_over_all_centers_i { -dE_i / dr_j }
    //       = -dE_total / dr_j     (correct!)
    // ================================================================
    for (int j = 0; j < nn; ++j)
    {
      size_t cell_b = 0, p_b = 0;
      buf.nbh.get(j, cell_b, p_b);
      // concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double>
      //   (
      //    locks[cell_b][p_b],
      //    cells[cell_b][field::fx][p_b],
      //    cells[cell_b][field::fy][p_b],
      //    cells[cell_b][field::fz][p_b],
      //    f[3 * (j + 1) + 0] * conv_f_factor,
      //    f[3 * (j + 1) + 1] * conv_f_factor,
      //    f[3 * (j + 1) + 2] * conv_f_factor
      //    );
      ONIKA_CU_ATOMIC_ADD( cells[cell_b][field::fx][p_b], f[3 * (j + 1) + 0] * conv_f_factor );
      ONIKA_CU_ATOMIC_ADD( cells[cell_b][field::fy][p_b], f[3 * (j + 1) + 1] * conv_f_factor );
      ONIKA_CU_ATOMIC_ADD( cells[cell_b][field::fz][p_b], f[3 * (j + 1) + 2] * conv_f_factor );      

    }

    // ================================================================
    // Virial contribution from E_0 (if requested)
    // ================================================================
    // if constexpr (compute_virial) {
    //   _vir.xx = v[0] * conv_e_factor;
    //   _vir.xy = v[1] * conv_e_factor;
    //   _vir.xz = v[2] * conv_e_factor;
    //   _vir.yx = v[3] * conv_e_factor;
    //   _vir.yy = v[4] * conv_e_factor;
    //   _vir.yz = v[5] * conv_e_factor;
    //   _vir.zx = v[6] * conv_e_factor;
    //   _vir.zy = v[7] * conv_e_factor;
    //   _vir.zz = v[8] * conv_e_factor;
    // }

    // ================================================================
    // Accumulate central atom's energy and force
    // ================================================================
    // concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double>
    //   (
    //    lock_a,
    //    fx,
    //    fy,
    //    fz,
    //    en,
    //    _fx,
    //    _fy,
    //    _fz,
    //    _en
    //    );
    ONIKA_CU_ATOMIC_ADD( fx, _fx );
    ONIKA_CU_ATOMIC_ADD( fy, _fy );
    ONIKA_CU_ATOMIC_ADD( fz, _fz );
    ONIKA_CU_ATOMIC_ADD( en, _en );
    
  }
};

template<class TypeFieldT>
struct alignas(onika::memory::DEFAULT_ALIGNMENT) DPMDForceOpBis
{

  DPMDThreadContextPT* m_thread_ctx = nullptr;

  TypeFieldT m_type_field = {};
  
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
   size_t n,
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

    static constexpr bool CPAA = true &&   gpu_device_execution();
    static constexpr bool LOCK = true && ! gpu_device_execution();
    
    Mat3dT _vir;
    double _en = 0.;
    double _fx = 0.;
    double _fy = 0.;
    double _fz = 0.;
    
    static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;
    size_t tid = omp_get_thread_num();
    DPMDThreadContextPT & dpmd_ctx = m_thread_ctx[tid];
    auto dpmd_ptr = dpmd_ctx.dpmd_model;

    // ================================================================
    // The buffer contains ALL neighbors within 2*rcut of atom 0.
    // n = total number of neighbors in the buffer.
    //
    // DeePMD array layout:
    //   index 0         : central atom (placed at origin)
    //   index 1 .. n    : neighbors from buffer (all treated as ghost)
    // ================================================================
    const int nn = static_cast<int>(n);
    const int np = nn + 1;  // total atoms in cluster

    std::vector<double> dcoord( 3 * np, 0.0 );
    std::vector<int>    atype( np );
    const std::vector<double> cell = {};  // no PBC for local cluster
    const std::vector<double> fparam = {};
    const std::vector<double> aparam = {};
    
    // Central atom at origin
    atype[0] = type;
    for (int i = 0; i < nn; ++i) {
      dcoord[3 * (i + 1) + 0] = buf.drx[i];
      dcoord[3 * (i + 1) + 1] = buf.dry[i];
      dcoord[3 * (i + 1) + 2] = buf.drz[i];
      //      atype[i + 1] = buf.ext.type[i];
      size_t cell_b = 0, p_b = 0;
      buf.nbh.get(i, cell_b, p_b);
      atype[i + 1] = cells[cell_b][m_type_field][p_b];
    }

    // ================================================================
    // Build FULL SYMMETRIC neighbor list.
    //
    // Atom 0 (local):
    //   lists all atoms j+1 where dist(0, j) < rcut
    //
    // Ghost atom j+1:
    //   lists atom 0   if dist(0, j) < rcut
    //   lists atom k+1 if dist(j, k) < rcut
    //
    // The 2*rcut gather sphere ensures that all atoms needed to
    // complete the environment of atom 0's rcut-neighbors are present.
    // ================================================================
    std::vector<std::vector<int>> nlist_vec(np);

    // Atom 0's neighbors: atoms within rcut of center
    for (int i = 0; i < nn; ++i) nlist_vec[0].push_back(i + 1);

    // Convert to LAMMPS-style InputNlist
    std::vector<int>  ilist(np), numneigh(np);
    std::vector<int*> firstneigh(np);
    deepmd::InputNlist nlist(np, &ilist[0], &numneigh[0], &firstneigh[0]);
    deepmd::convert_nlist(nlist, nlist_vec);

    // ================================================================
    // Call DeePMD compute.
    //
    // nloc   = 1   (only atom 0 is local)
    // nghost = nn  (atoms 1..np-1 are ghost)
    //
    // Result:
    //   e    = E_0       (atomic energy of central atom ONLY)
    //   f[k] = -dE_0/dr_k  for all atoms k in the cluster
    // ================================================================
    std::vector<double> e, f, v, atom_e, atom_v;
    const int nghost = nn;
    const int ago    = 0;
    dpmd_ptr->computew(e,
                       f,
                       v,
                       atom_e,
                       atom_v,
                       dcoord,
                       atype,
                       cell,
                       nghost,
                       nlist,
                       ago,
                       fparam,
                       aparam,
                       false);

    // int gpu_num = torch::cuda::device_count();
    // bool gpu_id = (gpu_num > 0) ? (gpu_rank % gpu_num) : 0;
    // bool gpu_enabled = torch::cuda::is_available();
    
    // lout << "GPU enabled = " << dpmd_ptr->gpu_enabled << std::endl;
    // lout << "GPU id      = " << dpmd_ptr->gpu_id << std::endl;
    // torch::Device device(torch::kCUDA, dpmd_ptr->gpu_id);
    
    // if (!gpu_enabled) {
    //   device = torch::Device(torch::kCPU);
    // }
    
    // ================================================================
    // Unit conversion
    //   DeePMD: energy in eV, forces in eV/Angstrom
    // ================================================================
    const double conv_e_factor = ONIKA_QUANTITY( 1. * eV ).convert();
    const double conv_f_factor = ONIKA_QUANTITY( 1. * eV/ang ).convert();

    // ================================================================
    // Energy: E_0 (atomic energy of the central atom only)
    // ================================================================
    _en = e[0] * conv_e_factor;

    // ================================================================
    // Force on central atom (index 0): -dE_0/dr_0
    // ================================================================
    _fx = f[0] * conv_f_factor;
    _fy = f[1] * conv_f_factor;
    _fz = f[2] * conv_f_factor;

    // ================================================================
    // Distribute partial forces to neighbors WITHIN rcut.
    //
    // f[3*(j+1) .. 3*(j+1)+2] = -dE_0 / dr_j
    //
    // Only neighbors within rcut get a non-negligible contribution:
    // E_0's descriptor smoothly vanishes at rcut, so atoms beyond
    // rcut have (essentially) zero force from E_0.
    // Those outer-shell atoms (rcut < r < 2*rcut) exist only to
    // provide complete environments for the inner neighbors.
    //
    // When every atom in the system has been processed as a center,
    // the total force on any atom j is:
    //   F_j = sum_over_all_centers_i { -dE_i / dr_j }
    //       = -dE_total / dr_j     (correct!)
    // ================================================================
    for (int j = 0; j < nn; ++j)
    {
      size_t cell_b = 0, p_b = 0;
      buf.nbh.get(j, cell_b, p_b);
      // concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double>
      //   (
      //    locks[cell_b][p_b],
      //    cells[cell_b][field::fx][p_b],
      //    cells[cell_b][field::fy][p_b],
      //    cells[cell_b][field::fz][p_b],
      //    f[3 * (j + 1) + 0] * conv_f_factor,
      //    f[3 * (j + 1) + 1] * conv_f_factor,
      //    f[3 * (j + 1) + 2] * conv_f_factor
      //    );
      ONIKA_CU_ATOMIC_ADD( cells[cell_b][field::fx][p_b], f[3 * (j + 1) + 0] * conv_f_factor );
      ONIKA_CU_ATOMIC_ADD( cells[cell_b][field::fy][p_b], f[3 * (j + 1) + 1] * conv_f_factor );
      ONIKA_CU_ATOMIC_ADD( cells[cell_b][field::fz][p_b], f[3 * (j + 1) + 2] * conv_f_factor );      

    }

    // ================================================================
    // Virial contribution from E_0 (if requested)
    // ================================================================
    // if constexpr (compute_virial) {
    //   _vir.xx = v[0] * conv_e_factor;
    //   _vir.xy = v[1] * conv_e_factor;
    //   _vir.xz = v[2] * conv_e_factor;
    //   _vir.yx = v[3] * conv_e_factor;
    //   _vir.yy = v[4] * conv_e_factor;
    //   _vir.yz = v[5] * conv_e_factor;
    //   _vir.zx = v[6] * conv_e_factor;
    //   _vir.zy = v[7] * conv_e_factor;
    //   _vir.zz = v[8] * conv_e_factor;
    // }

    // ================================================================
    // Accumulate central atom's energy and force
    // ================================================================
    // concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double>
    //   (
    //    lock_a,
    //    fx,
    //    fy,
    //    fz,
    //    en,
    //    _fx,
    //    _fy,
    //    _fz,
    //    _en
    //    );
    ONIKA_CU_ATOMIC_ADD( fx, _fx );
    ONIKA_CU_ATOMIC_ADD( fy, _fy );
    ONIKA_CU_ATOMIC_ADD( fz, _fz );
    ONIKA_CU_ATOMIC_ADD( en, _en );
    
  }
};
