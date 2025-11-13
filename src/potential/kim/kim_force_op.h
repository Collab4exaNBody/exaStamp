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

struct alignas(onika::memory::DEFAULT_ALIGNMENT) KimComputeBuffer
{
  alignas(onika::memory::DEFAULT_ALIGNMENT) int type[exanb::MAX_PARTICLE_NEIGHBORS];
};

// functor that populate compute buffer's extended storage for particle charges
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

struct alignas(onika::memory::DEFAULT_ALIGNMENT) KimForceOp
{

  KIMThreadContext* m_thread_ctx = nullptr;
  std::vector<int> kim_particle_codes;
  
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
    
    Mat3dT _vir;
    double _en = 0.;
    double _fx = 0.;
    double _fy = 0.;
    double _fz = 0.;
    
    static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;
    size_t tid = omp_get_thread_num();
    KIMThreadContext & kim_ctx = m_thread_ctx[tid];
    auto kimptr = kim_ctx.kim_model;
        
    // number of particles in this local cluster: 1 (center) + n neighbors
    const int np = static_cast<int>(n) + 1;
        
    // put central at origin; neighbors are already r_ij = (drx, dry, drz)
    std::vector<double> coords(3 * static_cast<size_t>(np), 0.0);
    std::vector<int> species_codes(np, 0);
    species_codes[0] = type;
    for (int i = 0; i < static_cast<int>(n); ++i) {
      coords[3 * ( i + 1 ) + 0] = buf.drx[i];
      coords[3 * ( i + 1 ) + 1] = buf.dry[i];
      coords[3 * ( i + 1 ) + 2] = buf.drz[i];
      species_codes[i+1] = kim_particle_codes[buf.ext.type[i]];
      //          std::cout << "cx,cy,cz = " << coords[3 * i + 0] << ","<< coords[3 * i + 1] << ","<< coords[3 * i + 2] << std::endl;
    }
        
    // contributing: only central particle is a contributing particle. Other particles just serve to compute energy and force on central particle.
    std::vector<int> contributing(np, 0);
    contributing[0] = 1;
        
    // species: reuse the code you already queried into particleSpecies_cluster_model[0]

    // This block is not to be done again : should crash at kim_init
    // int isSpeciesSupported;
    // int error = kimptr->GetSpeciesSupportAndCode(KIM::SPECIES_NAME::Ta,
    //                                              &isSpeciesSupported,
    //                                              &(species_codes[0]));
    // if (error) MY_ERROR("get_species_code");

    // Defining outputs
    double localEnergy = 0.0;
    std::vector<double> energies(static_cast<size_t>(np), 0.0);
    std::vector<double> forces(3 * static_cast<size_t>(np), 0.0);
    std::vector<double> virials(6 * static_cast<size_t>(np), 0.0);

    // lightweight neighbor-list payload: ONLY central has neighbors
    struct CentralOnlyNL {
      int np;
      int n;                              // number of neighbors of central
      std::vector<int> neighbors_indices; // length n, values 1..n
    } nl;
        
    nl.np = np;
    nl.n  = static_cast<int>(n);
    nl.neighbors_indices.resize(nl.n);
    for (int j = 0; j < nl.n; ++j) nl.neighbors_indices[j] = j + 1; // neighbors of central

    // prepare ComputeArguments
    KIM::ComputeArguments* computeArguments = nullptr;
    {
      int err = kimptr->ComputeArgumentsCreate(&computeArguments);
      if (err) { MY_ERROR("KIM::ComputeArgumentsCreate() failed."); }
    }

    // wire required argument pointers
    {
      int np_local = np;
      int err =
        computeArguments->SetArgumentPointer(
                                             KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &np_local) ||
        computeArguments->SetArgumentPointer(
                                             KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes, species_codes.data()) ||
        computeArguments->SetArgumentPointer(
                                             KIM::COMPUTE_ARGUMENT_NAME::particleContributing, contributing.data()) ||
        computeArguments->SetArgumentPointer(
                                             KIM::COMPUTE_ARGUMENT_NAME::coordinates, coords.data()) ||
        computeArguments->SetArgumentPointer(
                                             KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, &localEnergy) ||
        computeArguments->SetArgumentPointer(
                                             KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, energies.data());
      computeArguments->SetArgumentPointer(
                                           KIM::COMPUTE_ARGUMENT_NAME::partialForces, forces.data());
      computeArguments->SetArgumentPointer(
                                           KIM::COMPUTE_ARGUMENT_NAME::partialVirial, virials.data());          
      if (err) { kimptr->ComputeArgumentsDestroy(&computeArguments); MY_ERROR("Error in SetArgumentPointer."); }
    }

    // install a GetNeighborList callback (non-capturing lambda decays to function ptr)
    using GetNeighSig = int(void*, int, double const*, int, int, int*, int const**);
    GetNeighSig* get_neigh_cb = +[](void* dataObject,
                                    int numberOfNeighborLists,
                                    double const* /*cutoffs*/,
                                    int neighborListIndex,
                                    int particleNumber,
                                    int* numberOfNeighbors,
                                    int const** neighborsOfParticle) -> int
    {
      //          auto* nl = static_cast<LocalNeighList*>(dataObject);
      auto* nl = static_cast<CentralOnlyNL*>(dataObject);          
      if (numberOfNeighborLists != 1) return 1;
      if (neighborListIndex != 0)     return 1;
          
      // if (particleNumber < 0 || particleNumber >= nl->numberOfParticles) return 1;

      // *numberOfNeighbors   = nl->NNeighbors[particleNumber];
      // *neighborsOfParticle = &(nl->neighborList[particleNumber * nl->numberOfParticles]);

      if (particleNumber == 0) {
        *numberOfNeighbors   = nl->n;
        *neighborsOfParticle = nl->neighbors_indices.data();
      } else {
        // No neighbor list for non-central particles
        *numberOfNeighbors   = 0;
        *neighborsOfParticle = nullptr;
      }
      return 0;
    };

    {
      int err = computeArguments->SetCallbackPointer(
                                                     KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
                                                     KIM::LANGUAGE_NAME::cpp,
                                                     reinterpret_cast<KIM::Function*>(get_neigh_cb),
                                                     &nl);
      if (err) { kimptr->ComputeArgumentsDestroy(&computeArguments); MY_ERROR("Error in SetCallbackPointer(GetNeighborList)."); }
    }

    // call the model (use the 1-argument overload as in the official example)
    {
      int err = kimptr->Compute(computeArguments);
      if (err) { kimptr->ComputeArgumentsDestroy(&computeArguments); MY_ERROR("KIM::Model::Compute() failed."); }
    }
    
    double conv_energy_factor = ONIKA_CONST_QUANTITY( 1. * eV ).convert();
    Vec3d localForce = Vec3d{forces[0],forces[1],forces[2]} * conv_energy_factor;
    _fx = localForce.x;
    _fy = localForce.y;
    _fz = localForce.z;
    _en = (localEnergy * conv_energy_factor);

    // accumulate central-particle results (index 0)
    for (int j = 0; j < n; ++j)
      {
        size_t cell_b=0, p_b=0;
        buf.nbh.get(j,cell_b,p_b);
        auto& lock_b = locks[cell_b][p_b];
        lock_b.lock();
        cells[cell_b][field::fx][p_b] += forces[3*(j+1)+0] * conv_energy_factor;
        cells[cell_b][field::fy][p_b] += forces[3*(j+1)+1] * conv_energy_factor;
        cells[cell_b][field::fz][p_b] += forces[3*(j+1)+2] * conv_energy_factor;
        lock_b.unlock();
      }

    lock_a.lock();
    en += _en;
    fx += _fx;
    fy += _fy;
    fz += _fz;
    lock_a.unlock();

    // cleanup
    {
      int err = kimptr->ComputeArgumentsDestroy(&computeArguments);
      if (err) { MY_ERROR("KIM::ComputeArgumentsDestroy() failed."); }
    }
    // -----------------------------------------------------------------------------
    
    //    std::cout << "Inside Force Operator" << std::endl;
  }      
};
    
