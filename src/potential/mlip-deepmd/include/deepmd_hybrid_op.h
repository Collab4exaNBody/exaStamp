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

// DeePMD hybrid descriptor (host), one central atom per call: the sub-descriptors are evaluated on the neighbors of
// the largest cutoff, each with its own selection within its rcut (neighbors farther than a sub-descriptor's rcut in
// deepmd's hybrid list have a zero switch, as empty slots), their descriptors concatenated into the fitting network,
// then each sub-descriptor's backward part from its slice of dE/ddescriptor.
// Sub-descriptors of the cooperative functor (DeepmdForceOp) run on a copy of the neighbor buffer (their selection
// reorders it), attention ones (DeepmdAttenOp) on the buffer itself.

#include <vector>
#include "deepmd_force_op.h"
#include "deepmd_atten_op.h"

namespace exaStamp
{
  using namespace exanb;

  struct DeepmdHybridOp
  {
    DeepmdModelView top;                       // fitting networks, energy bias, species map, units
    const DeepmdForceOp* fops = nullptr;       // sub-descriptors in order: kind[k] 0 -> fops[idx[k]], 1 -> aops[idx[k]]
    const DeepmdAttenOp* aops = nullptr;
    const int* kind = nullptr;
    const int* idx = nullptr;
    const int* dim = nullptr;                  // descriptor size of each
    int nsub = 0;

    template <class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }
    template <class ComputeBufferT, class CellParticlesT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3d &virial, CellParticlesT cells) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }
    template <class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t n, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, CellParticlesT cells, GridCellLocksT locks, ParticleLockT &lock_a) const
    {
      FakeMat3d virial;
      this->operator()(n, buf, en, fx, fy, fz, type, virial, cells, locks, lock_a);
    }

    template <class ComputeBufferT, class CellParticlesT, class Mat3dT, class GridCellLocksT, class ParticleLockT>
    inline void operator()(size_t nn, ComputeBufferT &buf, double &en, double &fx, double &fy, double &fz, int type, Mat3dT &virial, CellParticlesT cells, GridCellLocksT, ParticleLockT &) const
    {
      static constexpr bool compute_virial = std::is_same_v<Mat3dT, Mat3d>;
      static thread_local std::vector<ComputeBufferT> bufs;
      static thread_local std::vector<DeepmdAttenWork> works;
      static thread_local std::vector<double> desc, ddesc, g0, g1;
      static thread_local DPMlpWork<DP_MAX_FIT_LAYERS, DP_MAX_FIT_WIDTH> fitw;
      if( int(bufs.size()) < nsub ) bufs.resize(nsub);
      if( int(works.size()) < nsub ) works.resize(nsub);
      const int ti = top.species_type[type];

      // forward parts, concatenated descriptor
      int nd = 0;
      for(int k = 0; k < nsub; k++) nd += dim[k];
      if( int(desc.size()) < nd ) { desc.resize(nd); ddesc.resize(nd); }
      int off = 0;
      for(int k = 0; k < nsub; k++)
      {
        const double* D = nullptr;
        if( kind[k] == 0 ) { bufs[k] = buf; fops[ idx[k] ].forward_desc(nn, bufs[k], type); D = bufs[k].ext.desc; }
        else { aops[ idx[k] ].forward(nn, buf, type, works[k]); D = works[k].desc.data(); }
        for(int c = 0; c < dim[k]; c++) desc[off + c] = D[c];
        off += dim[k];
      }

      // fitting network, dE/ddescriptor
      const DPNet& fit = top.fitting[ top.fitting_mixed_types ? 0 : ti ];
      const double e_atom = dp_mlp_forward( fit, top.layers, top.weights, desc.data(), fitw )[0] + top.energy_bias[ti];
      const size_t gw = std::max( nd , DP_MAX_FIT_WIDTH );
      if( g0.size() < gw ) { g0.resize(gw); g1.resize(gw); }
      const double one = 1.0;
      dp_mlp_backward( fit, top.layers, top.weights, desc.data(), fitw, &one, ddesc.data(), g0.data(), g1.data() );

      // backward parts
      Vec3d fi = { 0.0, 0.0, 0.0 };
      Mat3d vir = { 0.,0.,0., 0.,0.,0., 0.,0.,0. };
      off = 0;
      for(int k = 0; k < nsub; k++)
      {
        if( kind[k] == 0 )
        {
          auto& x = bufs[k].ext;
          for(int c = 0; c < dim[k]; c++) x.desc[c] = ddesc[off + c];
          fops[ idx[k] ].template backward_desc<compute_virial>(bufs[k], type, cells);
          fi = fi + Vec3d{ x.fi[0], x.fi[1], x.fi[2] };
          if constexpr ( compute_virial ) vir = vir + Mat3d{ x.vir[0], x.vir[1], x.vir[2], x.vir[3], x.vir[4], x.vir[5], x.vir[6], x.vir[7], x.vir[8] };
        }
        else
        {
          DeepmdAttenWork& W = works[k];
          for(int c = 0; c < dim[k]; c++) W.ddesc[c] = ddesc[off + c];
          aops[ idx[k] ].template backward<compute_virial>(buf, type, cells, W, fi, vir);
        }
        off += dim[k];
      }
      atomic_add_contribution(fx, fi.x);
      atomic_add_contribution(fy, fi.y);
      atomic_add_contribution(fz, fi.z);
      en += e_atom * top.conv_energy;
      if constexpr ( compute_virial ) { virial += vir; }
    }
  };

}

namespace exanb
{
  template <> struct ComputePairTraits<exaStamp::DeepmdHybridOp>
  {
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
    static inline constexpr bool CudaCompatible = false;   // host only
    static inline constexpr bool HasParticleContextStart = false;
    static inline constexpr bool HasParticleContextStop = false;
    static inline constexpr bool HasParticleContext = false;
    static inline constexpr bool RequiresNbhOptionalData = true;
    static inline constexpr bool ComputeEmptyNeighborhood = true;
  };
}
