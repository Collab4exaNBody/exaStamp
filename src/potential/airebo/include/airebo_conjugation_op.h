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
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>
#include <exanb/core/concurent_add_contributions.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include "airebo_params.h"
#include "airebo_force_op_cached.h"  // rebo_Sp

namespace exaStamp
{
  using namespace exanb;

  // ===========================================================================
  // Operator to compute per atom bond order
  // ===========================================================================

  struct AireboNconjOpExtStorage
  {
    int    itype = 0;    // central atom type, cached at ContextStart
    double nconj  = 0.0;
  };

  
  template<class NijcFieldT, class NijhFieldT, class NconjFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) AireboNconjOp
  {
    const AireboParamsRO* m_params  = nullptr;
    NijcFieldT  m_nijc_field  = {};   // read neighbour's nC
    NijhFieldT  m_nijh_field  = {};   // read neighbour's nH
    NconjFieldT m_nconj_field = {};   // write central atom's Nconj

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart ) const
    {
      ctx.ext.itype = cells[cell_a][field::type][p_a];
      ctx.ext.nconj  = 0.0;
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
    {
      cells[cell_a][m_nconj_field][p_a] = ctx.ext.nconj;
    }
    
    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator() (
      ComputeBufferT& ctx
      , const Vec3d& dr, double d2
      , CellParticlesT cells,size_t cell_b,size_t p_b
      , double /*scale */ ) const
    {
      // Only C neighbours contribute to Nconj
      const int type_b = cells[cell_b][field::type][p_b];
      if (type_b != 0) return;

      const int    type_a = ctx.ext.itype;
      const double r      = sqrt(d2);

      // w_ij = Sp(r, rcmin[i][j], rcmax[i][j])
      double dw;
      const double w_ij = rebo_Sp(r, m_params->rcmin[type_a][type_b],
                                      m_params->rcmax[type_a][type_b], dw);
      if (w_ij <= 0.0) return;

      // Nj_excl_i = nC[j] + nH[j] - w_ij
      // nC[j]+nH[j] includes i's contribution (w_ij); subtract it to get
      // j's coordination excluding i, as required by the LAMMPS formula.
      const double Nj = cells[cell_b][m_nijc_field][p_b]
                      + cells[cell_b][m_nijh_field][p_b]
                      - w_ij;

      double dSpN;
      const double SpN = rebo_Sp(Nj, m_params->Nmin, m_params->Nmax, dSpN);

      ctx.ext.nconj += w_ij * SpN;
    }   
  };

} // namespace exaStamp

// ComputePairTraits specialization must live in the exanb namespace
// where the primary template is defined.
namespace exanb
{
  template<class NijcFieldT, class NijhFieldT, class NconjFieldT>
  struct ComputePairTraits< exaStamp::AireboNconjOp<NijcFieldT, NijhFieldT, NconjFieldT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };
} // namespace exanb
