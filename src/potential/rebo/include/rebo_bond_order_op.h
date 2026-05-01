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

#include "rebo_params.h"
#include "rebo_force_op_cached.h"

namespace exaStamp
{
  using namespace exanb;

  // ===========================================================================
  // Operator to compute per atom bond order
  // ===========================================================================

  struct ReboBondOrderOpExtStorage
  {
    int    itype = 0;    // central atom type, cached at ContextStart
    double nijc  = 0.0;
    double nijh  = 0.0;
  };

  template<class NijcFieldT, class NijhFieldT>
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) ReboBondOrderOp
  {
    const ReboParamsRO* m_params = nullptr;
    NijcFieldT m_nijc_field = {};
    NijhFieldT m_nijh_field = {};

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStart ) const
    {
      ctx.ext.itype = cells[cell_a][field::type][p_a];
      ctx.ext.nijc  = 0.0;
      ctx.ext.nijh  = 0.0;
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
    {
      cells[cell_a][m_nijc_field][p_a] = ctx.ext.nijc;
      cells[cell_a][m_nijh_field][p_a] = ctx.ext.nijh;
    }
    
    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC ONIKA_ALWAYS_INLINE void operator() (
      ComputeBufferT& ctx
      , const Vec3d& dr, double d2
      , CellParticlesT cells,size_t cell_b,size_t p_b
      , double /*scale */ ) const
    {
      const int    type_b = cells[cell_b][field::type][p_b];
      const int    type_a = ctx.ext.itype;
      const double r      = sqrt(d2);

      // REBO switching weight w_ij = Sp(r, rcmin[i][j], rcmax[i][j])
      double dw;
      const double w = rebo_Sp(r,
                                m_params->rcmin[type_a][type_b],
                                m_params->rcmax[type_a][type_b], dw);

      if (type_b == 0) ctx.ext.nijc += w;   // C neighbour
      else             ctx.ext.nijh += w;    // H neighbour
    }   
  };

} // namespace exaStamp

// ComputePairTraits specialization must live in the exanb namespace
// where the primary template is defined.
namespace exanb
{
  template<class NijcFieldT, class NijhFieldT>
  struct ComputePairTraits< exaStamp::ReboBondOrderOp<NijcFieldT, NijhFieldT> >
  {
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;
    static inline constexpr bool CudaCompatible          = true;
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContext      = true;
    static inline constexpr bool HasParticleContextStop  = true;
  };
} // namespace exanb
