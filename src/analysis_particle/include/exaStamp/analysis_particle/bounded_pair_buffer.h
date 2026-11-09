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

#include <exanb/compute/compute_pair_buffer.h>
#include <onika/math/basic_types.h>
#include <cstddef>

// Compute-pair buffers that cannot overflow, for the structure analysis operators.
//
// exaNBody's default ComputePairBuffer2 has MAX_PARTICLE_NEIGHBORS (512) slots, and its overflow
// check (check_buffer_overflow) is compiled out when NDEBUG is defined. In a release build, an rcut
// that puts more than 512 neighbors around a particle writes past the end of the buffer and
// crashes with no message. The two buffers below make that impossible:
// - NearestNeighborsComputeBuffer<K> keeps only the K nearest neighbors, sorted by distance (for
//   operators that only look at the K nearest, e.g. PTM and CNA). Any rcut is safe.
// - OverflowFlagComputeBuffer<UseNeighbors> keeps every neighbor while there is room, and sets
//   buf.ext.overflow when one had to be dropped, so the operator can report a clear error.
// Both store at most MaxNeighbors-1 entries, so count never reaches the debug overflow check.

namespace exaStamp
{
  using namespace exanb;

  // Insertion into a buffer sorted by increasing d2. When the buffer is full, a neighbor farther
  // than the current farthest one is ignored, otherwise it replaces it.
  struct NearestNeighborsAppendFunc
  {
    template<class ComputeBufferT, class FieldArraysT, class NbhDataT=double>
    inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2,
                             FieldArraysT, size_t cell_b, size_t p_b, const NbhDataT& = 1.0 ) const noexcept
    {
      static constexpr int K = static_cast<int>( ComputeBufferT::MaxNeighbors ) - 1;
      const int n = tab.count;
      if( n == K && d2 >= tab.d2[K-1] ) return;
      int pos = ( n < K ) ? n : K-1;
      while( pos > 0 && tab.d2[pos-1] > d2 ) { tab.copy( pos-1 , pos ); --pos; }
      tab.d2[pos] = d2;
      tab.drx[pos] = dr.x;
      tab.dry[pos] = dr.y;
      tab.drz[pos] = dr.z;
      tab.nbh.set( pos , cell_b , p_b );
      if( n < K ) ++ tab.count;
    }
  };

  // K nearest neighbors + 1 spare slot (keeps count < MaxNeighbors)
  template<size_t K>
  using NearestNeighborsComputeBuffer = ComputePairBuffer2< false, false, NoExtraStorage, NearestNeighborsAppendFunc, K+1 >;

  struct PairBufferOverflowFlag
  {
    bool overflow = false;
  };

  // Default append while there is room; afterwards drop the neighbor and raise ext.overflow.
  // The operator functor must read and clear buf.ext.overflow for each particle.
  struct OverflowFlagAppendFunc
  {
    template<class ComputeBufferT, class FieldArraysT, class NbhDataT=double>
    inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2,
                             FieldArraysT cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data = 1.0 ) const noexcept
    {
      if( static_cast<size_t>( tab.count ) + 1 >= ComputeBufferT::MaxNeighbors ) { tab.ext.overflow = true; return; }
      DefaultComputePairBufferAppendFunc{}( tab, tab.count++, dr, d2, cells, cell_b, p_b, nbh_data );
    }
  };

  template<bool UseNeighbors>
  using OverflowFlagComputeBuffer = ComputePairBuffer2< false, UseNeighbors, PairBufferOverflowFlag, OverflowFlagAppendFunc >;

  // Largest number of neighbors an OverflowFlagComputeBuffer can hold
  template<bool UseNeighbors>
  static inline constexpr size_t overflow_flag_buffer_capacity = OverflowFlagComputeBuffer<UseNeighbors>::MaxNeighbors - 1;
}
