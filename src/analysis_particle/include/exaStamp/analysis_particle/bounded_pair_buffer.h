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

// Compute-pair buffers that cannot overflow (the default one has no overflow check in release builds):
// - NearestNeighborsComputeBuffer<K>: the K nearest neighbors, sorted by distance
// - OverflowFlagComputeBuffer: all neighbors while there is room, then sets buf.ext.overflow
// Both keep count < MaxNeighbors.

namespace exaStamp
{
  using namespace exanb;

  // sorted insertion, the farthest neighbor is dropped when full
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

  // the functor must read and clear buf.ext.overflow
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

  template<bool UseNeighbors>
  static inline constexpr size_t overflow_flag_buffer_capacity = OverflowFlagComputeBuffer<UseNeighbors>::MaxNeighbors - 1;
}
