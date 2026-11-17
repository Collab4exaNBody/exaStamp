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

// Host driver of the polynomial MLP kernels (polymlp_kernels.h): central atoms ordered by element, cut into chunks that fit the
// memory budget, the three stages per chunk.

#include <exaStamp/mlip_polymlp/polymlp_kernels.h>

#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_context.h>
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace exaStamp
{
  namespace mlip_polymlp
  {
    // forces, per-atom energies and virials added into out. gpu: the kernels run on the GPU (bv / out device
    // addressable); mem_gb: budget of the per-chunk buffers; pec: launcher (mlip_onika_launch.h or the emulator)
    template<class LaunchT>
    inline void polymlp_onika_compute( const PolymlpTables& T, PolymlpBuffers& buf, const mlip_utils::BondsView& bv, mlip_utils::OutputsView& out,
                                       bool gpu, double mem_gb, const LaunchT& pec )
    {
      using namespace polymlp_onika;
      const int inum = bv.inum;
      if( inum <= 0 ) return;

      // central atoms ordered by element (stable)
      buf.order.resize( inum );
      {
        std::vector<int> cnt( T.nel + 1 , 0 );
        for(int q=0;q<inum;q++) cnt[ bv.type[ bv.ilist[q] ] + 1 ] ++;
        for(int e=0;e<T.nel;e++) cnt[e+1] += cnt[e];
        for(int q=0;q<inum;q++) { const int i = bv.ilist[q]; buf.order[ cnt[ bv.type[i] ] ++ ] = i; }
      }

      const PolymlpView tv = make_view( T );
      const int Bmax = gpu ? ONIKA_CU_MAX_THREADS_PER_BLOCK : 1;
      size_t jmax = 0;
      if( !gpu ) for(int q=0;q<inum;q++) { const int i = bv.ilist[q]; jmax = std::max( jmax , size_t( bv.offset[i+1] - bv.offset[i] ) ); }
      const size_t nscr = gpu ? size_t(Bmax) * std::max( basis_record( tv ) , forces_record( tv ) ) : host_lane_scratch( tv, jmax );
      const size_t ns = std::max( 1 , T.nsrc_max );
      const size_t nf = std::max( 1 , T.nfeat_max );
      const size_t nchn = std::max( 1 , T.nch_max );
      const size_t per_atom = 4 * ns + 2 * nf + 4 * nchn + ( gpu ? nscr : 0 ) + 1;
      const double budget = ( mem_gb > 0 ? mem_gb : 1.0 ) * 1e9 / 8.0;
      const int chunk = std::max( 1 , std::min( inum , int( std::min( budget / double(per_atom) , 2.0e9 ) ) ) );
      const size_t cap = gpu ? size_t(chunk) : ( size_t(chunk) + kPolymlpLanes - 1 ) / kPolymlpLanes * kPolymlpLanes;
      buf.A.resize( 2 * ns * cap );
      buf.W.resize( 2 * ns * cap );
      buf.D.resize( nf * cap );
      buf.GD.resize( nf * cap );
      buf.C.resize( 2 * nchn * cap );
      buf.G.resize( 2 * nchn * cap );
      const size_t nslots = gpu ? cap : size_t( omp_get_max_threads() );
      buf.scratch.resize( nscr * nslots );

      ChunkView ch;
      static_assert( kPolymlpLanes == 8 , "host lanes: vs = log2" );
      if( gpu ) { ch.vs = 62; ch.vm = ~size_t(0); ch.V = cap; ch.sg = 0; ch.fg = 0; ch.cg = 0; }
      else { ch.vs = 3; ch.vm = kPolymlpLanes - 1; ch.V = kPolymlpLanes; ch.sg = ns * kPolymlpLanes; ch.fg = nf * kPolymlpLanes; ch.cg = nchn * kPolymlpLanes; }
      ch.order = buf.order.data();
      ch.Are = buf.A.data();
      ch.Aim = ch.Are + ns * cap;
      ch.Wre = buf.W.data();
      ch.Wim = ch.Wre + ns * cap;
      ch.D = buf.D.data();
      ch.GD = buf.GD.data();
      ch.Cre = buf.C.data();
      ch.Cim = ch.Cre + nchn * cap;
      ch.Gre = buf.G.data();
      ch.Gim = ch.Gre + nchn * cap;
      ch.scratch = buf.scratch.data();
      ch.nscratch = nscr;

      // POLYMLP_PROF=1: cumulative stage wall times (host, synchronous launches) printed to stderr after each call
      static const bool prof = std::getenv( "POLYMLP_PROF" ) != nullptr;
      static double tprof[3] = { 0.0, 0.0, 0.0 };
      using clk = std::chrono::steady_clock;
      auto t0 = clk::now();
      auto lap = [&](int k) { if( !prof ) return; const auto t1 = clk::now(); tprof[k] += std::chrono::duration<double>( t1 - t0 ).count(); t0 = t1; };
      for(int c0=0;c0<inum;c0+=chunk)
      {
        const int nc = std::min( chunk , inum - c0 );
        if( gpu ) pec.block( nc , BasisFunctor{ tv, bv, ch, c0 } , "polymlp_basis" );
        else pec.flat( nc , BasisHostFunctor{ tv, bv, ch, c0, jmax } , "polymlp_basis" );
        lap(0);
        if( gpu ) pec.flat( nc , ReadoutFunctor{ tv, bv, ch, out, c0 } , "polymlp_readout" );
        else pec.flat( ( nc + kPolymlpLanes - 1 ) / kPolymlpLanes , ReadoutBatchFunctor{ ReadoutFunctor{ tv, bv, ch, out, c0 }, nc } , "polymlp_readout" );
        lap(1);
        if( gpu ) pec.block( nc , ForcesFunctor{ tv, bv, ch, out, c0 } , "polymlp_forces" );
        else pec.flat( nc , ForcesHostFunctor{ tv, bv, ch, out, c0, jmax } , "polymlp_forces" );
        lap(2);
      }
      if( prof ) std::fprintf( stderr, "POLYMLP_PROF basis %.3f s readout %.3f s forces %.3f s\n", tprof[0], tprof[1], tprof[2] );
      if( gpu ) ONIKA_CU_DEVICE_SYNCHRONIZE();
    }

  }
}
