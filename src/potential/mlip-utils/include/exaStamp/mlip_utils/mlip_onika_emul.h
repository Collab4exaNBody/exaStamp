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

// Debug: run onika block kernels with GPU semantics on the host (a team of std::thread per block, synchronised by a
// barrier), to test the GPU code paths of the staged kernels without a GPU. Include this header in a dedicated .cpp,
// after the onika headers and before the kernel headers, with the kernel namespaces renamed by macros so that the
// emulated build does not clash with the regular one. The drivers then run with gpu = true and EmulLauncher{nthreads}.
// Slow: for small systems only.

#include <onika/cuda/cuda.h>
#include <onika/parallel/parallel_for.h>
#include <onika/parallel/block_parallel_for.h>

#include <barrier>
#include <cstddef>
#include <thread>
#include <vector>

namespace mlip_onika_emul
{
  inline int g_bs = 1;
  inline thread_local int t_tid = 0;
  inline std::barrier<>* g_bar = nullptr;
  inline void sync() { g_bar->arrive_and_wait(); }

  struct EmulLauncher
  {
    int nt = 32;
    template<class F> inline void block(size_t n, const F& f, const char*) const
    {
      if( n == 0 ) return;
      std::barrier<> bar( nt );
      g_bar = &bar;
      g_bs = nt;
      std::vector<std::thread> th;
      for(int t=0;t<nt;t++) th.emplace_back( [&,t]() { t_tid = t; for(size_t a=0;a<n;a++) { f(a); bar.arrive_and_wait(); } } );
      for(auto& x : th) x.join();
      g_bs = 1;
      g_bar = nullptr;
    }
    template<class F> inline void flat(size_t n, const F& f, const char*) const
    {
      for(size_t k=0;k<n;k++) f(k);
    }
  };
}

#undef ONIKA_CU_BLOCK_SIZE
#undef ONIKA_CU_THREAD_IDX
#undef ONIKA_CU_BLOCK_IDX
#undef ONIKA_CU_BLOCK_SYNC
#undef ONIKA_CU_VALUE_IF_CUDA
#define ONIKA_CU_BLOCK_SIZE ( ::mlip_onika_emul::g_bs )
#define ONIKA_CU_THREAD_IDX ( ::mlip_onika_emul::t_tid )
#define ONIKA_CU_BLOCK_IDX 0
#define ONIKA_CU_BLOCK_SYNC() ::mlip_onika_emul::sync()
#define ONIKA_CU_VALUE_IF_CUDA(a,b) (a)
