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

// Launcher of staged onika kernels: block(n, f, tag) runs f(0 .. n-1) a block per item, flat(n, f, tag) a thread per
// item, through onika (GPU, else OpenMP) with the execution contexts of the calling operator. Drivers take the launcher
// as a template parameter, so that a debug build can run them through another one (mlip_onika_emul.h).

#include <onika/parallel/parallel_for.h>
#include <onika/parallel/block_parallel_for.h>

#include <cstddef>

namespace exaStamp
{
  namespace mlip_utils
  {
    template<class OpT>
    struct OnikaLauncher
    {
      OpT* op = nullptr;
      template<class F> inline void block(size_t n, const F& f, const char* tag) const
      {
        if( n > 0 ) onika::parallel::block_parallel_for( n , f , op->parallel_execution_context(tag) );
      }
      template<class F> inline void flat(size_t n, const F& f, const char* tag) const
      {
        if( n > 0 ) onika::parallel::parallel_for( n , f , op->parallel_execution_context(tag) );
      }
    };
  }
}
