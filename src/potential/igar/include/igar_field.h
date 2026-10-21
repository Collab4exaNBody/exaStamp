/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/
#pragma once

#include <cstdint>
#include <onika/math/basic_types_def.h>
#include <onika/memory/allocator.h>

namespace exaStamp
{

  using namespace exanb;

  // IGAR scalar field over the whole simulation domain, replicated on every MPI process (static: never exchanged
  // nor migrated). Value n = i + dims.i * (j + dims.j * k) sits at the center of voxel (i,j,k), the voxels splitting
  // the domain bounds evenly. Stored either as bytes (8-bit images: value = byte * value_scale, exact) or as doubles.
  struct IgarField
  {
    IJK m_dims = {0, 0, 0};
    double m_value_scale = 1.0;
    onika::memory::CudaMMVector<uint8_t> m_u8;
    onika::memory::CudaMMVector<double> m_f64;

    inline size_t number_of_values() const { return size_t(m_dims.i) * size_t(m_dims.j) * size_t(m_dims.k); }
    inline bool is_uint8() const { return !m_u8.empty(); }
    inline bool empty() const { return m_u8.empty() && m_f64.empty(); }
  };

} // namespace exaStamp
