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

#include <onika/math/basic_types.h>
#include <vector>
#include <array>
#include <cstdint>

namespace exaStamp
{
  using namespace exanb;

  // This rank's tetrahedra (centroid in an owned cell, see compute_delaunay), with a compact local
  // numbering of the vertices they use.
  struct DelaunayTessellation
  {
    std::vector<Vec3d> vertices;
    std::vector<std::array<uint32_t,4>> tetrahedra;
    // flat particle index of each vertex (cell_particle_offset[cell] + particle)
    std::vector<uint32_t> vertex_particle_index;

    // particle id of each vertex, same on every rank
    std::vector<uint64_t> vertex_global_id;

    // false for ghost vertices (only possible with keep_ghost_tets)
    std::vector<uint8_t> vertex_is_owned;
  };
}
