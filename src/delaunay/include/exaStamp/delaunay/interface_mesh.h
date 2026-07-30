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

#include <vector>
#include <array>
#include <unordered_map>
#include <cstdint>

namespace exaStamp
{
  // DXA pipeline step (v) output (see compute_interface_mesh.cpp): the 2D surface separating
  // "good" from "bad" DelaunayTessellation tetrahedra (DXATetClassification) -- one triangle per
  // tetrahedron face where a good tet and a bad tet meet. Vertex indices reference
  // DelaunayTessellation::vertices directly (same numbering, no separate vertex list).
  //
  // Only faces where BOTH tets sharing them were actually kept in this rank's own
  // DelaunayTessellation are considered -- a face where the neighboring tet wasn't retained
  // locally (this rank's own tessellation boundary, a domain-decomposition artifact, same
  // "ghost-fringe-trust" concern compute_delaunay.cpp already reasons about) is deliberately left
  // out rather than guessed at, since its far side's good/bad status isn't actually known here.
  struct InterfaceMesh
  {
    std::vector<std::array<uint32_t,3>> triangles; // vertex order oriented so the right-hand-rule normal points outward, away from the bad tet's own interior (see compute_interface_mesh.cpp)
    std::vector<uint32_t> good_tet;                // per triangle: index into DelaunayTessellation::tetrahedra, the "good" side
    std::vector<uint32_t> bad_tet;                 // per triangle: index into DelaunayTessellation::tetrahedra, the "bad" side

    // interface mesh's own edge adjacency (v0<v1 key -> triangle indices sharing that edge).
    // Normally exactly 2 (a closed 2-manifold, since a dislocation line can't end inside a
    // perfect crystal) -- a count of 1 marks a domain-decomposition cutoff (the interface surface
    // is genuinely cut off by the edge of this rank's own kept-tet data, not a real edge of the
    // defect), and >2 would mark a branch point. Future steps (vi-ix: circuit search / sweep)
    // need this to walk the surface; not resolved into a strict half-edge (twin/next/prev)
    // structure yet since that traversal pattern isn't scoped out yet.
    std::unordered_map<uint64_t, std::vector<uint32_t>> edge_triangles;

    static inline uint64_t edge_key(uint32_t v0, uint32_t v1) noexcept
    {
      return v0 < v1 ? ( (uint64_t(v0)<<32) | uint64_t(v1) ) : ( (uint64_t(v1)<<32) | uint64_t(v0) );
    }
  };
}
