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
#include <unordered_map>
#include <cstdint>

namespace exaStamp
{
  using namespace exanb;

  // DXA pipeline step (iii) output (Stukowski, Bulatov, Arsenlis 2012): the ideal lattice vector
  // assigned to each edge of a DelaunayTessellation (see compute_dxa_edge_vectors.cpp), against a
  // single user-chosen target structure (e.g. BCC) -- an edge whose endpoints don't both locally
  // match that target structure (PTM's own per-particle classification), or whose bond direction
  // doesn't cleanly snap to one of the target lattice's ideal neighbor directions, is left
  // unresolved (this is exactly where step (iv) expects to find defects/grain boundaries).
  //
  // Edges are deduplicated across the tessellation (many tetrahedra share an edge): one entry per
  // unique (v0,v1) pair, v0<v1, in DelaunayTessellation's own (compacted) vertex numbering.
  struct DXAEdgeVectors
  {
    std::vector<std::array<uint32_t,2>> edges; // v0<v1
    std::vector<Vec3d> ideal_vector;           // target structure's own ideal/reference lattice coordinates (PTM's template units, not physical length) -- meaningful only where resolved[i]
    std::vector<bool> resolved;

    // (uint64_t(v0)<<32)|v1 -> index into the arrays above -- lets step (iv) look up a
    // tetrahedron's 6 edges by vertex pair without a linear scan.
    std::unordered_map<uint64_t,uint32_t> edge_index;

    // per-vertex (not per-edge): does this DelaunayTessellation vertex's own particle locally
    // match target_structure (ptm_type)? Parallel to DelaunayTessellation::vertices. This is what
    // compute_dxa_tet_classification actually uses to decide good/bad -- see that struct's own
    // comment for why (matches both reference DXA implementations checked against: "good" is a
    // per-atom structure-type question, not an edge-vector-closure one).
    std::vector<uint8_t> vertex_matches_target;

    static inline uint64_t key(uint32_t v0, uint32_t v1) noexcept
    {
      return v0 < v1 ? ( (uint64_t(v0)<<32) | uint64_t(v1) ) : ( (uint64_t(v1)<<32) | uint64_t(v0) );
    }
  };

  // DXA pipeline step (iv) output (see compute_dxa_tet_classification.cpp): whether each
  // DelaunayTessellation tetrahedron is part of the undistorted target-structure lattice
  // ("good") or not ("bad" -- a defect: dislocation core, grain boundary, stacking fault, second
  // phase, ...). "Good" means all 4 of its vertices individually match target_structure
  // (DXAEdgeVectors::vertex_matches_target) -- NOT an edge-vector-closure test (that was tried
  // first and found to overclassify a dislocation's ordinary elastic strain field as defective;
  // see compute_dxa_tet_classification.cpp's own comment). Parallel to
  // DelaunayTessellation::tetrahedra. 1.0=good, 0.0=bad (double, not bool, so it can be written
  // straight out as VTK CellData with no conversion -- same convention as ptm_type).
  struct DXATetClassification
  {
    std::vector<double> good;
  };

  // DXA pipeline steps (vi)-(vii) output (see compute_dxa_burgers_circuits.cpp): the confirmed
  // dislocation-core edges found by tracing a Burgers circuit -- the ring of tetrahedra sharing
  // each tessellation edge -- and summing their (already-resolved, step iii) ideal lattice
  // vectors around the ring. A closure failure (non-zero sum) above min_burgers_norm is a real
  // Burgers vector; this is also how "bad candidates" (isolated misclassified tetrahedra that
  // don't actually enclose a topological defect) get filtered out -- see
  // compute_dxa_burgers_circuits.cpp's DXATetClassification refinement.
  struct DXABurgersCircuits
  {
    std::vector<std::array<uint32_t,2>> edges; // v0<v1, subset of DXAEdgeVectors::edges
    std::vector<Vec3d> burgers_vector;         // target structure's own ideal/reference lattice coordinates, same units as DXAEdgeVectors::ideal_vector
  };
}
