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

#include <exaStamp/delaunay/dxa_lattice_correspondence.h>
#include <optional>
#include <unordered_map>
#include <array>
#include <cstdint>

namespace exaStamp
{
  using namespace exanb;

  // DXA elastic-mapping step (iii), the real mechanism -- port of OVITO's CrystalPathFinder
  // (ovito/src/ovito/crystalanalysis/modifier/dxa/CrystalPathFinder.{h,cpp}): finds an atom-to-atom
  // path from atom1 to atom2 that stays entirely within the "good crystal" region (only ever
  // stepping along a classified atom's own native lattice-correspondence neighbors, or -- via the
  // reverse-neighbor-search branch -- through a single disordered intermediate atom that a
  // classified neighbor reaches out to), bounded by max_path_length hops. Robust near defects in a
  // way a two-sided endpoint match isn't: the path only ever needs ONE endpoint (or an intermediate
  // stepping stone) to be classified, and routes around genuinely bad atoms rather than refusing
  // the whole edge outright.
  //
  // Returns the ideal vector connecting atom1->atom2 (as a ClusterVector: components expressed in
  // whichever cluster's frame the path's own concatenation happened to land in -- the caller
  // re-expresses this in a specific cluster via ClusterGraph::determine_transition, exactly as
  // ElasticMapping::assignIdealVectorsToEdges does, see compute_dxa_crystal_path_edge_vectors.cpp).
  // std::nullopt if no such path exists within max_path_length.
  //
  // visited_scratch is a caller-owned, atom-count-sized buffer reused across many calls (this
  // pipeline calls this once per tessellation edge -- often hundreds of thousands of times -- so
  // allocating a fresh visited-set per call would be a real cost, unlike OVITO's own reused
  // std::vector<bool> member). Guaranteed all-zero both before and after each call (this function
  // clears every flag it set before returning, mirroring CrystalPathFinder's own explicit cleanup).
  std::optional<ClusterVector> dxa_crystal_path_find(
      const DXALatticeCorrespondence& lc, DXALatticeClusters& clusters,
      int64_t atom1, int64_t atom2, int max_path_length,
      std::vector<uint8_t>& visited_scratch );

  // Port of ElasticMapping::assignVerticesToClusters() -- a SEPARATE, easy-to-miss step from
  // buildClusters/connectClusters (both entirely atom-classification-level; this one operates on
  // the DELAUNAY TESSELLATION's own vertex adjacency, not the classified-neighbor-only graph).
  // Propagates a cluster id to EVERY tessellation vertex, including atoms that
  // compute_dxa_lattice_correspondence never classified at all -- flood-filled outward from
  // already-clustered vertices via ordinary tessellation edges (BFS by hop count here; OVITO's own
  // "repeat until no change" scan converges to an equivalent, if not always identical in ties,
  // assignment). Critical distinction from DXALatticeClusters::atom_cluster: this is used ONLY to
  // decide whether an edge is even worth attempting to resolve, and to pick which cluster's frame
  // to express its ideal vector in -- dxa_crystal_path_find() itself still uses the RAW,
  // un-propagated atom_cluster internally (its own reverse-neighbor-search branch is exactly the
  // mechanism for stepping through a genuinely unclassified atom, and needs to know it really is
  // unclassified). Missing this step was found to be the actual cause of a real gap: without it,
  // every tessellation edge touching *any* unclassified atom was rejected outright before
  // CrystalPathFinder even got a chance to route around it -- unnecessarily, since the walk can
  // usually resolve those edges just fine.
  std::vector<int32_t> dxa_propagate_vertex_clusters(
      const std::vector<std::array<uint32_t,2>>& tessellation_edges,
      const std::vector<uint32_t>& vertex_particle_index,
      const DXALatticeClusters& clusters );

  // Output of compute_dxa_crystal_path_edge_vectors (ElasticMapping::assignIdealVectorsToEdges
  // port): per deduplicated Delaunay tessellation edge (v0<v1), the ideal vector connecting the two
  // vertices -- expressed in vertex v0's own cluster frame -- plus the cluster transition from v0's
  // cluster to v1's cluster (needed by the per-tetrahedron elastic-mapping compatibility test,
  // stage 3, for its own Frank-rotation check). Unresolved (resolved[i]==false) if either vertex's
  // atom isn't part of any cluster, or CrystalPathFinder couldn't find a path within
  // crystal_path_steps.
  struct DXACrystalPathEdgeVectors
  {
    std::vector<std::array<uint32_t,2>> edges;         // v0<v1, tessellation vertex indices
    std::vector<Vec3d> ideal_vector;                    // in v0's cluster frame, meaningful only where resolved[i]
    std::vector<int32_t> cluster_transition;             // index into ClusterGraph::transitions, v0's cluster -> v1's cluster; -1 if unresolved
    std::vector<bool> resolved;
    std::unordered_map<uint64_t,uint32_t> edge_index;    // (uint64_t(v0)<<32)|v1 -> index into the arrays above

    static inline uint64_t key(uint32_t v0, uint32_t v1) noexcept
    {
      return v0 < v1 ? ( (uint64_t(v0)<<32) | uint64_t(v1) ) : ( (uint64_t(v1)<<32) | uint64_t(v0) );
    }
  };

  // DXA pipeline step (iv), the real mechanism -- port of ElasticMapping::
  // isElasticMappingCompatible() (ovito/src/ovito/crystalanalysis/modifier/dxa/ElasticMapping.cpp):
  // a tetrahedron is "good" iff the elastic mapping (DXACrystalPathEdgeVectors) is self-consistent
  // across all 4 of its faces -- a genuine per-tetrahedron Burgers-circuit-closure test (the 3 edge
  // vectors around each face must sum to zero) AND a Frank-rotation/disclination test (the product
  // of the 3 faces' own cluster transitions must be the identity) -- NOT a per-vertex or
  // per-edge-count heuristic (the old compute_dxa_tet_classification's criterion). Requires all 6
  // of the tet's edges to be resolved; may mutate `clusters`' graph (cached transitions), same
  // reason dxa_crystal_path_find does.
  bool dxa_is_elastic_mapping_compatible(
      const DXACrystalPathEdgeVectors& edge_vectors, DXALatticeClusters& clusters,
      const std::array<uint32_t,4>& tet_vertices );

  // Output of compute_dxa_elastic_mapping_tet_classification: dxa_is_elastic_mapping_compatible()
  // applied to every tetrahedron of a DelaunayTessellation, parallel to
  // DelaunayTessellation::tetrahedra. 1.0=good (part of the undistorted lattice), 0.0=bad (a
  // defect) -- double, not bool, so write_delaunay_vtk can write it straight out as CellData with
  // no conversion, same convention as the old DXATetClassification.
  struct DXAElasticMappingTetClassification
  {
    std::vector<double> good;
  };
}
