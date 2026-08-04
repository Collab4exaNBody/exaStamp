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
#include <cstdint>

namespace exaStamp
{
  // Per-interface-mesh-triangle ownership: which LOCAL (per-rank, pre-MPI-stitch) dislocation
  // claimed each triangle during compute_dxa_circuit_sweep's own growth -- a direct, read-only
  // copy of that operator's own internal `facet_owner` array, resolved through the same-rank merge
  // chain (resolve_node) to the FINAL per-rank segment/dislocation numbering, matching
  // DXADislocationLines::dislocation_id at THIS rank BEFORE MPI stitching (compute_dxa_mpi_stitch_
  // lines' own dislocation_id renumbering happens later and is a different numbering -- see
  // compute_dxa_mpi_stitch_lines.cpp's own remap-broadcast mechanism for how the two get connected).
  // Parallel to InterfaceMesh::triangles. -1 = never claimed by any circuit (this triangle's own
  // bad tet is either not part of any traced dislocation, or belongs to noise/an untraced defect).
  struct DXATriangleOwnership
  {
    std::vector<int32_t> triangle_dislocation_id;
  };

  // Per-Delaunay-tetrahedron ownership: which LOCAL (per-rank, pre-MPI-stitch) dislocation each
  // "bad" (defective) tet belongs to, extending outward from the traced interface-mesh boundary
  // (DXATriangleOwnership) via a multi-source breadth-first flood through bad-tet-to-bad-tet
  // adjacency -- see compute_dxa_core_atoms.cpp's own header comment for the full mechanism and
  // why a plain unrestricted flood-fill from one dislocation's own boundary alone would incorrectly
  // conflate neighboring dislocations' own core regions wherever their bad-tet volumes touch (e.g.
  // near a real junction). Parallel to DelaunayTessellation::tetrahedra. -1 = good (never
  // defective) tet, OR a bad tet not reached (in the flood) by any traced dislocation's own
  // boundary claim.
  struct DXACoreTetOwnership
  {
    std::vector<int32_t> tet_dislocation_id;
  };

  // Sent by compute_dxa_mpi_stitch_lines back to EVERY rank (not just rank 0, unlike everything
  // else that operator produces): translates THIS rank's own local dislocation_id numbering
  // (matching compute_dxa_circuit_sweep's/compute_dxa_core_atoms' own local numbering, BEFORE MPI
  // stitching) into the FINAL, post-stitch dislocation id -- or -1 if this rank's own local
  // dislocation didn't survive into the final result (e.g. dropped as a redundant ghost-overlap
  // duplicate with no valid resolution). Indexed by local dislocation_id; size = however many local
  // dislocations THIS rank's own compute_dxa_circuit_sweep produced.
  struct DXALocalToFinalDislocationId
  {
    std::vector<int32_t> final_id;
  };
}
