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
#include <cstdint>

namespace exaStamp
{
  using namespace exanb;

  // DXA pipeline steps (viii)-(ix) output (see compute_dxa_dislocation_lines.cpp): a set of 1D
  // curves extracted from the disordered ("hole") atom population -- each curve is either an open
  // segment bounded by two junctions/endpoints, or a closed loop. Vertex indices reference
  // DelaunayTessellation::vertices directly (same numbering as InterfaceMesh/DXAEdgeVectors).
  struct DXADislocationLines
  {
    // one entry per extracted line. lines[i] is an ordered walk of DelaunayTessellation vertex
    // indices along that line; for a closed loop (is_loop[i]==1) the first and last entries are
    // the same vertex (an explicit closure), for an open segment they aren't. Populated by the
    // atom-graph-based extractors (compute_dxa_dislocation_lines, compute_dxa_mesh_dislocation_
    // lines); left empty by compute_dxa_circuit_sweep, whose own line vertices are synthetic
    // swept-circuit centroids with no corresponding tessellation vertex -- see line_positions
    // below for that case instead. A consumer should check line_positions first, falling back to
    // resolving `lines` through DelaunayTessellation::vertices only if it's empty.
    std::vector<std::vector<uint32_t>> lines;

    // same shape as `lines` (one entry per line), but real 3D positions directly rather than
    // vertex indices -- only populated by compute_dxa_circuit_sweep (each point is a swept
    // circuit's own center of mass, not any single atom's position). Empty for the other
    // extractors.
    std::vector<std::vector<Vec3d>> line_positions;

    // same shape as `line_positions` -- per-point weight, the sweeping circuit's own vertex count
    // (loop size) at the moment that point was recorded. Populated by compute_dxa_circuit_sweep
    // alongside line_positions (raw, one entry per sweep move); consumed and collapsed down to
    // match the reduced point count by smooth_dxa_dislocation_lines' coarsening pass (still
    // parallel to the coarsened line_positions afterward). This is OVITO's own
    // DislocationSegment::coreSize -- a wider/more stretched circuit (larger loop, typically near
    // a junction) means a noisier point, weighted accordingly during coarsening (see that
    // operator's own header comment). Empty for the other (non-circuit-sweep) extractors.
    std::vector<std::vector<int32_t>> core_size;

    std::vector<uint8_t> is_loop;

    // one Burgers vector per line (target structure's ideal/reference lattice coordinates, same
    // units as DXAEdgeVectors::ideal_vector / DXABurgersCircuits::burgers_vector) -- averaged from
    // nearby confirmed signal edges, see compute_dxa_dislocation_lines.cpp. Zero if no signal edge
    // was found near this particular line (shouldn't normally happen, but not fatal if it does).
    std::vector<Vec3d> burgers_vector;

    // DelaunayTessellation vertex indices of atoms where 3+ line segments meet (skeleton degree
    // >= 3) -- every open line's two endpoints are either one of these, or a genuine dead end
    // (skeleton degree 1, e.g. a domain-decomposition cutoff or an incomplete/too-aggressive
    // erosion, see min_core_depth in compute_dxa_dislocation_lines.cpp).
    std::vector<uint32_t> junction_vertices;

    // one entry per line (only populated by compute_dxa_circuit_sweep): which physical dislocation
    // this line/segment belongs to, after merging clean two-way (pass-through, not a real branch)
    // junction chains -- topologically, one physical dislocation can be swept out as several
    // separate segments (see compute_dxa_circuit_sweep.cpp's own header comment for exactly how
    // "clean chain" vs. "real branch" is distinguished). Segments sharing the same dislocation_id
    // should be summed together for length statistics; a segment whose id doesn't repeat elsewhere
    // is a complete dislocation on its own (or a real multi-way branch arm, kept separate on
    // purpose). Empty for the other extractors.
    std::vector<int32_t> dislocation_id;
  };
}
