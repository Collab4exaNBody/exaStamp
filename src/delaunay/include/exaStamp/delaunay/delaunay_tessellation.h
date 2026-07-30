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

  // This rank's share of a 3D Delaunay tessellation: only tetrahedra "owned" by this rank are
  // kept (centroid falls in one of this rank's owned, non-ghost cells -- see compute_delaunay.cpp),
  // so combining every rank's tessellation gives the whole system with neither gaps nor
  // duplicates, same owned/ghost partition already used everywhere else in this codebase.
  // vertices/tetrahedra use a compacted local numbering (only points touched by an owned
  // tetrahedron are kept), independent from the particle grid's own indexing.
  struct DelaunayTessellation
  {
    std::vector<Vec3d> vertices;
    std::vector<std::array<uint32_t,4>> tetrahedra;
    // vertex_particle_index[v] is the flat particle index (grid->cell_particle_offset_data()
    // convention: cell_particle_offset[cell]+particle) of vertices[v] -- lets a consumer (e.g.
    // write_delaunay_vtk's color_field) look up any per-particle grid field for a tet's vertices,
    // parallel to `vertices`.
    std::vector<uint32_t> vertex_particle_index;
  };
}
