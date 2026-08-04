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

#include <exaStamp/grains/grain_clusters.h>
#include <cstdint>

namespace exaStamp
{
  // Pure algorithm (no grid/MPI dependency, unit-testable in isolation): Node-Pair-Sampling graph
  // clustering + auto-threshold dendrogram cut + min-grain-size dissolution + (optional)
  // Dijkstra-by-distance orphan-atom adoption, following OVITO's GrainSegmentationEngine1/2 -- see
  // compute_grain_clusters.cpp's own header comment for the full algorithm writeup/citations.
  //
  // All flat arrays below are indexed by LOCAL atom index (0..n_atoms-1), same convention as
  // compute_grain_bond_misorientation's own output. `bond_id` holds GLOBAL particle ids (cross-rank
  // comparable); this function resolves them to local indices itself via `global_id`, silently
  // treating a bond whose target isn't present locally (e.g. a ghost just outside this rank's own
  // ghost halo) as absent -- this is a known, not-yet-addressed limitation: a grain that spans an
  // MPI rank boundary is NOT stitched into one id across ranks (each rank clusters its own local
  // view only), same scope gap already flagged for compute_dxa_lattice_clusters' own local BFS.
  void grain_segmentation_nps(
    size_t n_atoms,
    const double * struct_type,          // PTM_MATCH_* code per atom
    const double * orientation_mat3,     // 9 doubles/atom, row-major (Mat3d)
    const uint64_t * global_id,          // per atom
    const uint64_t * bond_id,            // max_neighbors slots/atom; GRAIN_BOND_EMPTY(UINT64_MAX) = unused slot
    const double * bond_distance,        // same shape, real Euclidean length
    const double * bond_disorientation,  // same shape, degrees; <0 = not a clustering candidate
    const int * bond_count,              // per atom, valid entries in the slots above
    int max_neighbors,
    bool auto_threshold,
    double manual_threshold_log,         // only used if !auto_threshold -- OVITO's own internal log-distance unit, NOT degrees (see this file's own header note on GraphClusteringManual)
    long min_grain_atom_count,
    bool orphan_adoption,
    unsigned int color_seed,
    GrainSegmentationResult & result );
}
