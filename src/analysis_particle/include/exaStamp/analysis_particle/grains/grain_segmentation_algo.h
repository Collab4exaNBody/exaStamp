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
#include <mpi.h>
#include <vector>
#include <cstddef>

// Grain segmentation building blocks, used by compute_grain_clusters.cpp

namespace exaStamp
{
  using namespace exanb;

  // weight of a crystalline bond, disorientation in degrees
  double grain_bond_weight( double disorientation_deg );

  struct NPSEdge
  {
    int u, v;
    double w;
  };

  struct NPSMerge
  {
    int a, b;          // merged nodes
    double distance;   // wnode[a]*wnode[b]/w(a,b)
    double merge_size; // harmonic mean of the two cluster sizes, in atoms
  };

  // Node-Pair-Sampling clustering (reciprocal nearest neighbor chain). node_weight: sum of all edge weights
  // of each node, including edges left out of `edges` (bonds to other ranks). node_size: atoms per node.
  // Duplicated edges are summed. Merges are returned sorted by distance.
  std::vector<NPSMerge> nps_dendrogram( const std::vector<double>& node_weight,
                                        const std::vector<double>& node_size,
                                        const std::vector<NPSEdge>& edges );

  // Automatic merge threshold: least absolute deviation fit of log(distance) vs log(merge_size), weighted
  // by merge_size; returns the largest log(distance) with residual < 1.5 * median |residual|. Collective:
  // each rank passes its share of the merges. -infinity without merges.
  double nps_auto_threshold( MPI_Comm comm, const std::vector<double>& merge_size, const std::vector<double>& distance );

  // random per-grain colors
  std::vector<Vec3d> grain_colors( long n_grains, unsigned int seed );
}
