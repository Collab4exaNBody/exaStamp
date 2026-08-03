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

// Plain-C++ port of OVITO's real DXA cluster-graph objects
// (ovito/src/ovito/crystalanalysis/objects/{Cluster,ClusterGraph,ClusterVector}.h): a "cluster" is
// a maximal contiguous group of atoms sharing one consistent discrete lattice-orientation labeling
// (see lattice_structure.h) -- roughly one grain. A "cluster transition" is the rotation/reflection
// relating two adjacent clusters' own frames (e.g. across a grain boundary or a low-angle
// sub-boundary), forming a graph whose nodes are clusters and whose edges let a vector be
// transported from one cluster's frame to another's.
//
// Simplified relative to OVITO's own version: uses plain std::vector storage instead of intrusive
// linked lists + a custom memory pool (an OVITO-wide performance idiom for its object graph, not
// something this host-only, moderate-cluster-count pipeline needs -- ponytail: a vector is simpler
// and this isn't a hot loop, unlike compute_ptm's per-atom inner loop). No `product`/
// `inverseProduct`/preferred-orientation machinery either, same reasoning as lattice_structure.h.
namespace exaStamp
{
  using namespace exanb;

  struct ClusterTransition
  {
    int cluster1 = -1; // transforms a vector FROM this cluster's frame...
    int cluster2 = -1; // ...TO this cluster's frame
    Mat3d tm = onika::math::make_identity_matrix();
    int reverse_index = -1; // index into ClusterGraph::transitions of the inverse transition
    int distance = 0;       // graph distance between cluster1 and cluster2 (0 for a self-transition)
    int area = 0;           // number of atom-atom bonds observed crossing this transition

    inline bool is_self_transition() const { return cluster1 == cluster2 && distance == 0; }

    inline Vec3d transform( const Vec3d& v ) const { return is_self_transition() ? v : (tm * v); }
  };

  struct Cluster
  {
    int id = 0;
    int structure_type = 0; // LatticeStructureType
    int atom_count = 0;
    Mat3d orientation = onika::math::make_identity_matrix(); // least-squares fit, cluster frame -> simulation frame (display/diagnostic only, not used by the elastic mapping itself)
    std::vector<int> transition_indices;                     // this cluster's own outgoing transitions (indices into ClusterGraph::transitions)
  };

  // Cluster 0 is always the reserved "null cluster" (unresolved/no match), matching OVITO's own
  // convention, and is created automatically by the default constructor.
  struct ClusterGraph
  {
    std::vector<Cluster> clusters;             // clusters[id] == the cluster with that id (dense, id is the index)
    std::vector<ClusterTransition> transitions; // both directions of every non-self transition, plus self-transitions

    ClusterGraph() { clusters.push_back( Cluster{ 0, 0, 0, onika::math::make_identity_matrix(), {} } ); }

    inline Cluster& cluster( int id ) { return clusters[id]; }
    inline const Cluster& cluster( int id ) const { return clusters[id]; }

    // Creates a new cluster and returns its id.
    int create_cluster( int structure_type );

    // Returns the (possibly newly-created) self-transition for a cluster.
    int self_transition( int cluster_id );

    // Creates a transition A->B (and its reverse B->A) unless an equal one already exists.
    // Returns the index of the A->B transition.
    int create_transition( int cluster_a, int cluster_b, const Mat3d& tm, int distance, double tm_epsilon = 1e-4 );

    // Finds the transformation matrix from cluster_a's frame to cluster_b's frame, searching the
    // graph up to distance 2 (OVITO's own hardcoded search depth -- see ClusterGraph::
    // determineClusterTransition) and caching the result as a new direct edge. Returns -1 if the
    // two clusters are in different, disconnected components.
    int determine_transition( int cluster_a, int cluster_b );
  };

  // A vector expressed in one cluster's own local (reference/ideal-lattice) frame -- the (0,0,0)
  // vector is the only one allowed to have no associated cluster.
  struct ClusterVector
  {
    Vec3d vec{0.0, 0.0, 0.0};
    int cluster = -1;
  };
}
