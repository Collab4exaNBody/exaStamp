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

#include <exaStamp/delaunay/cluster_graph.h>

#include <algorithm>
#include <cmath>

// Ported from OVITO's real DXA source, ovito/src/ovito/crystalanalysis/objects/ClusterGraph.cpp --
// see cluster_graph.h for what's simplified relative to the original.
namespace exaStamp
{
  using namespace exanb;

  namespace
  {
    inline bool mat3_close( const Mat3d& a, const Mat3d& b, double eps )
    {
      return std::abs(a.m11-b.m11)<=eps && std::abs(a.m12-b.m12)<=eps && std::abs(a.m13-b.m13)<=eps
          && std::abs(a.m21-b.m21)<=eps && std::abs(a.m22-b.m22)<=eps && std::abs(a.m23-b.m23)<=eps
          && std::abs(a.m31-b.m31)<=eps && std::abs(a.m32-b.m32)<=eps && std::abs(a.m33-b.m33)<=eps;
    }
  }

  int ClusterGraph::create_cluster( int structure_type )
  {
    const int id = static_cast<int>( clusters.size() );
    clusters.push_back( Cluster{ id, structure_type, 0, onika::math::make_identity_matrix(), {} } );
    return id;
  }

  int ClusterGraph::self_transition( int cluster_id )
  {
    Cluster& c = clusters[cluster_id];
    for( int ti : c.transition_indices )
    {
      if( transitions[ti].is_self_transition() ) { return ti; }
    }
    const int idx = static_cast<int>( transitions.size() );
    transitions.push_back( ClusterTransition{ cluster_id, cluster_id, onika::math::make_identity_matrix(), idx, 0, 0 } );
    c.transition_indices.push_back( idx );
    return idx;
  }

  int ClusterGraph::create_transition( int cluster_a, int cluster_b, const Mat3d& tm, int distance, double tm_epsilon )
  {
    if( cluster_a == cluster_b && mat3_close( tm, onika::math::make_identity_matrix(), tm_epsilon ) )
    {
      return self_transition( cluster_a );
    }

    Cluster& ca = clusters[cluster_a];
    for( int ti : ca.transition_indices )
    {
      if( transitions[ti].cluster2 == cluster_b && mat3_close( transitions[ti].tm, tm, tm_epsilon ) ) { return ti; }
    }

    const int idx_ab = static_cast<int>( transitions.size() );
    const int idx_ba = idx_ab + 1;
    transitions.push_back( ClusterTransition{ cluster_a, cluster_b, tm, idx_ba, distance, 0 } );
    transitions.push_back( ClusterTransition{ cluster_b, cluster_a, inverse(tm), idx_ab, distance, 0 } );
    clusters[cluster_a].transition_indices.push_back( idx_ab );
    clusters[cluster_b].transition_indices.push_back( idx_ba );
    return idx_ab;
  }

  int ClusterGraph::determine_transition( int cluster_a, int cluster_b )
  {
    if( cluster_a == cluster_b ) { return self_transition( cluster_a ); }

    const Cluster& ca = clusters[cluster_a];
    for( int ti : ca.transition_indices )
    {
      if( transitions[ti].cluster2 == cluster_b ) { return ti; }
    }

    // Hardcoded shortest-path search up to distance 2, exactly matching OVITO's own
    // ClusterGraph::determineClusterTransition (its _maximumClusterDistance is always 2).
    int shortest_distance = 3;
    int best_t1 = -1, best_t2 = -1;
    for( int t1 : ca.transition_indices )
    {
      const ClusterTransition& tr1 = transitions[t1];
      if( tr1.cluster2 == cluster_a ) { continue; } // skip self-transition
      const Cluster& mid = clusters[tr1.cluster2];
      for( int t2 : mid.transition_indices )
      {
        const ClusterTransition& tr2 = transitions[t2];
        if( tr2.cluster2 == cluster_b )
        {
          const int d = tr1.distance + tr2.distance;
          if( d < shortest_distance ) { shortest_distance = d; best_t1 = t1; best_t2 = t2; }
          break;
        }
      }
    }

    if( best_t1 < 0 ) { return -1; }
    const Mat3d combined = transitions[best_t2].tm * transitions[best_t1].tm;
    return create_transition( cluster_a, cluster_b, combined, shortest_distance );
  }
}
