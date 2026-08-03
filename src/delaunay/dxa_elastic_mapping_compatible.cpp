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

#include <exaStamp/delaunay/dxa_crystal_path.h>

#include <cmath>

// Port of ElasticMapping::isElasticMappingCompatible() (ovito/src/ovito/crystalanalysis/modifier/
// dxa/ElasticMapping.cpp) -- see dxa_crystal_path.h for the full rationale. Deliberately grid-
// independent (only needs DXACrystalPathEdgeVectors/DXALatticeClusters + a tet's 4 global vertex
// indices), same pattern as dxa_lattice_clusters_algo.cpp/dxa_crystal_path_finder.cpp.
namespace exaStamp
{
  using namespace exanb;

  namespace
  {
    // Same epsilons OVITO uses (Cluster.h/ClusterVector.h): CA_LATTICE_VECTOR_EPSILON=1e-3 for a
    // Burgers-vector-closure zero test, CA_TRANSITION_MATRIX_EPSILON=1e-4 for a transition-matrix
    // identity test.
    constexpr double CA_LATTICE_VECTOR_EPSILON = 1e-3;
    constexpr double CA_TRANSITION_MATRIX_EPSILON = 1e-4;

    inline bool vec3_is_zero( const Vec3d& v, double eps ) { return std::abs(v.x)<=eps && std::abs(v.y)<=eps && std::abs(v.z)<=eps; }

    inline bool mat3_close( const Mat3d& a, const Mat3d& b, double eps )
    {
      return std::abs(a.m11-b.m11)<=eps && std::abs(a.m12-b.m12)<=eps && std::abs(a.m13-b.m13)<=eps
          && std::abs(a.m21-b.m21)<=eps && std::abs(a.m22-b.m22)<=eps && std::abs(a.m23-b.m23)<=eps
          && std::abs(a.m31-b.m31)<=eps && std::abs(a.m32-b.m32)<=eps && std::abs(a.m33-b.m33)<=eps;
    }

    // Directed edge vector + its own cluster transition, in the SAME sense OVITO's
    // ElasticMapping::getEdgeClusterVector() resolves an edge stored canonically as (v0<v1) when
    // queried in either direction.
    struct DirectedEdge { Vec3d vec; int transition; }; // transition: index into ClusterGraph::transitions, vertex1->vertex2 of the QUERY direction

    bool get_directed_edge( const DXACrystalPathEdgeVectors& ev, DXALatticeClusters& clusters, uint32_t vertex1, uint32_t vertex2, DirectedEdge& out )
    {
      const auto it = ev.edge_index.find( DXACrystalPathEdgeVectors::key(vertex1,vertex2) );
      if( it == ev.edge_index.end() ) { return false; }
      const uint32_t e = it->second;
      if( !ev.resolved[e] ) { return false; }

      if( ev.edges[e][0] == vertex1 )
      {
        // Query direction matches storage direction (v0->v1) exactly.
        out.vec = ev.ideal_vector[e];
        out.transition = ev.cluster_transition[e];
      }
      else
      {
        // Query direction is reversed (v1->v0): negate the vector, transport it through the
        // stored transition into v1's own frame, and use the REVERSE transition.
        const auto& t = clusters.graph.transitions[ ev.cluster_transition[e] ];
        out.vec = t.transform( Vec3d{ -ev.ideal_vector[e].x, -ev.ideal_vector[e].y, -ev.ideal_vector[e].z } );
        out.transition = t.reverse_index;
      }
      return true;
    }
  }

  bool dxa_is_elastic_mapping_compatible(
      const DXACrystalPathEdgeVectors& edge_vectors, DXALatticeClusters& clusters,
      const std::array<uint32_t,4>& tet_vertices )
  {
    static constexpr int edge_lv[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };

    DirectedEdge edges[6];
    for(int i=0;i<6;i++)
    {
      const uint32_t v1 = tet_vertices[ edge_lv[i][0] ];
      const uint32_t v2 = tet_vertices[ edge_lv[i][1] ];
      if( !get_directed_edge( edge_vectors, clusters, v1, v2, edges[i] ) ) { return false; }
    }

    // Burgers-circuit closure test on each of the tet's 4 triangular faces -- the 3 directed edges
    // going around a face (0->1->2->0-equivalent, expressed via the 6 tet edges) must sum to zero
    // once transported into a common frame. Faces expressed as (edgeA, edgeB, edgeC) with
    // edgeC == edgeA followed by edgeB, exactly OVITO's own `circuits` table.
    static constexpr int circuits[4][3] = { {0,4,2}, {1,5,2}, {0,3,1}, {3,5,4} };

    for(int face=0; face<4; face++)
    {
      const DirectedEdge& e0 = edges[ circuits[face][0] ];
      const DirectedEdge& e1 = edges[ circuits[face][1] ];
      const DirectedEdge& e2 = edges[ circuits[face][2] ];

      const auto& t0 = clusters.graph.transitions[ e0.transition ];
      const auto& t0rev = clusters.graph.transitions[ t0.reverse_index ];
      Vec3d burgers = e0.vec + t0rev.transform( e1.vec ) - e2.vec;
      if( !vec3_is_zero( burgers, CA_LATTICE_VECTOR_EPSILON ) ) { return false; }
    }

    // Frank-rotation (disclination) test on the same 4 faces: the product of the 3 edges' own
    // cluster transitions must be the identity -- catches a self-consistent-looking Burgers
    // closure that's only achieved by an unphysical net rotation around the loop (a genuine
    // grain-boundary/disclination signal, or here, the general "this circuit doesn't actually
    // encircle a coherent patch of crystal" guard).
    for(int face=0; face<4; face++)
    {
      const auto& t0 = clusters.graph.transitions[ edges[ circuits[face][0] ].transition ];
      const auto& t1 = clusters.graph.transitions[ edges[ circuits[face][1] ].transition ];
      const auto& t2 = clusters.graph.transitions[ edges[ circuits[face][2] ].transition ];
      if( !t0.is_self_transition() || !t1.is_self_transition() || !t2.is_self_transition() )
      {
        const auto& t2rev = clusters.graph.transitions[ t2.reverse_index ];
        const Mat3d frankRotation = t2rev.tm * t1.tm * t0.tm;
        if( !mat3_close( frankRotation, onika::math::make_identity_matrix(), CA_TRANSITION_MATRIX_EPSILON ) ) { return false; }
      }
    }

    return true;
  }
}
