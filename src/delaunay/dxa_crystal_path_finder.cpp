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

#include <deque>

// Port of CrystalPathFinder::findPath() (ovito/src/ovito/crystalanalysis/modifier/dxa/
// CrystalPathFinder.cpp) -- see dxa_crystal_path.h for the full rationale. Deliberately grid-
// independent (only needs DXALatticeCorrespondence/DXALatticeClusters) so it's unit-testable
// standalone, same pattern as dxa_lattice_clusters_algo.cpp.
namespace exaStamp
{
  using namespace exanb;

  namespace
  {
    inline int find_neighbor_slot( const DXALatticeCorrespondence& lc, DXALatticeClusters& clusters, int64_t atom, int64_t target )
    {
      const int count = clusters.atom_neighbor_count[atom];
      for(int j=0;j<count;j++) { if( lc.neighbor(atom,j) == target ) { return j; } }
      return -1;
    }

    inline Vec3d neg( const Vec3d& v ) { return Vec3d{ -v.x, -v.y, -v.z }; }
  }

  std::optional<ClusterVector> dxa_crystal_path_find(
      const DXALatticeCorrespondence& lc, DXALatticeClusters& clusters,
      int64_t atom1, int64_t atom2, int max_path_length,
      std::vector<uint8_t>& visited_scratch )
  {
    const int32_t cluster1 = clusters.atom_cluster[atom1];
    const int32_t cluster2 = clusters.atom_cluster[atom2];

    // Direct-neighbor shortcut, both directions.
    if( cluster1 != 0 )
    {
      const int j = find_neighbor_slot( lc, clusters, atom1, atom2 );
      if( j >= 0 )
      {
        const LatticeStructure& ls = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[atom1] ) );
        return ClusterVector{ clusters.ideal_vector(atom1,j,ls), cluster1 };
      }
    }
    else if( cluster2 != 0 )
    {
      const int j = find_neighbor_slot( lc, clusters, atom2, atom1 );
      if( j >= 0 )
      {
        const LatticeStructure& ls = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[atom2] ) );
        return ClusterVector{ neg( clusters.ideal_vector(atom2,j,ls) ), cluster2 };
      }
    }

    if( max_path_length <= 1 ) { return std::nullopt; }

    struct PathNode { int64_t atom; ClusterVector vec; int distance; };
    std::vector<PathNode> nodes;
    nodes.reserve(64);
    nodes.push_back( PathNode{ atom1, ClusterVector{}, 0 } );
    visited_scratch[atom1] = 1;

    std::optional<ClusterVector> result;
    for(size_t qi=0; qi<nodes.size() && !result; qi++)
    {
      const int64_t cur = nodes[qi].atom;
      const int32_t curCluster = clusters.atom_cluster[cur];
      const int count = clusters.atom_neighbor_count[cur];

      for(int ni=0; ni<count; ni++)
      {
        const int64_t neighbor = lc.neighbor( cur, ni );
        if( neighbor < 0 ) { continue; }
        if( visited_scratch[neighbor] ) { continue; }
        if( nodes[qi].distance >= max_path_length-1 && neighbor != atom2 ) { continue; }

        ClusterVector step;
        if( curCluster != 0 )
        {
          const LatticeStructure& ls = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[cur] ) );
          step = ClusterVector{ clusters.ideal_vector(cur,ni,ls), curCluster };
        }
        else
        {
          // Reverse neighbor search: cur is a disordered stepping-stone atom appended into by a
          // classified neighbor during connectClusters -- look up how THAT neighbor sees cur.
          const int32_t neighborCluster = clusters.atom_cluster[neighbor];
          if( neighborCluster == 0 ) { continue; }
          const int j = find_neighbor_slot( lc, clusters, neighbor, cur );
          if( j < 0 ) { continue; }
          const LatticeStructure& ls2 = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[neighbor] ) );
          step = ClusterVector{ neg( clusters.ideal_vector(neighbor,j,ls2) ), neighborCluster };
        }

        ClusterVector pathVector = nodes[qi].vec;
        if( pathVector.cluster == step.cluster )
        {
          pathVector.vec = pathVector.vec + step.vec;
        }
        else if( pathVector.cluster != -1 )
        {
          const int ti = clusters.graph.determine_transition( step.cluster, pathVector.cluster );
          if( ti < 0 ) { continue; } // disconnected components -- can't concatenate, try another neighbor
          pathVector.vec = pathVector.vec + clusters.graph.transitions[ti].transform( step.vec );
        }
        else
        {
          pathVector = step;
        }

        if( neighbor == atom2 )
        {
          result = pathVector;
          break;
        }

        if( nodes[qi].distance < max_path_length-1 )
        {
          nodes.push_back( PathNode{ neighbor, pathVector, nodes[qi].distance+1 } );
          visited_scratch[neighbor] = 1;
        }
      }
    }

    for( const auto& node : nodes ) { visited_scratch[node.atom] = 0; }

    return result;
  }

  std::vector<int32_t> dxa_propagate_vertex_clusters(
      const std::vector<std::array<uint32_t,2>>& tessellation_edges,
      const std::vector<uint32_t>& vertex_particle_index,
      const DXALatticeClusters& clusters )
  {
    const size_t n_vertices = vertex_particle_index.size();
    std::vector<int32_t> vertex_cluster( n_vertices );
    for(size_t v=0; v<n_vertices; v++) { vertex_cluster[v] = clusters.atom_cluster[ vertex_particle_index[v] ]; }

    std::vector<std::vector<uint32_t>> adj( n_vertices );
    for(const auto& e : tessellation_edges) { adj[e[0]].push_back(e[1]); adj[e[1]].push_back(e[0]); }

    std::deque<uint32_t> queue;
    std::vector<uint8_t> queued( n_vertices, 0 );
    for(size_t v=0; v<n_vertices; v++) { if( vertex_cluster[v] != 0 ) { queue.push_back( static_cast<uint32_t>(v) ); queued[v] = 1; } }

    while( !queue.empty() )
    {
      const uint32_t v = queue.front(); queue.pop_front();
      for( uint32_t w : adj[v] )
      {
        if( vertex_cluster[w] == 0 )
        {
          vertex_cluster[w] = vertex_cluster[v];
          if( !queued[w] ) { queue.push_back(w); queued[w] = 1; }
        }
      }
    }

    return vertex_cluster;
  }
}
