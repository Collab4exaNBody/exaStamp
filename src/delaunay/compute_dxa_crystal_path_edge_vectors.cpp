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

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/math/basic_types.h>

#include <exaStamp/delaunay/delaunay_tessellation.h>
#include <exaStamp/delaunay/dxa_crystal_path.h>

#include <algorithm>

// DXA elastic-mapping step (iii), port of ElasticMapping::assignIdealVectorsToEdges()
// (ovito/src/ovito/crystalanalysis/modifier/dxa/ElasticMapping.cpp): assigns each Delaunay
// tessellation edge its ideal vector via CrystalPathFinder (dxa_crystal_path_find), then
// re-expresses it specifically in the lower-indexed vertex's own cluster frame -- the edge also
// records the cluster transition from v0's cluster to v1's cluster, needed by the per-tetrahedron
// elastic-mapping compatibility test (stage 3, not yet ported) for its own Frank-rotation check.
//
// Not yet wired into the rest of the pipeline (compute_dxa_tet_classification/
// compute_interface_mesh/compute_dxa_circuit_sweep still consume the older compute_dxa_edge_vectors
// path) -- see src/delaunay/README.md "Remaining work" item 0.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXACrystalPathEdgeVectors : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation      , delaunay_tessellation      , INPUT , REQUIRED );
    ADD_SLOT( DXALatticeCorrespondence  , dxa_lattice_correspondence , INPUT , REQUIRED );
    ADD_SLOT( DXALatticeClusters        , dxa_lattice_clusters       , INPUT_OUTPUT , REQUIRED ); // graph may gain cached transitions
    ADD_SLOT( DXACrystalPathEdgeVectors , dxa_crystal_path_edge_vectors , OUTPUT );
    ADD_SLOT( long                      , crystal_path_steps         , INPUT , 4 , DocString{"Max path length (hops) CrystalPathFinder is allowed to search through the good-crystal region -- OVITO's own default"} );
    ADD_SLOT( long                      , n_edges_resolved           , OUTPUT , DocString{"Number of tessellation edges resolved (out of dxa_crystal_path_edge_vectors->edges.size())"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXALatticeCorrespondence& lc = *dxa_lattice_correspondence;
      DXALatticeClusters& clusters = *dxa_lattice_clusters;
      DXACrystalPathEdgeVectors& result = *dxa_crystal_path_edge_vectors;

      result.edges.clear();
      result.ideal_vector.clear();
      result.cluster_transition.clear();
      result.resolved.clear();
      result.edge_index.clear();

      // deduplicate edges: every tetrahedron contributes its 6 vertex pairs, sort+unique -- same
      // technique as compute_dxa_edge_vectors.cpp.
      static constexpr int edge_lv[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
      std::vector<std::array<uint32_t,2>> all_edges;
      all_edges.reserve( mesh.tetrahedra.size() * 6 );
      for(const auto& tet : mesh.tetrahedra)
      {
        for(int e=0;e<6;e++)
        {
          uint32_t a = tet[ edge_lv[e][0] ];
          uint32_t b = tet[ edge_lv[e][1] ];
          all_edges.push_back( a<b ? std::array<uint32_t,2>{a,b} : std::array<uint32_t,2>{b,a} );
        }
      }
      std::sort( all_edges.begin(), all_edges.end() );
      all_edges.erase( std::unique( all_edges.begin(), all_edges.end() ), all_edges.end() );

      result.edges = std::move( all_edges );
      result.ideal_vector.resize( result.edges.size(), Vec3d{0,0,0} );
      result.cluster_transition.assign( result.edges.size(), -1 );
      result.resolved.assign( result.edges.size(), false );

      // ElasticMapping::assignVerticesToClusters() port -- propagates a cluster to EVERY
      // tessellation vertex (not just classified atoms) via ordinary tessellation-edge adjacency,
      // so a genuinely unclassified vertex still has a meaningful frame to gate on / express an
      // edge's ideal vector in. See dxa_propagate_vertex_clusters()'s own header comment: this is
      // used only for the gate check and re-expression target below, never inside
      // dxa_crystal_path_find() itself (which needs the RAW, un-propagated cluster for its own
      // reverse-neighbor-search branch).
      const std::vector<int32_t> vertex_cluster = dxa_propagate_vertex_clusters( result.edges, mesh.vertex_particle_index, clusters );

      const int max_path_length = static_cast<int>( *crystal_path_steps );
      std::vector<uint8_t> visited( lc.structure_type.size(), 0 );

      size_t n_resolved = 0;
      long n_unresolved_no_cluster = 0, n_unresolved_no_path = 0, n_unresolved_no_reexpress = 0, n_unresolved_no_edge_transition = 0;
      for(size_t e=0;e<result.edges.size();e++)
      {
        const uint32_t v0 = result.edges[e][0];
        const uint32_t v1 = result.edges[e][1];
        result.edge_index[ DXACrystalPathEdgeVectors::key(v0,v1) ] = static_cast<uint32_t>(e);

        const int64_t a0 = mesh.vertex_particle_index[v0];
        const int64_t a1 = mesh.vertex_particle_index[v1];

        const int32_t cluster0 = vertex_cluster[v0];
        const int32_t cluster1 = vertex_cluster[v1];
        if( cluster0 == 0 || cluster1 == 0 ) { ++n_unresolved_no_cluster; continue; }

        auto idealVector = dxa_crystal_path_find( lc, clusters, a0, a1, max_path_length, visited );
        if( !idealVector ) { ++n_unresolved_no_path; continue; }

        // Translate to vertex v0's own cluster frame.
        Vec3d localVec;
        if( idealVector->cluster == cluster0 )
        {
          localVec = idealVector->vec;
        }
        else
        {
          const int t = clusters.graph.determine_transition( idealVector->cluster, cluster0 );
          if( t < 0 ) { ++n_unresolved_no_reexpress; continue; }
          localVec = clusters.graph.transitions[t].transform( idealVector->vec );
        }

        // The edge's own transition, v0's cluster -> v1's cluster (may be disconnected in a
        // pathological multi-grain case even though a path vector was found through some other
        // route -- guard the same way OVITO does).
        const int edge_transition = clusters.graph.determine_transition( cluster0, cluster1 );
        if( edge_transition < 0 ) { ++n_unresolved_no_edge_transition; continue; }

        result.ideal_vector[e] = localVec;
        result.cluster_transition[e] = edge_transition;
        result.resolved[e] = true;
        ++n_resolved;
      }

      *n_edges_resolved = static_cast<long>( n_resolved );
      lout << "compute_dxa_crystal_path_edge_vectors: " << n_resolved << " / " << result.edges.size()
           << " edges resolved (unresolved breakdown: " << n_unresolved_no_cluster << " endpoint unclassified, "
           << n_unresolved_no_path << " no path found, " << n_unresolved_no_reexpress << " path found but couldn't re-express in v0's cluster, "
           << n_unresolved_no_edge_transition << " no direct cluster-to-cluster transition)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA elastic-mapping step (iii), the real mechanism: port of OVITO's ElasticMapping::
assignIdealVectorsToEdges() + CrystalPathFinder, assigning each Delaunay tessellation edge its ideal
vector via a bounded graph walk through the good-crystal region (compute_dxa_lattice_correspondence/
compute_dxa_lattice_clusters), not a per-endpoint angle-snap. See this file's own header comment and
src/delaunay/README.md "Remaining work" item 0.

Usage example:

compute_delaunay: {}
compute_dxa_lattice_correspondence: { rcut: 6.0 ang }
compute_dxa_lattice_clusters: {}
compute_dxa_crystal_path_edge_vectors: { crystal_path_steps: 4 }

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_dxa_crystal_path_edge_vectors)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_crystal_path_edge_vectors", make_simple_operator< ComputeDXACrystalPathEdgeVectors > );
  }

}
