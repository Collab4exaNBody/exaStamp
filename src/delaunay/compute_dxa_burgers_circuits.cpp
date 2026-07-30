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
#include <exaStamp/delaunay/dxa_edge_vectors.h>

#include <algorithm>
#include <array>
#include <deque>
#include <unordered_map>
#include <vector>

// DXA pipeline steps (vi)-(vii): Burgers circuit detection via a global elastic mapping.
//
// First attempt (superseded, see git history) traced a small circuit -- the ring of tetrahedra
// sharing one tessellation edge -- around each edge touching a "bad" tetrahedron. That doesn't
// work: on a real dislocation quadrupole test case, 93% of such rings included an atom PTM
// couldn't classify at all, because a genuine dislocation CORE is exactly where the crystal is
// too disordered for PTM to match anything -- and any loop that actually encircles a line defect
// is topologically forced to pass near its core. A small local circuit can never see the signal;
// it has to go around the defect at a safe distance, through the surrounding good lattice.
//
// Fix: build a global elastic mapping. Pick any vertex reachable via RESOLVED edges (step iii) as
// a root, and propagate an accumulated "ideal lattice position" (ideal_pos) to every vertex
// reachable from it by a spanning tree of resolved edges -- well-defined because a tree has no
// cycles. Then for every resolved edge (u,v) in the WHOLE mesh (not just tree edges), compare its
// own direct ideal vector against ideal_pos[v]-ideal_pos[u]: for a tree edge these trivially
// match (that's literally how ideal_pos was assigned). For any other resolved edge, a mismatch
// means the loop [root->u (tree) -> v (direct edge) -> root (tree, reversed)] has non-zero net
// circulation -- exactly the Burgers vector of whatever dislocation that loop encircles. Since the
// tree path can route arbitrarily far around the disordered core (using whatever resolved edges
// happen to be available anywhere in the good lattice), it never needs to touch the unresolved
// atoms right at the core -- unlike the small-ring approach, this one can actually see real
// dislocations.
//
// This is also exactly the mechanism for filtering "bad candidates": an isolated misclassified
// tetrahedron (thermal noise, not a real defect) simply has no non-tree resolved edge nearby with
// a non-zero residual -- min_burgers_norm separates genuine signal from noise. Any "bad"
// tetrahedron not connected (through bad-tet adjacency) to a tetrahedron touching a confirmed
// signal edge gets reclassified back to "good" in DXATetClassification, which is exactly what
// compute_interface_mesh needs to stop including it in the interface surface.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXABurgersCircuits : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation , INPUT , REQUIRED );
    ADD_SLOT( DXAEdgeVectors       , dxa_edge_vectors       , INPUT , REQUIRED );
    ADD_SLOT( DXATetClassification , dxa_tet_classification , INPUT_OUTPUT , DocString{"Step (iv) classification, refined in place: a 'bad' tetrahedron not connected to any confirmed dislocation core edge is reclassified 'good' (1.0)."} );
    ADD_SLOT( DXABurgersCircuits   , dxa_burgers_circuits   , OUTPUT );
    ADD_SLOT( double               , min_burgers_norm       , INPUT , 0.5 , DocString{"Minimum Burgers circuit closure-failure magnitude (same units as DXAEdgeVectors::ideal_vector, i.e. PTM's template units) to count as a genuine dislocation signal rather than numerical/classification noise."} );
    ADD_SLOT( long                 , n_core_edges           , OUTPUT , DocString{"Number of confirmed dislocation signal edges found (non-tree resolved edges whose direct ideal vector disagrees with the spanning tree's accumulated position)"} );
    ADD_SLOT( long                 , n_tets_reclassified     , OUTPUT , DocString{"Number of 'bad' tetrahedra reclassified 'good' (filtered out as noise, not connected to any confirmed signal edge)"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXAEdgeVectors& ev = *dxa_edge_vectors;
      DXATetClassification& cls = *dxa_tet_classification;
      const size_t n_tets = mesh.tetrahedra.size();
      const size_t n_verts = mesh.vertices.size();

      // vertex -> tets containing it: used both to find which bad tets a signal edge's endpoints
      // touch, and (below) for the flood-fill's bad-tet adjacency
      std::vector<std::vector<uint32_t>> vertex_tets( n_verts );
      for(uint32_t t=0; t<n_tets; t++)
      {
        const auto& tet = mesh.tetrahedra[t];
        for(int i=0;i<4;i++) { vertex_tets[ tet[i] ].push_back(t); }
      }

      // resolved-edge adjacency graph over vertices, for the spanning-tree BFS
      std::vector<std::vector<uint32_t>> resolved_adj( n_verts );
      for(size_t i=0;i<ev.edges.size();i++)
      {
        if( ! ev.resolved[i] ) { continue; }
        const uint32_t a = ev.edges[i][0], b = ev.edges[i][1];
        resolved_adj[a].push_back(b);
        resolved_adj[b].push_back(a);
      }

      // BFS spanning forest: one tree per connected component of the resolved-edge graph,
      // propagating an accumulated ideal-lattice position from an arbitrary root
      std::vector<Vec3d> ideal_pos( n_verts, Vec3d{0.,0.,0.} );
      std::vector<int8_t> visited( n_verts, 0 );
      std::deque<uint32_t> bfs_queue;
      for(uint32_t root=0; root<n_verts; root++)
      {
        if( visited[root] || resolved_adj[root].empty() ) { continue; }
        visited[root] = 1;
        ideal_pos[root] = Vec3d{0.,0.,0.};
        bfs_queue.push_back(root);
        while( ! bfs_queue.empty() )
        {
          const uint32_t u = bfs_queue.front(); bfs_queue.pop_front();
          for(uint32_t v : resolved_adj[u])
          {
            if( visited[v] ) { continue; }
            visited[v] = 1;
            const auto it = ev.edge_index.find( DXAEdgeVectors::key(u,v) );
            const Vec3d step = ( u < v ) ? ev.ideal_vector[it->second] : -ev.ideal_vector[it->second];
            ideal_pos[v] = ideal_pos[u] + step;
            bfs_queue.push_back(v);
          }
        }
      }

      // every resolved edge, tree or not: compare its own direct ideal vector against the tree's
      // accumulated position difference. Tree edges trivially match (that's how ideal_pos got
      // assigned); any mismatch elsewhere is a genuine circulation signal.
      DXABurgersCircuits& result = *dxa_burgers_circuits;
      result.edges.clear();
      result.burgers_vector.clear();
      const double min_norm = *min_burgers_norm;

      std::vector<uint8_t> confirmed( n_tets, 0 );
      std::deque<uint32_t> to_visit;

      for(size_t i=0;i<ev.edges.size();i++)
      {
        if( ! ev.resolved[i] ) { continue; }
        const uint32_t a = ev.edges[i][0], b = ev.edges[i][1]; // a<b
        const Vec3d direct = ev.ideal_vector[i]; // a->b
        const Vec3d residual = direct - ( ideal_pos[b] - ideal_pos[a] );
        if( norm(residual) <= min_norm ) { continue; }

        result.edges.push_back( { a, b } );
        result.burgers_vector.push_back( residual );

        for(uint32_t t : vertex_tets[a]) { if( cls.good[t]==0.0 && ! confirmed[t] ) { confirmed[t]=1; to_visit.push_back(t); } }
        for(uint32_t t : vertex_tets[b]) { if( cls.good[t]==0.0 && ! confirmed[t] ) { confirmed[t]=1; to_visit.push_back(t); } }
      }

      // flood-fill: grow each signal edge's seed tet(s) through bad-tet adjacency to the whole
      // connected bad-tet cluster -- a real dislocation's entire tube counts as confirmed, not
      // just the specific tets immediately touching a signal edge. Deliberately VERTEX-adjacency
      // here (any two bad tets sharing at least one vertex), not the stricter face-adjacency
      // compute_interface_mesh itself needs: where a dislocation's cross-section pinches down to
      // a thin/irregular sliver, consecutive bad tets along the tube can share only an edge or a
      // vertex, not a full face, in a irregular Delaunay tessellation -- face-only adjacency was
      // failing to flood-fill through those pinches, incorrectly reclassifying genuine core tets
      // back to "good" and splitting one continuous dislocation into disconnected mesh islands
      // (each individually closed, but the tube visually broken into segments).
      std::vector<std::vector<uint32_t>> bad_adjacency( n_tets );
      for(uint32_t v=0; v<n_verts; v++)
      {
        std::vector<uint32_t> bad_here;
        for(uint32_t t : vertex_tets[v]) { if( cls.good[t]==0.0 ) { bad_here.push_back(t); } }
        for(size_t i=0;i<bad_here.size();i++)
        {
          for(size_t j=i+1;j<bad_here.size();j++)
          {
            bad_adjacency[bad_here[i]].push_back(bad_here[j]);
            bad_adjacency[bad_here[j]].push_back(bad_here[i]);
          }
        }
      }
      while( ! to_visit.empty() )
      {
        const uint32_t t = to_visit.front(); to_visit.pop_front();
        for(uint32_t nb : bad_adjacency[t]) { if( ! confirmed[nb] ) { confirmed[nb] = 1; to_visit.push_back(nb); } }
      }

      size_t n_reclassified = 0;
      for(uint32_t t=0; t<n_tets; t++)
      {
        if( cls.good[t]==0.0 && ! confirmed[t] ) { cls.good[t] = 1.0; ++n_reclassified; }
      }

      *n_core_edges = static_cast<long>( result.edges.size() );
      *n_tets_reclassified = static_cast<long>( n_reclassified );
      lout << "compute_dxa_burgers_circuits: " << result.edges.size() << " confirmed signal edges, "
           << n_reclassified << " bad tetrahedra reclassified good (filtered as noise)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline steps (vi)-(vii): detects genuine dislocation signal via a global elastic mapping.
Builds a spanning tree over all vertices connected by resolved edges (compute_dxa_edge_vectors),
propagating an accumulated ideal-lattice position from an arbitrary root. For every resolved edge
in the mesh, compares its own direct ideal vector against the tree's accumulated position
difference -- a mismatch above min_burgers_norm is a genuine Burgers-vector signal (the tree path
can route arbitrarily far around a disordered dislocation core, so this sees real defects that a
small local circuit can't). Any "bad" tetrahedron (compute_dxa_tet_classification) not connected
to a confirmed signal edge is reclassified "good" in place -- this is how isolated misclassified
tetrahedra (not real defects) get filtered out before compute_interface_mesh builds the final
interface surface.

Usage example:

compute_delaunay: {}
compute_dxa_edge_vectors: { target_structure: BCC }
compute_dxa_tet_classification: {}
compute_dxa_burgers_circuits: { min_burgers_norm: 0.5 }
compute_interface_mesh: {}

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_dxa_burgers_circuits)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_burgers_circuits", make_simple_operator< ComputeDXABurgersCircuits > );
  }

}
