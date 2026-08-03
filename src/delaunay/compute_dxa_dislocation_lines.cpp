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
#include <exaStamp/delaunay/dxa_dislocation_lines.h>

#include <algorithm>
#include <array>
#include <functional>
#include <set>
#include <string>
#include <vector>

// DXA pipeline steps (viii)-(ix): extract 1D dislocation line geometry from the disordered
// ("hole") atom population, and detect junctions where multiple lines meet. Purely a
// DXAEdgeVectors + DXABurgersCircuits post-process (no grid access needed), same "plain operator"
// pattern as compute_dxa_tet_classification.
//
// Approach (deliberately atom-graph-based, not a surface-mesh sweep -- see rationale below):
//   1. Build the "hole graph": adjacency among non-crystalline (vertex_matches_target==0) atoms,
//      from DXAEdgeVectors' own already-deduplicated edge list.
//   2. Topological thinning: the hole-atom population is NOT a 1D chain -- it's a genuine tube
//      SURFACE (atoms wrap around the tube's own circumference too, not just along its length), so
//      a naive fixed-depth erosion (this operator's first version) either does nothing (radius ~1
//      almost everywhere -> erosion round 1 removes nearly everything) or, with zero erosion,
//      treats nearly every atom as a "junction" (high degree from the tube's own circumferential
//      neighbors, not a real dislocation branch point) -- confirmed on the real quadrupole test
//      case (measured both failure modes, see src/delaunay/README.md). Fixed with real graph
//      thinning instead: repeatedly find any atom with degree > 2 in the *current* (shrinking)
//      graph whose removal would NOT disconnect its neighbors (i.e. not an articulation point,
//      Tarjan's algorithm) and remove it; converges to a graph where every remaining atom either
//      has degree <= 2 (an ordinary curve point) or is a genuine articulation point (a real
//      junction that can't be thinned away without breaking connectivity) -- this naturally adapts
//      to locally-varying tube radius, unlike a single fixed depth.
//   3. Trace the thinned graph: junctions/endpoints are vertices with degree != 2; walk from each
//      along its own incident edges until hitting another such vertex (an open line), or -- for
//      any edges left over, forming a pure cycle of degree-2 vertices with no junction at all --
//      walk all the way around back to the start (a closed loop).
//   4. Per-line Burgers vector: average the burgers_vector of every DXABurgersCircuits "confirmed
//      signal edge" (crystalline-crystalline, from the elastic-mapping spanning tree) whose
//      midpoint's nearest skeleton vertex belongs to that line.
//
// Why atom-graph-based rather than a surface-mesh "sweep" (DXA1.3.6's DXATracing.cpp,
// burgersSearchWalkEdge -- the real reference algorithm's own approach, sweeping an elementary
// Burgers circuit stepwise around the interface mesh's tube): the atomistic interface mesh's own
// vertices already ARE the disordered atoms (see compute_atomistic_interface_mesh.cpp), so their
// own adjacency graph already encodes almost all the topology a mesh-sweep would have to
// rediscover, at a fraction of the implementation cost (no half-edge walk, no per-step elastic-
// mapping re-evaluation) -- traded off against a less rigorous notion of "skeleton" than a true
// medial axis (graph thinning approximates it well for tube-like topology, but isn't one).
// ponytail: no minimum-segment-length pruning -- a short noise "hair" sticking off a real junction
// (e.g. from classification noise) shows up as its own short line rather than being merged away;
// add a length filter if real data shows this is common.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXADislocationLines : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation  , INPUT , REQUIRED );
    ADD_SLOT( DXAEdgeVectors       , dxa_edge_vectors       , INPUT , REQUIRED );
    ADD_SLOT( DXABurgersCircuits   , dxa_burgers_circuits   , INPUT , REQUIRED );
    ADD_SLOT( DXADislocationLines  , dxa_dislocation_lines  , OUTPUT );
    ADD_SLOT( long                 , n_lines_extracted      , OUTPUT , DocString{"Number of dislocation lines extracted (open segments + closed loops)"} );
    ADD_SLOT( long                 , n_junctions            , OUTPUT , DocString{"Number of skeleton atoms where 3+ lines meet"} );
    ADD_SLOT( long                 , n_thinning_removed     , OUTPUT , DocString{"Number of hole atoms removed by topological thinning (started as surface/circumference atoms, ended up not part of the 1D skeleton)"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXAEdgeVectors& ev = *dxa_edge_vectors;
      const DXABurgersCircuits& bc = *dxa_burgers_circuits;
      const size_t n_vertices = mesh.vertices.size();

      // step 1: hole-graph adjacency
      std::vector<std::vector<uint32_t>> hole_adj( n_vertices );
      long n_hole_atoms = 0;
      for(uint32_t v=0; v<n_vertices; v++) { if( !ev.vertex_matches_target[v] ) { ++n_hole_atoms; } }
      for(const auto& e : ev.edges)
      {
        const uint32_t a = e[0], b = e[1];
        if( !ev.vertex_matches_target[a] && !ev.vertex_matches_target[b] ) { hole_adj[a].push_back(b); hole_adj[b].push_back(a); }
      }

      // step 2: topological thinning via iterative articulation-point-safe removal
      std::vector<uint8_t> alive( n_vertices, 0 );
      for(uint32_t v=0; v<n_vertices; v++) { alive[v] = !ev.vertex_matches_target[v] ? 1 : 0; }

      auto current_degree = [&]( uint32_t v ) -> int
      {
        int d = 0;
        for(uint32_t w : hole_adj[v]) { if( alive[w] ) { ++d; } }
        return d;
      };

      long n_removed = 0;
      bool progress = true;
      while( progress )
      {
        progress = false;

        // Tarjan's articulation points, restricted to the currently-alive subgraph (possibly
        // several disjoint components).
        std::vector<int32_t> disc( n_vertices, -1 ), low( n_vertices, -1 ), parent( n_vertices, -1 );
        std::vector<uint8_t> is_articulation( n_vertices, 0 );
        int32_t timer = 0;

        std::function<void(uint32_t)> dfs = [&]( uint32_t root )
        {
          // explicit stack: (vertex, neighbor-iterator-index) to avoid deep recursion. Indexes
          // into `stack` directly (never holds a reference/structured-binding across a push_back
          // on the same vector -- reallocation would dangle it) for manifest correctness.
          std::vector<std::pair<uint32_t,size_t>> stack;
          disc[root] = low[root] = timer++;
          stack.push_back( { root, 0 } );
          int32_t root_children = 0;

          while( !stack.empty() )
          {
            const size_t top = stack.size() - 1;
            const uint32_t u = stack[top].first;
            const size_t idx = stack[top].second;

            if( idx < hole_adj[u].size() )
            {
              const uint32_t v = hole_adj[u][idx];
              stack[top].second = idx + 1; // completed BEFORE any push_back below -- safe
              if( !alive[v] ) { continue; }
              if( disc[v] < 0 )
              {
                parent[v] = static_cast<int32_t>(u);
                disc[v] = low[v] = timer++;
                if( u == root ) { ++root_children; }
                stack.push_back( { v, 0 } ); // may reallocate -- no reference held across this
              }
              else if( static_cast<int32_t>(v) != parent[u] )
              {
                low[u] = std::min( low[u], disc[v] );
              }
            }
            else
            {
              stack.pop_back();
              if( !stack.empty() )
              {
                const uint32_t p = stack.back().first;
                low[p] = std::min( low[p], low[u] );
                if( static_cast<int32_t>(p) != static_cast<int32_t>(root) && low[u] >= disc[p] ) { is_articulation[p] = 1; }
              }
            }
          }
          if( root_children > 1 ) { is_articulation[root] = 1; }
        };

        for(uint32_t v=0; v<n_vertices; v++) { if( alive[v] && disc[v] < 0 ) { dfs(v); } }

        for(uint32_t v=0; v<n_vertices; v++)
        {
          if( !alive[v] || is_articulation[v] ) { continue; }
          if( current_degree(v) > 2 )
          {
            alive[v] = 0;
            ++n_removed;
            progress = true;
            break; // recompute articulation points fresh after every single removal
          }
        }
      }

      // step 3: trace the thinned graph
      std::vector<std::vector<uint32_t>> skel_adj( n_vertices );
      for(uint32_t v=0; v<n_vertices; v++)
      {
        if( !alive[v] ) { continue; }
        for(uint32_t w : hole_adj[v]) { if( alive[w] ) { skel_adj[v].push_back(w); } }
      }

      DXADislocationLines& result = *dxa_dislocation_lines;
      result.lines.clear();
      result.is_loop.clear();
      result.junction_vertices.clear();

      auto ekey = []( uint32_t a, uint32_t b ) -> std::array<uint32_t,2> { return a<b ? std::array<uint32_t,2>{a,b} : std::array<uint32_t,2>{b,a}; };
      std::set<std::array<uint32_t,2>> visited_edges;

      for(uint32_t v=0; v<n_vertices; v++)
      {
        if( alive[v] && skel_adj[v].size() >= 3 ) { result.junction_vertices.push_back(v); }
      }

      // open segments/dangling ends: walk from every vertex whose skeleton-degree != 2
      for(uint32_t v=0; v<n_vertices; v++)
      {
        if( !alive[v] || skel_adj[v].size() == 2 ) { continue; }
        for(uint32_t w : skel_adj[v])
        {
          const auto k0 = ekey(v,w);
          if( visited_edges.count(k0) ) { continue; }
          visited_edges.insert(k0);

          std::vector<uint32_t> path { v, w };
          uint32_t prev = v, cur = w;
          while( skel_adj[cur].size() == 2 )
          {
            const uint32_t nxt = ( skel_adj[cur][0] == prev ) ? skel_adj[cur][1] : skel_adj[cur][0];
            const auto k = ekey(cur,nxt);
            if( visited_edges.count(k) ) { break; } // guards against a degenerate 2-cycle
            visited_edges.insert(k);
            path.push_back(nxt);
            prev = cur; cur = nxt;
          }
          result.lines.push_back( std::move(path) );
          result.is_loop.push_back(0);
        }
      }

      // closed loops: whatever's left is a pure cycle of degree-2 vertices with no junction at all
      for(uint32_t v=0; v<n_vertices; v++)
      {
        if( !alive[v] || skel_adj[v].size() != 2 ) { continue; }
        for(uint32_t w : skel_adj[v])
        {
          const auto k0 = ekey(v,w);
          if( visited_edges.count(k0) ) { continue; }
          visited_edges.insert(k0);

          std::vector<uint32_t> path { v, w };
          const uint32_t start = v;
          uint32_t prev = v, cur = w;
          while( cur != start )
          {
            const uint32_t nxt = ( skel_adj[cur][0] == prev ) ? skel_adj[cur][1] : skel_adj[cur][0];
            const auto k = ekey(cur,nxt);
            visited_edges.insert(k);
            path.push_back(nxt);
            prev = cur; cur = nxt;
          }
          result.lines.push_back( std::move(path) );
          result.is_loop.push_back(1);
        }
      }

      // step 4: per-line Burgers vector, from nearby DXABurgersCircuits signal edges
      result.burgers_vector.assign( result.lines.size(), Vec3d{0.,0.,0.} );
      std::vector<int> nearby_count( result.lines.size(), 0 );

      std::vector<std::vector<uint32_t>> vertex_to_lines( n_vertices );
      for(size_t li=0; li<result.lines.size(); li++)
      {
        for(uint32_t v : result.lines[li]) { vertex_to_lines[v].push_back( static_cast<uint32_t>(li) ); }
      }

      std::vector<uint32_t> skeleton_points;
      for(uint32_t v=0; v<n_vertices; v++) { if( alive[v] && !skel_adj[v].empty() ) { skeleton_points.push_back(v); } }

      long n_signal_edges_unassigned = 0;
      for(size_t i=0; i<bc.edges.size(); i++)
      {
        const uint32_t a = bc.edges[i][0], b = bc.edges[i][1];
        const Vec3d midpoint = (mesh.vertices[a] + mesh.vertices[b]) * 0.5;

        uint32_t nearest_v = 0; double nearest_d2 = -1.0;
        for(uint32_t v : skeleton_points)
        {
          const Vec3d d = mesh.vertices[v] - midpoint;
          const double d2 = d.x*d.x + d.y*d.y + d.z*d.z;
          if( nearest_d2 < 0.0 || d2 < nearest_d2 ) { nearest_d2 = d2; nearest_v = v; }
        }
        if( nearest_d2 < 0.0 || vertex_to_lines[nearest_v].empty() ) { ++n_signal_edges_unassigned; continue; }

        for(uint32_t li : vertex_to_lines[nearest_v])
        {
          result.burgers_vector[li] = result.burgers_vector[li] + bc.burgers_vector[i];
          ++nearby_count[li];
        }
      }
      for(size_t li=0; li<result.lines.size(); li++)
      {
        if( nearby_count[li] > 0 ) { result.burgers_vector[li] = result.burgers_vector[li] / static_cast<double>(nearby_count[li]); }
      }

      *n_lines_extracted = static_cast<long>( result.lines.size() );
      *n_junctions = static_cast<long>( result.junction_vertices.size() );
      *n_thinning_removed = n_removed;

      lout << "compute_dxa_dislocation_lines: " << n_hole_atoms << " hole atoms thinned to "
           << (n_hole_atoms - n_removed) << " skeleton atoms (" << n_removed << " removed), "
           << result.lines.size() << " lines extracted ("
           << std::count(result.is_loop.begin(), result.is_loop.end(), uint8_t(1)) << " closed loops), "
           << result.junction_vertices.size() << " junctions, "
           << n_signal_edges_unassigned << "/" << bc.edges.size() << " signal edges unassigned to any line" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline steps (viii)-(ix): extracts 1D dislocation line geometry from the disordered atom
population and detects junctions, via topological graph thinning + tracing (not a surface-mesh
sweep -- see this file's own header comment for the rationale and its tradeoffs). Produces
DXADislocationLines: an ordered vertex-index polyline per line (open segments and closed loops),
one averaged Burgers vector per line (from nearby DXABurgersCircuits signal edges), and the set of
junction atoms (thinned-graph degree >= 3).

Usage example:

compute_dxa_edge_vectors: { struct_field: cna_type, target_structure: BCC }
compute_dxa_tet_classification: {}
compute_dxa_burgers_circuits: { min_burgers_norm: 0.5 }
compute_dxa_dislocation_lines: {}

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_dxa_dislocation_lines)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_dislocation_lines", make_simple_operator< ComputeDXADislocationLines > );
  }

}
