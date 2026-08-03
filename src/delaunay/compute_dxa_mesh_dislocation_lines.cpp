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
#include <exaStamp/delaunay/interface_mesh.h>
#include <exaStamp/delaunay/dxa_edge_vectors.h>
#include <exaStamp/delaunay/dxa_dislocation_lines.h>

#include <algorithm>
#include <array>
#include <functional>
#include <set>
#include <string>
#include <vector>

// DXA pipeline steps (viii)-(ix), take 2: dislocation line tracing seeded from confirmed Burgers
// circuits (compute_dxa_mesh_burgers_circuits), not the raw hole-atom population
// (compute_dxa_dislocation_lines, superseded -- see its own file header and
// src/delaunay/README.md).
//
// OVITO's documented algorithm advances ("sweeps") a confirmed seed circuit step by step along the
// interface mesh, recomputing it at each step and taking its center of mass as the next line
// vertex. A literal implementation of that advancing-front walk (locally replacing one loop vertex
// with an adjacent one via shared triangles, re-checking Burgers non-zero-ness at each step) is a
// substantial piece of computational geometry in its own right. This operator instead reuses this
// codebase's own already-validated topological-thinning machinery (compute_dxa_dislocation_lines'
// core algorithm -- iteratively remove any vertex with graph-degree > 2 whose removal wouldn't
// disconnect its neighbors, i.e. not a Tarjan articulation point, until only curve points and real
// junctions remain), but applied to a MUCH smaller and already-validated graph: only the vertices
// touched by a CONFIRMED signal edge (real Burgers-vector-carrying circuits, not raw noisy atom
// connectivity), pulling in the interface mesh's own full edge connectivity among just those
// vertices to keep the graph properly connected for thinning. ponytail: this approximates "sweep,
// center-of-mass per step" as "thin the validated core region, use its own vertex positions
// directly" -- geometrically similar in spirit (both end up tracing the middle of the tube), but
// not a literal port of the advancing-front algorithm. Revisit if results still don't look right.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXAMeshDislocationLines : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation     , INPUT , REQUIRED );
    ADD_SLOT( InterfaceMesh        , interface_mesh            , INPUT , REQUIRED );
    ADD_SLOT( DXABurgersCircuits   , dxa_mesh_burgers_circuits , INPUT , REQUIRED );
    ADD_SLOT( DXADislocationLines  , dxa_dislocation_lines     , OUTPUT , DocString{"Same slot name as the superseded compute_dxa_dislocation_lines' own output -- deliberate, so write_dxa_dislocation_lines auto-wires without a rebind; don't run both operators in the same pipeline"} );
    ADD_SLOT( long                 , n_mesh_lines_extracted    , OUTPUT , DocString{"Number of dislocation lines extracted (open segments + closed loops)"} );
    ADD_SLOT( long                 , n_mesh_junctions          , OUTPUT , DocString{"Number of skeleton atoms where 3+ lines meet"} );
    ADD_SLOT( long                 , n_mesh_thinning_removed   , OUTPUT , DocString{"Number of core-region atoms removed by topological thinning"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const InterfaceMesh& iface = *interface_mesh;
      const DXABurgersCircuits& bc = *dxa_mesh_burgers_circuits;
      const size_t n_vertices = mesh.vertices.size();

      // core vertex set: only atoms touched by a confirmed signal edge
      std::vector<uint8_t> is_core( n_vertices, 0 );
      for(const auto& e : bc.edges) { is_core[e[0]] = 1; is_core[e[1]] = 1; }

      // core adjacency: the interface mesh's own full edge connectivity, restricted to core
      // vertices -- keeps the graph properly connected (a confirmed signal edge is just one
      // arbitrary "closing" edge of its own fundamental cycle, not necessarily directly touching
      // its geometric neighbors along the tube; the surrounding mesh edges provide that).
      std::vector<std::vector<uint32_t>> core_adj( n_vertices );
      for(const auto& kv : iface.edge_triangles)
      {
        const uint32_t a = static_cast<uint32_t>( kv.first >> 32 );
        const uint32_t b = static_cast<uint32_t>( kv.first & 0xffffffffu );
        if( is_core[a] && is_core[b] ) { core_adj[a].push_back(b); core_adj[b].push_back(a); }
      }

      long n_core_atoms = 0;
      for(uint32_t v=0; v<n_vertices; v++) { if( is_core[v] ) { ++n_core_atoms; } }

      // topological thinning via iterative articulation-point-safe removal (same technique as
      // compute_dxa_dislocation_lines.cpp -- see that file for the detailed rationale).
      std::vector<uint8_t> alive( n_vertices, 0 );
      for(uint32_t v=0; v<n_vertices; v++) { alive[v] = is_core[v]; }

      auto current_degree = [&]( uint32_t v ) -> int
      {
        int d = 0;
        for(uint32_t w : core_adj[v]) { if( alive[w] ) { ++d; } }
        return d;
      };

      long n_removed = 0;
      bool progress = true;
      while( progress )
      {
        progress = false;

        std::vector<int32_t> disc( n_vertices, -1 ), low( n_vertices, -1 ), parent( n_vertices, -1 );
        std::vector<uint8_t> is_articulation( n_vertices, 0 );
        int32_t timer = 0;

        std::function<void(uint32_t)> dfs = [&]( uint32_t root )
        {
          std::vector<std::pair<uint32_t,size_t>> stack;
          disc[root] = low[root] = timer++;
          stack.push_back( { root, 0 } );
          int32_t root_children = 0;

          while( !stack.empty() )
          {
            const size_t top = stack.size() - 1;
            const uint32_t u = stack[top].first;
            const size_t idx = stack[top].second;

            if( idx < core_adj[u].size() )
            {
              const uint32_t v = core_adj[u][idx];
              stack[top].second = idx + 1;
              if( !alive[v] ) { continue; }
              if( disc[v] < 0 )
              {
                parent[v] = static_cast<int32_t>(u);
                disc[v] = low[v] = timer++;
                if( u == root ) { ++root_children; }
                stack.push_back( { v, 0 } );
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
            break;
          }
        }
      }

      // trace the thinned graph
      std::vector<std::vector<uint32_t>> skel_adj( n_vertices );
      for(uint32_t v=0; v<n_vertices; v++)
      {
        if( !alive[v] ) { continue; }
        for(uint32_t w : core_adj[v]) { if( alive[w] ) { skel_adj[v].push_back(w); } }
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
            if( visited_edges.count(k) ) { break; }
            visited_edges.insert(k);
            path.push_back(nxt);
            prev = cur; cur = nxt;
          }
          result.lines.push_back( std::move(path) );
          result.is_loop.push_back(0);
        }
      }

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

      // per-line Burgers vector: average every confirmed signal edge whose midpoint's nearest
      // skeleton point belongs to that line -- should be tight/direct here since skeleton points
      // are themselves confirmed-signal-adjacent by construction.
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

      *n_mesh_lines_extracted = static_cast<long>( result.lines.size() );
      *n_mesh_junctions = static_cast<long>( result.junction_vertices.size() );
      *n_mesh_thinning_removed = n_removed;

      lout << "compute_dxa_mesh_dislocation_lines: " << n_core_atoms << " core atoms (from "
           << bc.edges.size() << " signal edges) thinned to " << (n_core_atoms - n_removed)
           << " skeleton atoms (" << n_removed << " removed), " << result.lines.size() << " lines extracted ("
           << std::count(result.is_loop.begin(), result.is_loop.end(), uint8_t(1)) << " closed loops), "
           << result.junction_vertices.size() << " junctions, "
           << n_signal_edges_unassigned << "/" << bc.edges.size() << " signal edges unassigned to any line" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline steps (viii)-(ix), corrected: extracts 1D dislocation line geometry seeded from
confirmed Burgers circuits (compute_dxa_mesh_burgers_circuits), not the raw hole-atom population --
see this file's own header comment for the approach and its relation to OVITO's documented
advancing-front sweep. Produces DXADislocationLines: an ordered vertex-index polyline per line,
one averaged Burgers vector per line, and the set of junction atoms.

Usage example:

compute_atomistic_interface_mesh: { target_structure: BCC }
compute_dxa_mesh_burgers_circuits: { min_burgers_norm: 0.5 }
compute_dxa_mesh_dislocation_lines: {}

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_dxa_mesh_dislocation_lines)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_mesh_dislocation_lines", make_simple_operator< ComputeDXAMeshDislocationLines > );
  }

}
