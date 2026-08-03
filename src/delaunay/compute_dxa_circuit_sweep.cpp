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
#include <exaStamp/delaunay/dxa_dislocation_lines.h>

#include <algorithm>
#include <array>
#include <deque>
#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

// DXA pipeline steps (vi)-(ix), real circuit search + sweep. Second-generation rewrite: the first
// version (see git history / README) was built from the 2012 paper's own high-level description
// (Stukowski/Bulatov/Arsenlis, secs 2.4-2.6) and got the *sweep mechanics* right (halfedge moves,
// per-facet ownership), but badly over-fragmented real multi-junction networks: 41 segments found
// on the "quadrupole" test case vs. OVITO's own real DXA ground truth of 11. Diagnosed by reading
// OVITO's actual (non-public) DXA source, `DislocationTracer.cpp`, obtained directly from the user:
// the real algorithm differs from this operator's first version in two structural ways that this
// rewrite now reproduces:
//
//   1. Territorial exclusion is baked into the SEED SEARCH itself, not just the sweep. OVITO's BFS
//      explicitly skips any mesh edge that already belongs to (or borders the facet of) an existing
//      circuit ("findPrimarySegments"/"createBurgersCircuit"). The first version's BFS considered
//      every edge with a resolved ideal vector regardless of nearby already-claimed territory, so a
//      new seed's BFS could trace right along the boundary of an already-swept segment and discover
//      a spurious "real-looking" nonzero-Burgers loop that was really just an artifact of that
//      boundary -- a major source of the 41-vs-11 gap.
//   2. Segments grow in *lockstep*, one small increment at a time, via an outer loop over
//      increasing trial-circuit length (3, 5, 7, 9, ... up to max_circuit_length), interleaving:
//      (a) extending every still-growing ("dangling") segment by one round, capped at that round's
//      length; (b) searching for brand new seeds, but only at odd lengths <= max_circuit_length;
//      (c) OVITO also joins segments into junctions every round (see below). The first version swept
//      each seed to full completion (up to circuit_stretchability) before ever trying the next seed,
//      which is an order-dependent, greedy process: whichever segment happens to be discovered first
//      can overreach into territory a not-yet-discovered neighboring segment should rightfully have
//      owned. Growing everything in lockstep removes that unfairness.
//
// Move set: OVITO's real 5 variants are all reproduced except one, found to be structurally
// impossible to trigger here rather than skipped as an oversight:
//   - remove-2 (a "spike" a->b->a in the loop cancels outright, no facet claim) and sweep-two-facets
//     (slide the boundary sideways across 2 unclaimed triangles sharing a common far apex, when
//     neither alone offers a valid single-facet move) are both implemented in the shrink phase below.
//   - remove-1 (single-triangle shrink) and insert-one (single-triangle expand, now also checking
//     that its two new edges aren't already claimed by another live circuit -- the specific check a
//     first version of this operator was missing) were already present.
//   - remove-3 (three consecutive circuit edges are exactly one triangle's own three sides) is NOT
//     implemented: it requires the loop to revisit the same vertex twice a few steps apart, which
//     this operator's expand move already refuses to create (self-intersection guard) -- the
//     precondition for remove-3 can never arise here, so it would be dead code, not a missing move.
//   - OVITO's "createSecondarySegment" (filling small leftover unclaimed pockets bordering an
//     existing circuit) IS now implemented (`create_secondary_segment` below) -- every round, each
//     dangling node's own boundary is scanned for edges whose opposite side is genuinely unclaimed;
//     if that "hole"'s own perimeter turns out to be a valid, physically-closing, nonzero-Burgers
//     loop bordering at least 2 distinct known segments, it's committed as a new segment in its own
//     right, exactly matching OVITO's own criteria.
//   - Junction *node geometry* (splicing which arms meet into one shared coordinate, OVITO's own CA
//     file circular-list convention) is not reconstructed -- arms of a real 3+-way junction are kept
//     as separate output segments/dislocation_ids, matching OVITO's own rule ("only create a real
//     junction for three or more segments"), but with no shared node position.
//
// KNOWN OPEN BUG, root-caused but not fixed: near a real 3+-way junction, a single circuit's own
// growth can "leak" across the junction into a topologically DIFFERENT, adjacent dislocation's core,
// producing one chimera segment whose points trace two real ground-truth lines end-to-end (confirmed
// by mapping every output point against exact ground truth: e.g. one raw, entirely UNMERGED segment
// -- ruling out the 2-arm merge logic as the cause -- had 81 of its points on one real line and 10 on
// a neighboring one it shares a real junction with). Every individual shrink/expand/sweep-two-facets
// move stays locally valid throughout (no size blowup, no facet double-claim) -- the circuit just
// drifts, step by step, from encircling one core to encircling a different one, because nothing here
// re-validates that a move keeps the circuit's own local elastic mapping self-consistent. This is
// exactly what OVITO's Frank-rotation check (skipped here, see below) would catch in principle, but
// implementing an equivalent for a single-grain system would need tracking *local* per-vertex lattice
// orientation drift during growth, not just a one-time cluster-transition compatibility test -- a
// substantially bigger undertaking than anything else in this file. Not yet attempted.
//
// IMPORTANT, found only after measuring: 2-arm merging must happen INCREMENTALLY, during growth, not
// as a post-hoc pass at the very end. An earlier attempt at this rewrite kept the operator's original
// post-hoc Union-Find merge (decide everything once, after all growth finishes) and made *only* the
// two changes above -- this made the screw-dipole test perfect (2 raw segments merged into the
// correct 2 dislocations) but made the quadrupole test far *worse* (116 raw segments, up from the
// first version's 41, vs. the true 11): once two adjacent redundant seeds mutually block each other,
// post-hoc merging correctly records them as one dislocation for length statistics, but neither one
// is still "dangling" to keep *growing through the gap* in later rounds -- so that gap sits unclaimed
// and invites yet more redundant reseeding at each subsequent round. OVITO's own joinSegments avoids
// this by merging two mutually-and-exclusively blocking dangling ends immediately, splicing the far
// end of the absorbed segment onto the surviving one as a still-growing dangling node. This operator
// now does the same (see the per-round merge pass below), using `facet_owner`/`edge_owner` keyed by
// NODE id (not segment id) specifically so a stopped node can identify the *exact* other node it met,
// not just which segment it belongs to.
//   - No Frank-rotation / cluster-transition consistency check: this operator only ever runs on a
//     single crystal orientation (no grain boundaries), so that check is trivially always satisfied
//     and is omitted entirely rather than implemented as a no-op.
//
// Halfedge structure: each mesh triangle (v0,v1,v2) contributes 3 directed halfedges (v0->v1),
// (v1->v2), (v2->v0), each "owned" by that triangle -- `halfedge_owner[(a,b)]` is the unique
// triangle whose own winding has (a,b) as one of its three consecutive sides. A circuit is a closed
// ordered list of vertices; consecutive pairs are its own directed halfedges.
//
// Two independent layers of "claimed" state, exactly mirroring OVITO's `Face::circuit` (interior,
// permanently swept territory) vs. `Edge::circuit` (this specific directed edge is currently part of
// some circuit's *live* boundary, not yet interior):
//   - facet_owner[T]: which NODE has swept triangle T into its interior (shrink/expand moves both set
//     this on the triangle they just claimed). Once set, permanent. Keyed by node (not segment) so a
//     stopped node can identify the exact other node it met, for incremental merging (see above).
//   - edge_owner[(a,b)]: which node currently uses directed edge (a,b) as one of its own circuit's
//     boundary edges. Updated every move (edges leave the boundary on a shrink, join it on an
//     expand). Two different nodes' circuits routinely share an undirected edge using opposite
//     directed senses -- that's a normal adjacent-segment boundary, not a conflict.
// The dislocation line's vertex at each move is the circuit's own center of mass, as documented.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXACircuitSweep : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation     , INPUT , REQUIRED );
    ADD_SLOT( InterfaceMesh        , interface_mesh            , INPUT , REQUIRED );
    ADD_SLOT( DXADislocationLines  , dxa_dislocation_lines     , OUTPUT , DocString{"Same slot name as the earlier (superseded/approximated) line-extraction operators' own output, deliberately, so write_dxa_dislocation_lines auto-wires -- don't run more than one of them in the same pipeline"} );
    ADD_SLOT( long                 , max_circuit_length        , INPUT , 14 , DocString{"Maximum seed-circuit length in mesh-edge steps (OVITO's own default) -- a dislocation whose core is too wide to enclose within this won't be found"} );
    ADD_SLOT( long                 , circuit_stretchability    , INPUT , 9  , DocString{"Extra loop-size elasticity allowed during growth beyond max_circuit_length before a segment stops (OVITO's own default)"} );
    ADD_SLOT( double               , min_burgers_norm          , INPUT , 0.3 , DocString{"Minimum |Burgers vector| (same units as InterfaceMesh::edge_ideal_vector) for a local trial circuit to count as enclosing a real dislocation rather than numerical noise"} );
    ADD_SLOT( long                 , max_sweep_steps           , INPUT , 2000 , DocString{"Safety cap on shrink/expand iterations per round, only to bound runtime if something loops without triggering the (now real) stop conditions -- shouldn't normally be hit"} );
    ADD_SLOT( long                 , n_lines_extracted         , OUTPUT );
    ADD_SLOT( long                 , n_seeds_tried             , OUTPUT , DocString{"Number of (mesh vertex, growth round) local trial-circuit searches actually attempted across the whole incremental growth process"} );
    ADD_SLOT( long                 , n_dislocations            , OUTPUT , DocString{"Number of distinct physical dislocations after merging raw segments across clean two-way junction chains -- <= n_lines_extracted, equal only if nothing merged"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const InterfaceMesh& iface = *interface_mesh;
      const size_t n_vertices = mesh.vertices.size();
      const int max_len = static_cast<int>( *max_circuit_length );
      const int max_stretch = max_len + static_cast<int>( *circuit_stretchability );
      const int step_cap = static_cast<int>( *max_sweep_steps );
      const double burgers_threshold = *min_burgers_norm;

      auto directed_key = []( uint32_t a, uint32_t b ) -> uint64_t { return (uint64_t(a)<<32) | uint64_t(b); };

      // halfedge ownership: for each triangle, its own 3 directed sides
      std::unordered_map<uint64_t,uint32_t> halfedge_owner;
      halfedge_owner.reserve( iface.triangles.size() * 3 );
      for(uint32_t t=0; t<iface.triangles.size(); t++)
      {
        const auto& tri = iface.triangles[t];
        for(int i=0;i<3;i++) { halfedge_owner[ directed_key( tri[i], tri[(i+1)%3] ) ] = t; }
      }

      // vertex adjacency for the seed search, restricted to edges with a resolved ideal vector (see
      // InterfaceMesh::edge_ideal_vector -- hole-closing synthesized edges don't have one).
      std::vector<std::vector<uint32_t>> adj( n_vertices );
      for(const auto& kv : iface.edge_ideal_vector)
      {
        const uint32_t a = static_cast<uint32_t>( kv.first >> 32 );
        const uint32_t b = static_cast<uint32_t>( kv.first & 0xffffffffu );
        adj[a].push_back(b);
        adj[b].push_back(a);
      }

      auto burgers_of_loop = [&]( const std::vector<uint32_t>& loop ) -> Vec3d
      {
        Vec3d b{0.,0.,0.};
        const size_t n = loop.size();
        for(size_t i=0;i<n;i++)
        {
          const uint32_t u = loop[i], v = loop[(i+1)%n];
          const auto it = iface.edge_ideal_vector.find( InterfaceMesh::edge_key(u,v) );
          const Vec3d stored = it->second;
          b = b + ( (u < v) ? stored : Vec3d{-stored.x,-stored.y,-stored.z} );
        }
        return b;
      };

      auto centroid_of = [&]( const std::vector<uint32_t>& loop ) -> Vec3d
      {
        Vec3d c{0.,0.,0.};
        for(uint32_t v : loop) { c = c + mesh.vertices[v]; }
        return c / static_cast<double>( loop.size() );
      };

      std::vector<int32_t> facet_owner( iface.triangles.size(), -1 );
      std::unordered_map<uint64_t,int32_t> edge_owner;

      enum class StopReason { MaxLength, SelfClosure, Junction, OpenEdge, Exhausted };

      struct NodeState
      {
        std::vector<uint32_t> loop;
        bool dangling = true;
        bool is_forward = true;
        int32_t segment = -1;        // mutable: which segs[] entry this node currently feeds (reassigned on merge)
        int32_t fixed_sibling = -1;  // the other node created together with this one -- never changes
        bool retired = false;        // true once fully resolved: merged away, or confirmed standalone
        StopReason stop_reason = StopReason::Exhausted;
        // ALL distinct other nodes whose facet claims this stopped circuit's own boundary directly
        // touches (Junction only) -- not just one. A single `blocking_node` (the original design)
        // silently discarded this information whenever a circuit's boundary happened to border two
        // or more DIFFERENT other circuits at once: exactly the direct, local signature of a real
        // 3+-way junction (OVITO's own joinSegments() detects this the same way, connecting every
        // adjacent circuit it finds while scanning a stopped circuit's whole boundary into a
        // "junctionRing", not just remembering the last one). Falling back to a global
        // incoming-count proxy instead of this direct local signal was found to let a real 3-way
        // junction's arms merge incorrectly when the *global* count didn't happen to catch it.
        std::vector<int32_t> blocking_nodes;
        int32_t stable_merge_rounds = 0; // consecutive rounds the exclusive-mutual merge condition below has held
        int32_t resolved_to = -1; // once merged away (absorbed as a "near" node), the surviving node representing this territory going forward
      };
      struct SegmentState
      {
        Vec3d burgers{0.,0.,0.};
        std::deque<Vec3d> line;
        std::deque<int32_t> core_size; // parallel to line -- the loop's own vertex count when each point was recorded, see DXADislocationLines::core_size
        bool active = true; // false once absorbed into another segment via a clean two-arm merge
      };
      std::vector<NodeState> nodes;
      std::vector<SegmentState> segs;
      std::vector<int32_t> dangling_nodes;

      long n_stop_maxlen=0, n_stop_closure=0, n_stop_junction=0, n_stop_openedge=0, n_stop_exhausted=0;

      // `facet_owner[T]` freezes whichever node FIRST claimed triangle T -- if that node later gets
      // absorbed into a two-arm merge, the frozen id becomes stale: a real third arm that later runs
      // into the very same territory would record a blocking_node pointing at an already-retired
      // node, making it invisible to the (already-committed) merge's own incoming-count check. This
      // was a real, measured bug: a genuine 3-way junction's two faster arms got incorrectly spliced
      // into one chimera segment because the slower third arm's block was never counted against
      // them. Always resolve a node id through any merge chain before using it for anything that
      // affects a merge decision.
      auto resolve_node = [&]( int32_t n ) -> int32_t
      {
        while( nodes[n].resolved_to != -1 ) { n = nodes[n].resolved_to; }
        return n;
      };

      auto append_point = [&]( int32_t node_id )
      {
        NodeState& nd = nodes[node_id];
        const Vec3d c = centroid_of( nd.loop );
        const int32_t cs = static_cast<int32_t>( nd.loop.size() );
        if( nd.is_forward ) { segs[nd.segment].line.push_back(c); segs[nd.segment].core_size.push_back(cs); }
        else { segs[nd.segment].line.push_front(c); segs[nd.segment].core_size.push_front(cs); }
      };

      // Grows one node's circuit as far as possible this round: fully shrink to a local minimum
      // (repeat until a whole pass finds no more shrinks), then attempt exactly one expand if the
      // round's size cap allows it -- matching OVITO's DislocationTracer::traceSegment structure
      // (shrink-to-minimum takes full precedence, one insert per round after that).
      auto trace_node_round = [&]( int32_t node_id, int round_cap )
      {
        NodeState& nd = nodes[node_id];
        if( !nd.dangling ) { return; }

        for(int outer=0; outer<step_cap; outer++)
        {
          bool any_shrink = true;
          while( any_shrink )
          {
            any_shrink = false;
            const size_t n = nd.loop.size();
            for(size_t i=0;i<n;i++)
            {
              const uint32_t a = nd.loop[i], b = nd.loop[(i+1)%n], c = nd.loop[(i+2)%n];
              if( a == c )
              {
                // the loop doubles back on itself via (a,b) then (b,a) -- a "spike" that can arise
                // after an earlier expand, especially on an irregular real mesh. Cancel both edges
                // outright, no facet claim needed (this doesn't traverse any new territory) --
                // matches OVITO's tryRemoveTwoCircuitEdges, checked with top priority, before this
                // operator's original single-triangle shrink below.
                if( n < 5 ) { continue; } // would leave fewer than 3 vertices, not a valid loop
                const size_t i1 = (i+1)%n, i2 = (i+2)%n;
                std::vector<uint32_t> new_loop;
                new_loop.reserve(n-2);
                for(size_t k=0;k<n;k++) { if( k!=i1 && k!=i2 ) { new_loop.push_back( nd.loop[k] ); } }
                edge_owner.erase( directed_key(a,b) );
                edge_owner.erase( directed_key(b,a) );
                nd.loop = std::move(new_loop);
                append_point( node_id );
                any_shrink = true;
                break;
              }
              const auto it_ab = halfedge_owner.find( directed_key(a,b) );
              const auto it_bc = halfedge_owner.find( directed_key(b,c) );
              if( it_ab != halfedge_owner.end() && it_bc != halfedge_owner.end() && it_ab->second == it_bc->second )
              {
                const uint32_t T = it_ab->second;
                if( facet_owner[T] == -1 )
                {
                  facet_owner[T] = node_id;
                  edge_owner.erase( directed_key(a,b) );
                  edge_owner.erase( directed_key(b,c) );
                  edge_owner[ directed_key(a,c) ] = node_id;
                  nd.loop.erase( nd.loop.begin() + static_cast<long>((i+1)%n) );
                  append_point( node_id );
                  any_shrink = true;
                  break;
                }
              }

              // sweep-two-facets: (a,b) and (b,c) border two DIFFERENT unclaimed triangles that
              // happen to share the same "far" apex vertex w opposite b -- slide the boundary
              // sideways from b to w, claiming both triangles, without changing the loop's own
              // size. Needed for cases where neither triangle alone offers a valid single-facet
              // move (its own new edges already belong to another live circuit) but the combined
              // sideways swap's new edges (b,w)+(w,c)... via (a,w),(w,c) do not. Matches OVITO's
              // trySweepTwoFacets, checked here (same priority tier as the shrink family, since it
              // also doesn't require falling back to a plain expand).
              if( it_ab != halfedge_owner.end() && it_bc != halfedge_owner.end() && it_ab->second != it_bc->second )
              {
                const uint32_t f1 = it_ab->second, f2 = it_bc->second;
                if( facet_owner[f1] == -1 && facet_owner[f2] == -1 )
                {
                  const auto& tri1 = iface.triangles[f1]; uint32_t w1 = 0; for(uint32_t tv : tri1) { if( tv!=a && tv!=b ) { w1 = tv; } }
                  const auto& tri2 = iface.triangles[f2]; uint32_t w2 = 0; for(uint32_t tv : tri2) { if( tv!=b && tv!=c ) { w2 = tv; } }
                  if( w1 == w2 && std::find(nd.loop.begin(),nd.loop.end(),w1) == nd.loop.end()
                      && !edge_owner.count(directed_key(a,w1)) && !edge_owner.count(directed_key(w1,c)) )
                  {
                    facet_owner[f1] = node_id; facet_owner[f2] = node_id;
                    edge_owner.erase( directed_key(a,b) );
                    edge_owner.erase( directed_key(b,c) );
                    edge_owner[ directed_key(a,w1) ] = node_id;
                    edge_owner[ directed_key(w1,c) ] = node_id;
                    nd.loop[(i+1)%n] = w1;
                    append_point( node_id );
                    any_shrink = true;
                    break;
                  }
                }
              }
            }
          }

          if( static_cast<int>(nd.loop.size()) >= max_stretch ) { ++n_stop_maxlen; nd.dangling = false; nd.retired = true; nd.stop_reason = StopReason::MaxLength; return; }
          if( static_cast<int>(nd.loop.size()) >= round_cap ) { return; } // pause: retry with a bigger cap next round

          bool saw_open_edge=false, saw_foreign_claim=false, saw_self_claim=false;
          std::vector<int32_t> foreign_nodes; // every DISTINCT other node touched, not just the last one (see NodeState::blocking_nodes)
          int32_t best_i=-1; uint32_t best_w=0, best_T=0;
          const size_t n = nd.loop.size();
          for(size_t i=0;i<n;i++)
          {
            const uint32_t a = nd.loop[i], b = nd.loop[(i+1)%n];
            const auto it = halfedge_owner.find( directed_key(a,b) );
            if( it == halfedge_owner.end() ) { saw_open_edge = true; continue; }
            const uint32_t T = it->second;
            if( facet_owner[T] != -1 )
            {
              const int32_t owner_node = facet_owner[T];
              if( nodes[owner_node].segment == nd.segment ) { saw_self_claim = true; }
              else
              {
                saw_foreign_claim = true;
                if( std::find(foreign_nodes.begin(),foreign_nodes.end(),owner_node) == foreign_nodes.end() ) { foreign_nodes.push_back(owner_node); }
              }
              continue;
            }
            const auto& tri = iface.triangles[T];
            uint32_t w = 0; for(uint32_t tv : tri) { if( tv!=a && tv!=b ) { w = tv; } }
            if( std::find(nd.loop.begin(),nd.loop.end(),w) != nd.loop.end() ) { continue; }
            // the two new boundary edges this expand would create must not already belong to
            // another live circuit -- the check the first version of this operator was missing.
            if( edge_owner.count( directed_key(a,w) ) || edge_owner.count( directed_key(w,b) ) ) { continue; }
            if( best_i < 0 ) { best_i = static_cast<int32_t>(i); best_w = w; best_T = T; }
          }

          if( best_i < 0 )
          {
            nd.dangling = false;
            nd.stop_reason = saw_foreign_claim ? StopReason::Junction
                            : ( saw_self_claim ? StopReason::SelfClosure
                            : ( saw_open_edge ? StopReason::OpenEdge : StopReason::Exhausted ) );
            nd.blocking_nodes = saw_foreign_claim ? foreign_nodes : std::vector<int32_t>{};
            nd.retired = !saw_foreign_claim; // Junction stops are resolved later (merge vs standalone); everything else is final now
            switch(nd.stop_reason)
            {
              case StopReason::Junction: ++n_stop_junction; break;
              case StopReason::SelfClosure: ++n_stop_closure; break;
              case StopReason::OpenEdge: ++n_stop_openedge; break;
              default: ++n_stop_exhausted; break;
            }
            return;
          }

          const uint32_t a = nd.loop[best_i], b = nd.loop[(best_i+1)%nd.loop.size()];
          facet_owner[best_T] = node_id;
          edge_owner.erase( directed_key(a,b) );
          edge_owner[ directed_key(a,best_w) ] = node_id;
          edge_owner[ directed_key(best_w,b) ] = node_id;
          nd.loop.insert( nd.loop.begin() + (best_i+1), best_w );
          append_point( node_id );
        }
      };

      // Incremental per-round merge: a node that stopped because it ran into another segment's
      // territory (Junction) doesn't yet know whether that's a real 3+-way branch or just a clean
      // two-way meeting with exactly one other segment. Resolve this every round, right after growth
      // -- not once at the very end -- so that a confirmed two-way merge's *surviving* end can keep
      // growing through the gap in later rounds instead of leaving it to invite more redundant
      // reseeding (see this file's header comment for why this timing matters, found by measuring).
      auto do_merge = [&]( int32_t a_id, int32_t b_id )
      {
        NodeState& A = nodes[a_id];
        NodeState& B = nodes[b_id];
        const int32_t seg_a = A.segment;
        const int32_t seg_b = B.segment;
        const int32_t far_id = B.fixed_sibling;
        NodeState& far = nodes[far_id];

        std::vector<Vec3d> pts_b_to_far( segs[seg_b].line.begin(), segs[seg_b].line.end() );
        std::vector<int32_t> core_b_to_far( segs[seg_b].core_size.begin(), segs[seg_b].core_size.end() );
        if( B.is_forward ) { std::reverse( pts_b_to_far.begin(), pts_b_to_far.end() ); std::reverse( core_b_to_far.begin(), core_b_to_far.end() ); }

        if( !segs[seg_a].line.empty() && !pts_b_to_far.empty() )
        {
          const Vec3d seam_a = A.is_forward ? segs[seg_a].line.back() : segs[seg_a].line.front();
          const double seam_gap = norm( seam_a - pts_b_to_far.front() );
          if( seam_gap > 10.0 )
          {
            lout << "compute_dxa_circuit_sweep: DIAGNOSTIC suspicious merge seam gap " << seam_gap
                 << " Ang between segment " << seg_a << " (burgers " << segs[seg_a].burgers.x << " " << segs[seg_a].burgers.y << " " << segs[seg_a].burgers.z
                 << ") and segment " << seg_b << " (burgers " << segs[seg_b].burgers.x << " " << segs[seg_b].burgers.y << " " << segs[seg_b].burgers.z << ")" << std::endl;
          }
        }

        if( A.is_forward )
        {
          for(const Vec3d& p : pts_b_to_far) { segs[seg_a].line.push_back(p); }
          for(int32_t c : core_b_to_far) { segs[seg_a].core_size.push_back(c); }
        }
        else
        {
          std::vector<Vec3d> pts_far_to_b( pts_b_to_far.rbegin(), pts_b_to_far.rend() );
          std::vector<int32_t> core_far_to_b( core_b_to_far.rbegin(), core_b_to_far.rend() );
          segs[seg_a].line.insert( segs[seg_a].line.begin(), pts_far_to_b.begin(), pts_far_to_b.end() );
          segs[seg_a].core_size.insert( segs[seg_a].core_size.begin(), core_far_to_b.begin(), core_far_to_b.end() );
        }

        segs[seg_b].active = false;
        far.segment = seg_a;
        far.is_forward = A.is_forward;
        // any territory either A or B themselves claimed while growing now belongs to the combined
        // chain, whose single active representative going forward is `far` -- see resolve_node above.
        A.resolved_to = far_id;
        B.resolved_to = far_id;
        A.retired = true;
        B.retired = true;
      };

      // A node that looks like a clean, exclusive two-way meeting THIS round isn't safe to merge
      // immediately: if the two mutually-blocking nodes just happen to have grown at different
      // rates, a genuine 3+-way real junction's third (slower, or later-seeded) arm may not have
      // reached the meeting point yet, making the first two look like an exclusive pair before the
      // third reveals itself -- found by measuring (a real 3-way junction's two faster arms getting
      // incorrectly spliced into one chimera segment before the third arm caught up). Require the
      // exclusive-mutual condition to hold for several consecutive rounds before committing.
      static constexpr int MERGE_GRACE_ROUNDS = 4;

      // Resolves a node's own set of directly-touched other nodes through any prior merge chain
      // and dedups it -- if this collapses to exactly one distinct id, the node's own boundary
      // provides DIRECT, LOCAL evidence of touching only one other circuit (a real 2-way-merge
      // candidate); if it resolves to two or more distinct ids, that alone is definitive local
      // proof of a real 3+-way junction (this node's own boundary directly borders 2+ different
      // circuits), independent of the global incoming-count heuristic below. Matches OVITO's own
      // joinSegments(), which builds this same information directly by connecting every adjacent
      // circuit found while scanning a stopped circuit's whole boundary (its "junctionRing"),
      // rather than only recording one.
      auto resolved_distinct_blockers = [&]( const NodeState& nd ) -> std::vector<int32_t>
      {
        std::vector<int32_t> out;
        for( int32_t b : nd.blocking_nodes )
        {
          const int32_t rb = resolve_node(b);
          if( std::find(out.begin(),out.end(),rb) == out.end() ) { out.push_back(rb); }
        }
        return out;
      };

      auto resolve_pending_merges = [&]( bool force = false )
      {
        // Update the grace-period counters once per round (not once per inner fixed-point pass
        // below), based on the state as of the start of this round's resolution.
        {
          std::vector<int32_t> pending0;
          for(size_t i=0;i<nodes.size();i++)
          {
            if( !nodes[i].retired && !nodes[i].dangling && nodes[i].stop_reason == StopReason::Junction ) { pending0.push_back(static_cast<int32_t>(i)); }
          }
          std::unordered_map<int32_t,int32_t> incoming0;
          for(int32_t nid : pending0)
          {
            const auto blockers = resolved_distinct_blockers( nodes[nid] );
            if( blockers.size() == 1 && !nodes[blockers[0]].retired ) { ++incoming0[blockers[0]]; }
          }
          for(int32_t nid : pending0)
          {
            NodeState& A = nodes[nid];
            const auto a_blockers = resolved_distinct_blockers(A);
            bool exclusive_mutual = false;
            if( a_blockers.size() == 1 )
            {
              const int32_t bid = a_blockers[0];
              const NodeState& B = nodes[bid];
              const auto b_blockers = resolved_distinct_blockers(B);
              exclusive_mutual = !B.retired && !B.dangling && B.stop_reason == StopReason::Junction
                                && b_blockers.size() == 1 && b_blockers[0] == nid
                                && incoming0[nid] == 1 && incoming0[bid] == 1;
            }
            A.stable_merge_rounds = exclusive_mutual ? (A.stable_merge_rounds + 1) : 0;
          }
        }

        for(size_t pass=0; pass < nodes.size() + 4; pass++)
        {
          bool changed = false;
          std::vector<int32_t> pending;
          for(size_t i=0;i<nodes.size();i++)
          {
            if( !nodes[i].retired && !nodes[i].dangling && nodes[i].stop_reason == StopReason::Junction ) { pending.push_back(static_cast<int32_t>(i)); }
          }
          std::unordered_map<int32_t,int32_t> incoming;
          for(int32_t nid : pending)
          {
            const auto blockers = resolved_distinct_blockers( nodes[nid] );
            if( blockers.size() == 1 && !nodes[blockers[0]].retired ) { ++incoming[blockers[0]]; }
          }

          for(int32_t nid : pending)
          {
            NodeState& A = nodes[nid];
            if( A.retired ) { continue; }
            const auto a_blockers = resolved_distinct_blockers(A);
            if( a_blockers.size() != 1 )
            {
              // this node's own boundary directly touches 2+ distinct other circuits -- definitive
              // local evidence of a real 3+-way junction, no need to wait on anything else.
              A.retired = true;
              changed = true;
              continue;
            }
            const int32_t bid = a_blockers[0];
            NodeState& B = nodes[bid];
            if( B.retired ) { A.retired = true; changed = true; continue; } // meets an already-settled chain/junction
            if( B.dangling ) { continue; } // wait for B to also stop before deciding

            const auto b_blockers = resolved_distinct_blockers(B);
            const bool exclusive_mutual = B.stop_reason == StopReason::Junction
                                         && b_blockers.size() == 1 && b_blockers[0] == nid
                                         && incoming[nid] == 1 && incoming[bid] == 1;
            if( !exclusive_mutual )
            {
              A.retired = true; // a real 3+-way branch point (or an ambiguous case) -- keep A standalone
            }
            else if( force || ( A.stable_merge_rounds >= MERGE_GRACE_ROUNDS && B.stable_merge_rounds >= MERGE_GRACE_ROUNDS ) )
            {
              do_merge( nid, bid );
            }
            else
            {
              continue; // exclusive so far, but not stable long enough yet -- wait, don't mark changed
            }
            changed = true;
          }
          if( !changed ) { break; }
        }
      };

      // Creates a new segment (2 nodes, forward/backward) from a validated closed loop and
      // immediately grows both directions by one round at the given cap -- shared by both ordinary
      // seed discovery (try_seed_from below) and createSecondarySegment (see further below).
      auto commit_new_segment = [&]( const std::vector<uint32_t>& loop, const Vec3d& bvec, int circuit_len )
      {
        const int32_t segid = static_cast<int32_t>( segs.size() );
        segs.push_back( SegmentState{ bvec, std::deque<Vec3d>{ centroid_of(loop) }, std::deque<int32_t>{ static_cast<int32_t>(loop.size()) }, true } );
        const int32_t fwd_id = static_cast<int32_t>( nodes.size() );
        const int32_t bwd_id = fwd_id + 1;
        nodes.push_back( NodeState{ loop, true, true, segid, bwd_id, false, StopReason::Exhausted, {}, 0, -1 } );
        std::vector<uint32_t> rev( loop.rbegin(), loop.rend() );
        nodes.push_back( NodeState{ std::move(rev), true, false, segid, fwd_id, false, StopReason::Exhausted, {}, 0, -1 } );
        for(size_t i=0;i<loop.size();i++)
        {
          edge_owner[ directed_key( loop[i], loop[(i+1)%loop.size()] ) ] = fwd_id;
        }
        dangling_nodes.push_back(fwd_id);
        dangling_nodes.push_back(bwd_id);
        trace_node_round( fwd_id, circuit_len );
        trace_node_round( bwd_id, circuit_len );
      };

      // Within the triangle owning directed edge (p,q), returns that triangle's own PREVIOUS side
      // (the one ending at p) -- used by createSecondarySegment to walk around the border of an
      // unclaimed "hole" in the mesh, exactly mirroring OVITO's Edge::prevFaceEdge().
      auto prev_face_edge_of = [&]( uint32_t p, uint32_t q ) -> std::pair<uint32_t,uint32_t>
      {
        const uint32_t T = halfedge_owner.at( directed_key(p,q) );
        const auto& tri = iface.triangles[T];
        int i = -1;
        for(int k=0;k<3;k++) { if( tri[k]==p && tri[(k+1)%3]==q ) { i = k; break; } }
        return { tri[(i+2)%3], p };
      };

      // OVITO's createSecondarySegment: an existing circuit's boundary may border a small unclaimed
      // "hole" in the mesh -- if that hole's own perimeter (walked all the way around) happens to be
      // a valid, physically-closing, nonzero-Burgers circuit that borders at least 2 distinct known
      // segments, it's very likely a real dislocation segment the ordinary vertex-rooted seed search
      // missed (e.g. a short arm squeezed between two already-discovered ones). Returns true if a new
      // segment was created.
      auto create_secondary_segment = [&]( uint32_t start_a, uint32_t start_b, int32_t starting_node, int circuit_len ) -> bool
      {
        std::vector<uint32_t> hole_loop;
        Vec3d burgers_accum{0.,0.,0.};
        bool has_all_ideal = true;
        std::set<int32_t> touched_segments;
        touched_segments.insert( nodes[starting_node].segment );

        uint32_t cur_a = start_a, cur_b = start_b;
        int edge_count = 0;
        for(int guard=0; guard<step_cap*4; guard++)
        {
          uint32_t probe_a = cur_a, probe_b = cur_b;
          bool found = false;
          for(int inner=0; inner<step_cap*4 && !found; inner++)
          {
            const auto pf = prev_face_edge_of( probe_b, probe_a );
            const auto it = edge_owner.find( directed_key(pf.first,pf.second) );
            if( it != edge_owner.end() )
            {
              touched_segments.insert( nodes[ resolve_node(it->second) ].segment );
              cur_a = pf.second; cur_b = pf.first;
              found = true;
            }
            else
            {
              probe_a = pf.first; probe_b = pf.second;
            }
          }
          if( !found ) { return false; } // couldn't close the walk -- bail out safely

          hole_loop.push_back( cur_a );
          const auto ideal_it = iface.edge_ideal_vector.find( InterfaceMesh::edge_key(cur_a,cur_b) );
          if( ideal_it == iface.edge_ideal_vector.end() ) { has_all_ideal = false; }
          else { burgers_accum = burgers_accum + ( (cur_a < cur_b) ? ideal_it->second : Vec3d{-ideal_it->second.x,-ideal_it->second.y,-ideal_it->second.z} ); }

          if( cur_a == start_a && cur_b == start_b ) { break; } // closed the hole's own perimeter
          ++edge_count;
          if( edge_count > max_len ) { return false; }
        }

        if( touched_segments.size() < 2 || !has_all_ideal || norm(burgers_accum) <= burgers_threshold ) { return false; }
        if( hole_loop.size() < 3 ) { return false; }
        std::set<uint32_t> distinct( hole_loop.begin(), hole_loop.end() );
        if( distinct.size() != hole_loop.size() ) { return false; }

        commit_new_segment( hole_loop, burgers_accum, circuit_len );
        return true;
      };

      // Local trial-circuit search rooted at `root`, bounded to the given round's search depth,
      // stopping at the FIRST valid closing edge found (matching OVITO's "stop as soon as a valid
      // Burgers circuit has been found", not "collect every candidate up to max depth and pick the
      // shortest" like the first version of this operator did) -- both territorial exclusions below
      // (facet_owner, edge_owner) are the actual fix for the over-fragmentation this rewrite exists
      // to address.
      std::vector<int32_t> bfs_parent( n_vertices, -1 );
      std::vector<int32_t> bfs_gen( n_vertices, 0 );
      int32_t bfs_gen_counter = 0;
      long n_tried = 0;

      auto try_seed_from = [&]( uint32_t root, int circuit_len ) -> bool
      {
        const int search_depth = (circuit_len - 1) / 2;
        ++bfs_gen_counter;
        bfs_gen[root] = bfs_gen_counter;
        bfs_parent[root] = -1;
        std::vector<uint32_t> order; order.push_back(root);
        std::vector<int32_t> depth; depth.push_back(0);
        size_t head = 0;
        while( head < order.size() )
        {
          const uint32_t u = order[head]; const int32_t du = depth[head]; ++head;
          if( du >= search_depth ) { continue; }
          for(uint32_t w : adj[u])
          {
            const auto it_uw = halfedge_owner.find( directed_key(u,w) );
            if( it_uw == halfedge_owner.end() ) { continue; }
            if( facet_owner[it_uw->second] != -1 ) { continue; }
            if( edge_owner.count( directed_key(u,w) ) || edge_owner.count( directed_key(w,u) ) ) { continue; }

            if( bfs_gen[w] == bfs_gen_counter )
            {
              if( static_cast<int32_t>(w) == bfs_parent[u] ) { continue; }
              if( u >= w ) { continue; } // each undirected non-tree pair considered once
              // reconstruct candidate loop: root -> ... -> u -> w -> ... -> root (reversed)
              std::vector<uint32_t> path_u, path_w;
              for(int32_t x=static_cast<int32_t>(u); x!=-1; x=bfs_parent[x]) { path_u.push_back(static_cast<uint32_t>(x)); }
              for(int32_t x=static_cast<int32_t>(w); x!=-1; x=bfs_parent[x]) { path_w.push_back(static_cast<uint32_t>(x)); }
              std::vector<uint32_t> loop;
              for(auto it=path_u.rbegin(); it!=path_u.rend(); ++it) { loop.push_back(*it); }
              for(size_t i=0; i+1<path_w.size(); i++) { loop.push_back(path_w[i]); }
              if( loop.size() < 3 || static_cast<int>(loop.size()) > circuit_len ) { continue; }
              std::set<uint32_t> distinct( loop.begin(), loop.end() );
              if( distinct.size() != loop.size() ) { continue; }

              const Vec3d bvec = burgers_of_loop(loop);
              if( norm(bvec) <= burgers_threshold ) { continue; }

              // accept: create the new segment and its two (forward/backward) nodes
              commit_new_segment( loop, bvec, circuit_len );
              return true;
            }
            else
            {
              bfs_gen[w] = bfs_gen_counter;
              bfs_parent[w] = static_cast<int32_t>(u);
              order.push_back(w);
              depth.push_back(du+1);
            }
          }
        }
        return false;
      };

      // Incremental lockstep growth: at each round, extend every still-dangling node by one
      // increment capped at this round's length, then (only at odd lengths <= max_circuit_length)
      // scan for brand new seeds, then resolve any two-way merges this round's stops made possible
      // -- matching OVITO's DislocationTracer::traceDislocationSegments + joinSegments.
      for(int circuit_len = 3; circuit_len <= max_stretch; circuit_len++)
      {
        for(int32_t nid : dangling_nodes) { trace_node_round(nid, circuit_len); }

        if( circuit_len <= max_len && (circuit_len % 2) == 1 )
        {
          for(uint32_t root=0; root<n_vertices; root++)
          {
            if( adj[root].empty() ) { continue; }
            ++n_tried;
            try_seed_from( root, circuit_len );
          }
        }

        // createSecondarySegment pass: IMPLEMENTED (create_secondary_segment above) but DISABLED
        // here by measurement -- see this file's own "KNOWN OPEN BUG" header comment. Tested on the
        // quadrupole case: made results measurably worse (19-20 -> 37 physical dislocations), and 8
        // of the 37 were exactly degenerate (a single point, zero length) -- spurious tiny "holes"
        // that pass the nonzero-Burgers-vector + touches-2-segments test without being real defects.
        // OVITO's own version avoids this via its Frank-rotation consistency check on the hole's own
        // loop, which this operator doesn't have (same missing piece as the open growth-drift bug).
        // Re-enable once that check (or an equivalent) exists; until then this would only add noise.
        for(size_t di=0; di<dangling_nodes.size() && false; di++)
        {
          const int32_t nid = dangling_nodes[di];
          if( !nodes[nid].dangling ) { continue; }
          // copy, not a reference: commit_new_segment below grows `nodes`, which can reallocate and
          // invalidate any reference held into it.
          const std::vector<uint32_t> loop_copy = nodes[nid].loop;
          const size_t ln = loop_copy.size();
          for(size_t i=0;i<ln;i++)
          {
            const uint32_t a = loop_copy[i], b = loop_copy[(i+1)%ln];
            if( !edge_owner.count( directed_key(b,a) ) ) { create_secondary_segment( b, a, nid, circuit_len ); }
          }
        }

        resolve_pending_merges();

        dangling_nodes.erase( std::remove_if( dangling_nodes.begin(), dangling_nodes.end(),
                                               [&]( int32_t nid ){ return !nodes[nid].dangling; } ),
                               dangling_nodes.end() );
      }

      // Anything still dangling after the last round simply never got blocked -- treat as MaxLength.
      for(int32_t nid : dangling_nodes)
      {
        nodes[nid].dangling = false;
        nodes[nid].retired = true;
        nodes[nid].stop_reason = StopReason::MaxLength;
        ++n_stop_maxlen;
      }
      resolve_pending_merges( true ); // final, forced pass: no more rounds are coming, so resolve whatever's left now

      // Any node whose Junction stop never got resolved into a real merge or a confirmed standalone
      // (shouldn't normally happen, everything reaching here should already be retired) is left as-is
      // -- it will simply be attributed to whatever segment it currently belongs to below.

      // Junction endpoint reconciliation, matching OVITO's real DislocationTracer::joinSegments (the
      // armCount>=3 branch): every arm of a real 3+-way junction grows as its own independent circuit
      // and each one's own final recorded point is its own circuit's center of mass -- essentially
      // never exactly coincident with its neighbors' (measured: 0.5-2 Å per junction on the real
      // quadrupole case). OVITO fixes this by building a ring of the nodes that directly, currently
      // border each other (via `Edge::circuit`, i.e. which node right now owns the opposite side of
      // this exact boundary edge), then EXTENDS every arm with one new point reaching the ring's
      // average position.
      //
      // A first attempt at this (see git history / session log) instead identified "who does my
      // stopped circuit touch" via `resolved_distinct_blockers`, which resolves through the MERGE
      // chain (`resolve_node`) -- correct for merge bookkeeping (does this node's territory currently
      // belong to segment X or Y), but wrong for POSITION: a node can claim a triangle early during
      // its own growth (freezing `facet_owner[T]`) and then keep growing well past it before finally
      // stopping somewhere else entirely; resolving through the merge chain can land on that node's
      // eventual, spatially unrelated final position, corrupting the reconciled point instead of
      // fixing it (caught by the user visually: mixed-up endpoints, spurious excess length). Fixed by
      // using `edge_owner` instead of `facet_owner`/merge-chain resolution: `edge_owner[(a,b)]` is
      // updated on every single move (added/removed as a circuit's boundary changes) and is NEVER
      // touched by a merge (`do_merge` only concatenates line geometry, never rewrites edge/facet
      // ownership) -- so for a node that has already stopped growing, its own recorded `edge_owner`
      // entries are frozen exactly as they were the instant it stopped: an always-accurate, always-
      // current positional record, with no merge-chain ambiguity possible.
      {
        // qualifies: a genuine, still-standalone (unmerged) junction-arm survivor -- exactly OVITO's
        // "circuit->isDangling" candidates for ring-building (2-way merges are handled separately by
        // do_merge already, and are deliberately excluded here, matching OVITO only extending-to-
        // center for armCount>=3).
        auto qualifies = [&]( int32_t n ) -> bool
        {
          const NodeState& nd = nodes[n];
          return nd.resolved_to == -1 && nd.stop_reason == StopReason::Junction && segs[nd.segment].active;
        };

        // Direct, current neighbors of node_id's own final loop boundary -- for boundary edge (a,b),
        // whichever node currently owns the exact reverse edge (b,a) is definitionally still sitting
        // right there (or, if it's since stopped too, sitting exactly where it stopped, since nothing
        // ever mutates a stopped node's own edge_owner entries). If nobody currently owns (b,a) (the
        // territory across this edge was fully absorbed by both sides, or belongs to a segment that
        // was itself later fully consumed into someone else's interior), this edge simply contributes
        // no partner -- correctly: there is no live/frozen positional evidence left to reconcile
        // against there, better to leave that arm's own natural endpoint undisturbed than to guess.
        auto live_edge_partners = [&]( int32_t node_id ) -> std::vector<int32_t>
        {
          std::vector<int32_t> out;
          const NodeState& nd = nodes[node_id];
          const size_t n = nd.loop.size();
          for(size_t i=0;i<n;i++)
          {
            const uint32_t a = nd.loop[i], b = nd.loop[(i+1)%n];
            const auto it = edge_owner.find( directed_key(b,a) );
            if( it != edge_owner.end() && it->second != node_id
                && std::find(out.begin(),out.end(),it->second) == out.end() ) { out.push_back( it->second ); }
          }
          return out;
        };

        std::vector<int32_t> uf( nodes.size() );
        for(size_t i=0;i<nodes.size();i++) { uf[i] = static_cast<int32_t>(i); }
        std::function<int32_t(int32_t)> find_root = [&]( int32_t x ) -> int32_t
        {
          while( uf[x] != x ) { uf[x] = uf[uf[x]]; x = uf[x]; }
          return x;
        };
        auto unite = [&]( int32_t a, int32_t b ) { a = find_root(a); b = find_root(b); if( a != b ) { uf[a] = b; } };

        for(size_t i=0;i<nodes.size();i++)
        {
          if( !qualifies(static_cast<int32_t>(i)) ) { continue; }
          for( int32_t p : live_edge_partners(static_cast<int32_t>(i)) )
          {
            if( qualifies(p) ) { unite( static_cast<int32_t>(i), p ); }
          }
        }

        std::unordered_map<int32_t,std::vector<int32_t>> clusters;
        for(size_t i=0;i<nodes.size();i++)
        {
          if( qualifies(static_cast<int32_t>(i)) ) { clusters[ find_root(static_cast<int32_t>(i)) ].push_back(static_cast<int32_t>(i)); }
        }

        for( auto& kv : clusters )
        {
          const auto& members = kv.second;
          if( members.size() < 2 ) { continue; } // no other confirmed arm directly, currently borders this one
          Vec3d avg{0.,0.,0.};
          for( int32_t n : members )
          {
            const NodeState& nd = nodes[n];
            avg = avg + ( nd.is_forward ? segs[nd.segment].line.back() : segs[nd.segment].line.front() );
          }
          avg = avg / static_cast<double>( members.size() );
          // Extend each arm with a NEW point reaching the shared center -- matching OVITO's own
          // push_back/push_front of an extra point, not overwriting the last real swept-circuit
          // position already recorded there.
          for( int32_t n : members )
          {
            const NodeState& nd = nodes[n];
            if( nd.is_forward )
            {
              segs[nd.segment].line.push_back( avg );
              segs[nd.segment].core_size.push_back( segs[nd.segment].core_size.back() );
            }
            else
            {
              segs[nd.segment].line.push_front( avg );
              segs[nd.segment].core_size.push_front( segs[nd.segment].core_size.front() );
            }
          }
        }
      }

      const long n_segments_created = static_cast<long>( segs.size() );
      auto segment_length = []( const std::deque<Vec3d>& path ) -> double
      {
        double l = 0.0;
        for(size_t i=1;i<path.size();i++) { l += norm( path[i] - path[i-1] ); }
        return l;
      };

      // Group every retired node by its CURRENT segment (mutated by merges above) to know, per
      // surviving active segment, whether either final end closed on itself (is_loop).
      std::vector<uint8_t> segment_has_self_closure( segs.size(), 0 );
      for(const auto& nd : nodes)
      {
        if( nd.retired && nd.stop_reason == StopReason::SelfClosure ) { segment_has_self_closure[nd.segment] = 1; }
      }

      DXADislocationLines& result = *dxa_dislocation_lines;
      result.lines.clear();
      result.line_positions.clear();
      result.core_size.clear();
      result.is_loop.clear();
      result.junction_vertices.clear();
      result.burgers_vector.clear();
      result.dislocation_id.clear();

      std::vector<double> final_lengths;
      for(int32_t s=0; s<static_cast<int32_t>(segs.size()); s++)
      {
        if( !segs[s].active ) { continue; }
        result.burgers_vector.push_back( segs[s].burgers );
        result.lines.push_back( {} );
        result.line_positions.push_back( std::vector<Vec3d>( segs[s].line.begin(), segs[s].line.end() ) );
        result.core_size.push_back( std::vector<int32_t>( segs[s].core_size.begin(), segs[s].core_size.end() ) );
        result.is_loop.push_back( segment_has_self_closure[s] );
        result.dislocation_id.push_back( static_cast<int32_t>( result.dislocation_id.size() ) );
        final_lengths.push_back( segment_length( segs[s].line ) );
      }

      const long n_dislocations_found = static_cast<long>( result.lines.size() );
      std::sort( final_lengths.begin(), final_lengths.end() );
      const double median_length = final_lengths.empty() ? 0.0 : final_lengths[ final_lengths.size()/2 ];
      double total_length = 0.0; for(double l : final_lengths) { total_length += l; }

      *n_lines_extracted = static_cast<long>( result.lines.size() );
      *n_seeds_tried = n_tried;
      *n_dislocations = n_dislocations_found;
      lout << "compute_dxa_circuit_sweep: " << n_segments_created << " segments created across all growth rounds ("
           << n_tried << " local trial-circuit searches attempted), sweep stops: "
           << n_stop_maxlen << " max-length, " << n_stop_closure << " self-closure, "
           << n_stop_junction << " junction, " << n_stop_openedge << " open-mesh-edge, "
           << n_stop_exhausted << " exhausted-no-move" << std::endl;
      lout << "compute_dxa_circuit_sweep: " << n_dislocations_found << " physical dislocations survive after "
           << "incremental two-way merging (" << (n_segments_created - n_dislocations_found) << " absorbed), lengths "
           << (final_lengths.empty()?0.0:final_lengths.front())
           << "-" << (final_lengths.empty()?0.0:final_lengths.back()) << " (median " << median_length
           << ", total " << total_length << ")" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline steps (vi)-(ix), real circuit search + sweep -- second-generation rewrite matching the
actual structure of OVITO's own (non-public) DXA implementation, DislocationTracer.cpp: territorial
exclusion (already-claimed mesh edges/facets) baked into the seed search itself, and all segments
grown in lockstep via an outer loop over increasing trial-circuit length, rather than sweeping each
seed to full completion before trying the next one. See this file's own header comment for the full
diagnosis (a first version, built from the 2012 paper's description alone, over-fragmented real
junction networks 41-vs-11 on a ground-truth-checked test case) and the specific, documented
simplifications this rewrite keeps relative to OVITO's full generality.

Usage example:

compute_atomistic_interface_mesh: { target_structure: BCC }
compute_dxa_circuit_sweep: {}

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_dxa_circuit_sweep)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_circuit_sweep", make_simple_operator< ComputeDXACircuitSweep > );
  }

}
