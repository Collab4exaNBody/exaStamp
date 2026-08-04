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

#include <exaStamp/grains/grain_segmentation_algo.h>
#include <exaStamp/grains/grain_quat_math.h>
#include <onika/log.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace exaStamp
{
  static constexpr uint64_t GRAIN_BOND_EMPTY = std::numeric_limits<uint64_t>::max();

  // ---- simple union-find (path compression + union-by-size), same role as OVITO's DisjointSet ----
  namespace
  {
    struct UnionFind
    {
      std::vector<int> parent, size;
      explicit UnionFind( size_t n ) : parent(n), size(n,1) { for(size_t i=0;i<n;i++) parent[i]=static_cast<int>(i); }
      int find( int a ) { while( parent[a]!=a ) { parent[a]=parent[parent[a]]; a=parent[a]; } return a; }
      int merge( int a, int b )
      {
        a=find(a); b=find(b); if(a==b) return a;
        if( size[a] < size[b] ) { std::swap(a,b); }
        parent[b]=a; size[a]+=size[b]; return a;
      }
    };

    // ponytail: per-node std::unordered_map adjacency instead of OVITO's own intrusive red-black
    // tree + preallocated edge buffer -- same Node-Pair-Sampling algorithm/result, just simpler
    // machinery; atom coordination numbers here (~12-18) make this a non-issue at this project's
    // scale, add a flatter structure only if this ever needs to scale far beyond a few 1e5 atoms.
    struct GraphNode
    {
      std::unordered_map<int,double> adj;
      double wnode = 0.0;
      bool alive = true;
    };

    struct DendrogramEdge
    {
      int parent, child;
      double distance;       // NPS's own d = wnode[v]/weight(a,v) metric (dimensionless, NOT degrees)
      double disorientation; // degrees, cluster-level, informational only (matches OVITO's own dead-but-harmless field)
      long size = 0;
      double merge_size = 0.0;
    };

    void recompute_wnode( std::vector<GraphNode>& graph, int a )
    {
      double s = 0.0; for( const auto& kv : graph[a].adj ) { s += kv.second; }
      graph[a].wnode = s;
    }

    // Merges child into parent (the higher-degree node survives, same union-by-degree rationale as
    // OVITO's own Graph::contract_edge): transfers/sums child's own edges onto parent, drops child.
    int contract_edge( std::vector<GraphNode>& graph, int a, int b )
    {
      const int parent = ( graph[a].adj.size() >= graph[b].adj.size() ) ? a : b;
      const int child  = ( parent == a ) ? b : a;
      for( const auto& kv : graph[child].adj )
      {
        const int v = kv.first; const double w = kv.second;
        if( v == parent ) { continue; }
        auto it = graph[parent].adj.find(v);
        const double new_w = ( it != graph[parent].adj.end() ) ? ( it->second + w ) : w;
        graph[parent].adj[v] = new_w;
        graph[v].adj[parent] = new_w;
        graph[v].adj.erase(child);
        recompute_wnode( graph, v );
      }
      graph[parent].adj.erase(child);
      graph[child].adj.clear(); // dead now -- must not leave stale entries a re-visited chain index could read
      graph[child].wnode = 0.0;
      graph[child].alive = false;
      recompute_wnode( graph, parent );
      return parent;
    }

    // argmin_v( wnode[v] / weight(a,v) ), ties broken by smallest index -- OVITO's own metric.
    // Returns (inf,-1) for a dead (already-merged-away) node -- defensive: a stale index can still
    // sit on the chain from a push that predates its own death (see contract_edge's own comment).
    std::pair<double,int> nearest_neighbor( const std::vector<GraphNode>& graph, int a )
    {
      if( !graph[a].alive ) { return { std::numeric_limits<double>::infinity(), -1 }; }
      double best_d = std::numeric_limits<double>::infinity(); int best_v = -1;
      for( const auto& kv : graph[a].adj )
      {
        const int v = kv.first; const double w = kv.second;
        const double d = graph[v].wnode / w;
        if( d < best_d || ( d == best_d && v < best_v ) ) { best_d = d; best_v = v; }
      }
      return { best_d, best_v };
    }

    inline double quat_norm( const double* q ) { return std::sqrt( q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] ); }

    // OVITO's own cluster-merge accumulation (GrainSegmentationEngine1::calculate_disorientation):
    // maps child's own running quaternion sum onto parent's current orientation neighborhood, then
    // accumulates (weighted by the child's own pre-mapping norm, i.e. its own accumulated atom
    // count) into parent's running sum. Returns the disorientation between the two (degrees,
    // informational only -- matches OVITO's own harmless-dead-field, see grain_quat_math.h).
    double merge_orientation( int structure_type, double* qsum_parent, const double* qsum_child )
    {
      const double np = quat_norm(qsum_parent);
      double qtarget[4] = { qsum_parent[0]/np, qsum_parent[1]/np, qsum_parent[2]/np, qsum_parent[3]/np };
      const double nc = quat_norm(qsum_child);
      double q[4] = { qsum_child[0]/nc, qsum_child[1]/nc, qsum_child[2]/nc, qsum_child[3]/nc };
      const double disorientation = gb_map_quaternion_onto_target( structure_type, qtarget, q );
      qsum_parent[0] += q[0]*nc; qsum_parent[1] += q[1]*nc; qsum_parent[2] += q[2]*nc; qsum_parent[3] += q[3]*nc;
      return disorientation;
    }

    // Robust (IRLS/least-absolute-deviations) log-log power-law fit + inlier-based threshold pick,
    // following OVITO's ThresholdSelection::Regressor/calculate_threshold exactly (100 fixed IRLS
    // iterations, cutoff=1.5). x=log(merge_size), y=log(distance), weighted by merge_size.
    double auto_select_threshold( const std::vector<DendrogramEdge>& dendrogram )
    {
      const size_t n = dendrogram.size();
      if( n == 0 ) { return -std::numeric_limits<double>::infinity(); }
      std::vector<double> x(n), y(n), w(n);
      for(size_t i=0;i<n;i++) { x[i]=std::log(dendrogram[i].merge_size); y[i]=std::log(dendrogram[i].distance); w[i]=dendrogram[i].merge_size; }

      double gradient = 0.0, intercept = 0.0;
      for(int iter=0; iter<100; iter++)
      {
        double sw=0, swx=0, swy=0, swxx=0, swxy=0;
        for(size_t i=0;i<n;i++) { sw+=w[i]; swx+=w[i]*x[i]; swy+=w[i]*y[i]; swxx+=w[i]*x[i]*x[i]; swxy+=w[i]*x[i]*y[i]; }
        const double denom = sw*swxx - swx*swx;
        if( std::fabs(denom) > 1e-300 )
        {
          gradient = ( sw*swxy - swx*swy ) / denom;
          intercept = ( swy - gradient*swx ) / sw;
        }
        for(size_t i=0;i<n;i++)
        {
          const double residual = y[i] - ( gradient*x[i] + intercept );
          w[i] = dendrogram[i].merge_size / std::max( 1e-4, std::fabs(residual) );
        }
      }

      std::vector<double> abs_residuals(n);
      for(size_t i=0;i<n;i++) { abs_residuals[i] = std::fabs( y[i] - (gradient*x[i]+intercept) ); }
      std::vector<double> sorted_res = abs_residuals; std::sort( sorted_res.begin(), sorted_res.end() );
      const double mad = sorted_res[ sorted_res.size()/2 ];

      double threshold = -std::numeric_limits<double>::infinity();
      const double cutoff = 1.5;
      for(size_t i=0;i<n;i++)
      {
        const double residual = y[i] - ( gradient*x[i] + intercept );
        if( residual < cutoff*mad ) { threshold = std::max( threshold, y[i] ); }
      }
      return threshold;
    }
  }

  void grain_segmentation_nps(
    size_t n_atoms,
    const double * struct_type,
    const double * orientation_mat3,
    const uint64_t * global_id,
    const uint64_t * bond_id,
    const double * bond_distance,
    const double * bond_disorientation,
    const int * bond_count,
    int max_neighbors,
    bool auto_threshold,
    double manual_threshold_log,
    long min_grain_atom_count,
    bool orphan_adoption,
    unsigned int color_seed,
    GrainSegmentationResult & result )
  {
    result = GrainSegmentationResult{};
    result.atom_grain_id.assign( n_atoms, 0 );
    if( n_atoms == 0 ) { return; }

    // global_id -> local index, needed since bond targets are cross-rank-comparable global ids
    std::unordered_map<uint64_t,int> id_to_local;
    id_to_local.reserve( n_atoms*2 );
    for(size_t i=0;i<n_atoms;i++) { id_to_local[ global_id[i] ] = static_cast<int>(i); }

    // Symmetric adjacency over EVERY recorded bond (not just crystalline candidates), built once and
    // reused for orphan adoption below. This matters specifically at a periodic domain boundary: the
    // caller's own bond_id/bond_count/bond_distance arrays are only ever populated from an OWNED
    // atom's own perspective (ghost atoms never get their own bond list computed upstream) -- reading
    // them directly and walking "from an assigned atom to ITS OWN recorded neighbors" therefore can
    // never propagate FROM a ghost, even when that ghost (a periodic image of a real, already-
    // clustered atom) is the geometrically closest assigned neighbor to an orphan on the opposite
    // face. Symmetrizing here once removes that directional gap for every consumer, matching how the
    // crystalline `graph` above is already symmetric by construction.
    std::vector<std::vector<std::pair<int,double>>> all_bonds( n_atoms );
    for(size_t i=0;i<n_atoms;i++)
    {
      const int cnt = bond_count[i];
      for(int k=0;k<cnt;k++)
      {
        const size_t slot = i*static_cast<size_t>(max_neighbors) + static_cast<size_t>(k);
        const uint64_t nid = bond_id[slot];
        if( nid == GRAIN_BOND_EMPTY ) { continue; }
        auto it = id_to_local.find(nid);
        if( it == id_to_local.end() ) { continue; }
        const int j = it->second;
        if( j == static_cast<int>(i) ) { continue; }
        const double dist = bond_distance[slot];
        all_bonds[i].push_back( { j, dist } );
        all_bonds[j].push_back( { static_cast<int>(i), dist } );
      }
    }

    // per-atom running quaternion sum (starts as each atom's own unit quaternion)
    std::vector<std::array<double,4>> qsum( n_atoms );
    for(size_t i=0;i<n_atoms;i++)
    {
      Mat3d R { orientation_mat3[9*i+0],orientation_mat3[9*i+1],orientation_mat3[9*i+2]
              , orientation_mat3[9*i+3],orientation_mat3[9*i+4],orientation_mat3[9*i+5]
              , orientation_mat3[9*i+6],orientation_mat3[9*i+7],orientation_mat3[9*i+8] };
      gb_matrix_to_quat( R, qsum[i].data() );
    }

    // ---- build the crystalline candidate graph (undirected; both directions of a bond agree on
    // weight since disorientation is a true symmetric distance, see grain_quat_math.h) ----
    std::vector<GraphNode> graph( n_atoms );
    for(size_t i=0;i<n_atoms;i++)
    {
      const int cnt = bond_count[i];
      for(int k=0;k<cnt;k++)
      {
        const size_t slot = i*static_cast<size_t>(max_neighbors) + static_cast<size_t>(k);
        const double deg = bond_disorientation[slot];
        if( deg < 0.0 ) { continue; }
        const uint64_t nid = bond_id[slot];
        auto it = id_to_local.find(nid);
        if( it == id_to_local.end() ) { continue; } // neighbor outside this rank's own local view
        const int j = it->second;
        if( j == static_cast<int>(i) ) { continue; }
        const double w = std::exp( -deg*deg/3.0 );
        graph[i].adj[j] = w;
        graph[j].adj[static_cast<int>(i)] = w;
      }
    }
    for(size_t i=0;i<n_atoms;i++) { recompute_wnode( graph, static_cast<int>(i) ); }

    // ---- Node-Pair-Sampling: reciprocal-nearest-neighbor chain agglomerative clustering ----
    std::unordered_set<int> active;
    for(size_t i=0;i<n_atoms;i++) { if( !graph[i].adj.empty() ) { active.insert( static_cast<int>(i) ); } }

    std::vector<DendrogramEdge> dendrogram;
    std::vector<int> chain;
    while( !active.empty() )
    {
      chain.clear();
      chain.push_back( *active.begin() );
      while( !chain.empty() )
      {
        const int a = chain.back(); chain.pop_back();
        auto [d,b] = nearest_neighbor( graph, a );
        if( b < 0 )
        {
          active.erase(a);
          continue;
        }
        if( !chain.empty() )
        {
          const int c = chain.back(); chain.pop_back();
          if( c == b )
          {
            const int structure_type = static_cast<int>( struct_type[a] );
            const int parent = contract_edge( graph, a, b );
            const int child  = ( parent == a ) ? b : a;
            active.erase(child);
            const double disorientation = merge_orientation( structure_type, qsum[parent].data(), qsum[child].data() );
            dendrogram.push_back( { parent, child, d, disorientation, 0, 0.0 } );
            chain.push_back( parent );
          }
          else { chain.push_back(c); chain.push_back(a); chain.push_back(b); }
        }
        else { chain.push_back(a); chain.push_back(b); }
      }
    }

    if( dendrogram.empty() ) { return; } // nothing crystalline enough to cluster at all

    std::sort( dendrogram.begin(), dendrogram.end(), []( const DendrogramEdge& p, const DendrogramEdge& q ){ return p.distance < q.distance; } );

    // size/merge_size bookkeeping (replay through a fresh union-find), same as OVITO's own pass
    {
      UnionFind uf( n_atoms );
      for( auto& e : dendrogram )
      {
        const long sa = uf.size[ uf.find(e.parent) ];
        const long sb = uf.size[ uf.find(e.child) ];
        e.size = std::min(sa,sb);
        e.merge_size = 2.0 / ( 1.0/static_cast<double>(sa) + 1.0/static_cast<double>(sb) );
        uf.merge( e.parent, e.child );
      }
    }

    const double threshold_log = auto_threshold ? auto_select_threshold( dendrogram ) : manual_threshold_log;
    result.merge_threshold_log = threshold_log;

    // ---- cut the dendrogram: apply every merge whose own log(distance) doesn't exceed threshold ----
    UnionFind uf( n_atoms );
    for( const auto& e : dendrogram )
    {
      if( std::log(e.distance) > threshold_log ) { break; }
      uf.merge( e.parent, e.child );
    }

    // ---- assign contiguous 1-based grain ids to every union-find root whose final cluster size
    // clears min_grain_atom_count, largest-to-smallest (matches OVITO's own grain-table ordering) ----
    std::unordered_map<int,long> root_size;
    for(size_t i=0;i<n_atoms;i++) { ++root_size[ uf.find(static_cast<int>(i)) ]; }
    std::vector<std::pair<int,long>> roots( root_size.begin(), root_size.end() );
    std::sort( roots.begin(), roots.end(), []( const auto& p, const auto& q ){ return p.second > q.second; } );

    // qsum for a merged cluster only ever lives at whichever index the graph-level contract_edge
    // chose as `parent` through the whole chain -- NOT necessarily uf.find()'s own root (union-by-
    // size in UnionFind and union-by-degree in the graph can pick different survivors). Resolve by
    // walking the SAME merge history uf just replayed, tracking which graph-parent index ends up
    // representing each final uf-root.
    std::unordered_map<int,int> uf_root_to_qsum_index;
    {
      UnionFind uf2( n_atoms );
      for( const auto& e : dendrogram )
      {
        if( std::log(e.distance) > threshold_log ) { break; }
        uf2.merge( e.parent, e.child );
        uf_root_to_qsum_index[ uf2.find(e.parent) ] = e.parent; // e.parent is always where qsum accumulated
      }
    }

    result.n_grains = 0;
    for( const auto& [root,size] : roots )
    {
      if( size < min_grain_atom_count ) { continue; }
      ++result.n_grains;
      const long grain_id = result.n_grains;
      const int qsum_idx = uf_root_to_qsum_index.count(root) ? uf_root_to_qsum_index[root] : root;
      const double n = quat_norm( qsum[qsum_idx].data() );
      for(int c=0;c<4;c++) { result.grain_orientation.push_back( n>0. ? qsum[qsum_idx][c]/n : (c==0?1.:0.) ); }
      result.grain_structure_type.push_back( static_cast<int>( struct_type[qsum_idx] ) );
      result.grain_size.push_back( size );
      for(size_t i=0;i<n_atoms;i++) { if( uf.find(static_cast<int>(i)) == root ) { result.atom_grain_id[i] = static_cast<int32_t>(grain_id); } }
    }

    // ---- per-grain random color, fixed seed matching OVITO exactly ----
    {
      std::default_random_engine rng( color_seed );
      std::uniform_real_distribution<double> uniform(0.,1.);
      result.grain_color.resize( result.n_grains );
      for( long g=0; g<result.n_grains; g++ )
      {
        const double hue = uniform(rng), sat = 1.0-0.8*uniform(rng), val = 1.0-0.5*uniform(rng);
        // HSV->RGB, standard 6-sector formula
        const double h6 = hue*6.0; const int sector = static_cast<int>(std::floor(h6)) % 6;
        const double f = h6 - std::floor(h6);
        const double p = val*(1.-sat), q = val*(1.-f*sat), t = val*(1.-(1.-f)*sat);
        double r=val,gc=val,b=val;
        switch( sector<0 ? sector+6 : sector )
        {
          case 0: r=val; gc=t; b=p; break;
          case 1: r=q; gc=val; b=p; break;
          case 2: r=p; gc=val; b=t; break;
          case 3: r=p; gc=q; b=val; break;
          case 4: r=t; gc=p; b=val; break;
          case 5: r=val; gc=p; b=q; break;
        }
        result.grain_color[g] = Vec3d{r,gc,b};
      }
    }

    // ---- orphan atom adoption: multi-source Dijkstra over ALL recorded bonds (not just crystalline
    // candidates), real distance as edge cost, matching OVITO's mergeOrphanAtoms exactly ----
    if( orphan_adoption && result.n_grains > 0 )
    {
      struct PQNode { int cluster; int atom; double length; };
      struct PQCompare { bool operator()( const PQNode& a, const PQNode& b ) const { return a.length > b.length; } };
      std::priority_queue<PQNode, std::vector<PQNode>, PQCompare> pq;

      for(size_t i=0;i<n_atoms;i++)
      {
        if( result.atom_grain_id[i] == 0 ) { continue; } // seed only from already-assigned atoms
        for( const auto& [j,dist] : all_bonds[i] )
        {
          if( result.atom_grain_id[j] != 0 ) { continue; }
          pq.push( { result.atom_grain_id[i], j, dist } );
        }
      }

      while( !pq.empty() )
      {
        const PQNode node = pq.top(); pq.pop();
        if( result.atom_grain_id[ node.atom ] != 0 ) { continue; } // already settled via a cheaper path
        result.atom_grain_id[ node.atom ] = static_cast<int32_t>( node.cluster );
        result.grain_size[ node.cluster - 1 ]++;
        ++result.n_orphans_adopted;
        for( const auto& [j,dist] : all_bonds[ node.atom ] )
        {
          if( result.atom_grain_id[j] != 0 ) { continue; }
          pq.push( { node.cluster, j, node.length + dist } );
        }
      }
    }
  }
}
