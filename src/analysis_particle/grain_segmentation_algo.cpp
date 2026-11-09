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
#include <exaStamp/analysis_particle/grains/grain_segmentation_algo.h>

#include <unordered_map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace exaStamp
{
  double grain_bond_weight( double disorientation_deg )
  {
    // avoids endless loops on ideal crystals
    if( disorientation_deg < 1e-5 ) { disorientation_deg = 0.0; }
    return std::exp( -disorientation_deg * disorientation_deg / 3.0 );
  }

  namespace
  {
    // union-find with atom counts
    struct WeightedUnionFind
    {
      std::vector<int> parent;
      std::vector<double> size;
      explicit WeightedUnionFind( const std::vector<double>& node_size ) : parent(node_size.size()), size(node_size)
      {
        for(size_t i=0;i<parent.size();i++) parent[i]=static_cast<int>(i);
      }
      int find( int a ) { while( parent[a]!=a ) { parent[a]=parent[parent[a]]; a=parent[a]; } return a; }
      void merge( int a, int b )
      {
        a=find(a); b=find(b); if(a==b) return;
        if( size[a] < size[b] ) { std::swap(a,b); }
        parent[b]=a; size[a]+=size[b];
      }
    };

    struct GraphNode
    {
      std::unordered_map<int,double> adj;
      double wnode = 0.0;
    };

    // the node with more edges survives; parallel edges are summed, node weights accumulate and the
    // neighbors' weights do not change
    int contract_edge( std::vector<GraphNode>& graph, int a, int b )
    {
      const int parent = ( graph[a].adj.size() >= graph[b].adj.size() ) ? a : b;
      const int child  = ( parent == a ) ? b : a;
      for( const auto& [v,w] : graph[child].adj )
      {
        if( v == parent ) { continue; }
        double& pw = graph[parent].adj[v];
        pw += w;
        graph[v].adj[parent] = pw;
        graph[v].adj.erase(child);
      }
      graph[parent].adj.erase(child);
      graph[parent].wnode += graph[child].wnode;
      graph[child].adj.clear();
      graph[child].wnode = 0.0;
      return parent;
    }

    // argmin_v wnode[v]/w(a,v) (ties: smallest index), returns wnode[a]*wnode[v]/w(a,v)
    std::pair<double,int> nearest_neighbor( const std::vector<GraphNode>& graph, int a )
    {
      double best_d = std::numeric_limits<double>::infinity(); int best_v = -1;
      for( const auto& [v,w] : graph[a].adj )
      {
        const double d = graph[v].wnode / w;
        if( d < best_d || ( d == best_d && v < best_v ) ) { best_d = d; best_v = v; }
      }
      return { best_d * graph[a].wnode, best_v };
    }

    unsigned long long global_count_le( MPI_Comm comm, const std::vector<double>& values, double t )
    {
      unsigned long long c = 0;
      for( double v : values ) { if( v <= t ) ++c; }
      MPI_Allreduce( MPI_IN_PLACE, &c, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm );
      return c;
    }

    // k-th smallest value (0-based) over all ranks, by bisection on global counts: once no double lies
    // between lo and hi, hi is the answer
    double distributed_kth( MPI_Comm comm, const std::vector<double>& values, unsigned long long k )
    {
      double lo = std::numeric_limits<double>::infinity(), hi = -std::numeric_limits<double>::infinity();
      for( double v : values ) { lo = std::min(lo,v); hi = std::max(hi,v); }
      MPI_Allreduce( MPI_IN_PLACE, &lo, 1, MPI_DOUBLE, MPI_MIN, comm );
      MPI_Allreduce( MPI_IN_PLACE, &hi, 1, MPI_DOUBLE, MPI_MAX, comm );
      if( global_count_le( comm, values, lo ) > k ) { return lo; }
      for(int it=0; it<2000; it++)
      {
        const double mid = lo + 0.5*(hi-lo);
        if( !( mid > lo && mid < hi ) ) { break; }
        if( global_count_le( comm, values, mid ) > k ) { hi = mid; } else { lo = mid; }
      }
      return hi;
    }
  }

  std::vector<NPSMerge> nps_dendrogram( const std::vector<double>& node_weight,
                                        const std::vector<double>& node_size,
                                        const std::vector<NPSEdge>& edges )
  {
    const size_t n = node_weight.size();
    std::vector<GraphNode> graph( n );
    for(size_t i=0;i<n;i++) { graph[i].wnode = node_weight[i]; }
    for( const auto& e : edges )
    {
      if( e.u == e.v ) { continue; }
      graph[e.u].adj[e.v] += e.w;
      graph[e.v].adj[e.u] += e.w;
    }

    // nodes with edges, smallest index first
    std::set<int> active;
    for(size_t i=0;i<n;i++) { if( !graph[i].adj.empty() ) { active.insert( static_cast<int>(i) ); } }

    std::vector<NPSMerge> dendrogram;
    std::vector<int> chain;
    while( !active.empty() )
    {
      chain.clear();
      chain.push_back( *active.begin() );
      while( !chain.empty() )
      {
        const int a = chain.back(); chain.pop_back();
        const auto [d,b] = nearest_neighbor( graph, a );
        if( b < 0 ) { active.erase(a); continue; } // component fully merged
        if( !chain.empty() )
        {
          const int c = chain.back(); chain.pop_back();
          if( c == b )
          {
            const int parent = contract_edge( graph, a, b );
            const int child  = ( parent == a ) ? b : a;
            active.erase(child);
            dendrogram.push_back( { parent, child, d, 0.0 } );
            chain.push_back( parent );
          }
          else { chain.push_back(c); chain.push_back(a); chain.push_back(b); }
        }
        else { chain.push_back(a); chain.push_back(b); }
      }
    }

    std::sort( dendrogram.begin(), dendrogram.end(), []( const NPSMerge& p, const NPSMerge& q ){ return p.distance < q.distance; } );
    WeightedUnionFind uf( node_size );
    for( auto& e : dendrogram )
    {
      const double sa = uf.size[ uf.find(e.a) ];
      const double sb = uf.size[ uf.find(e.b) ];
      e.merge_size = 2.0 / ( 1.0/sa + 1.0/sb );
      uf.merge( e.a, e.b );
    }
    return dendrogram;
  }

  double nps_auto_threshold( MPI_Comm comm, const std::vector<double>& merge_size, const std::vector<double>& distance )
  {
    const size_t n = merge_size.size();
    unsigned long long n_global = n;
    MPI_Allreduce( MPI_IN_PLACE, &n_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm );
    if( n_global == 0 ) { return -std::numeric_limits<double>::infinity(); }

    std::vector<double> x(n), y(n), w(n);
    for(size_t i=0;i<n;i++) { x[i]=std::log(merge_size[i]); y[i]=std::log(distance[i]); w[i]=merge_size[i]; }

    // least absolute deviations by iteratively reweighted least squares
    double gradient = 0.0, intercept = 0.0;
    for(int iter=0; iter<100; iter++)
    {
      double s[5] = {0,0,0,0,0}; // sw, swx, swy, swxx, swxy
      for(size_t i=0;i<n;i++) { s[0]+=w[i]; s[1]+=w[i]*x[i]; s[2]+=w[i]*y[i]; s[3]+=w[i]*x[i]*x[i]; s[4]+=w[i]*x[i]*y[i]; }
      MPI_Allreduce( MPI_IN_PLACE, s, 5, MPI_DOUBLE, MPI_SUM, comm );
      const double denom = s[0]*s[3] - s[1]*s[1];
      if( std::fabs(denom) > 1e-300 )
      {
        gradient = ( s[0]*s[4] - s[1]*s[2] ) / denom;
        intercept = ( s[2] - gradient*s[1] ) / s[0];
      }
      for(size_t i=0;i<n;i++)
      {
        const double r = std::fabs( y[i] - ( gradient*x[i] + intercept ) );
        w[i] = merge_size[i] / std::max( 1e-4, r );
      }
    }

    std::vector<double> abs_res(n);
    for(size_t i=0;i<n;i++) { abs_res[i] = std::fabs( y[i] - (gradient*x[i]+intercept) ); }
    // median (mean of the two middle values for an even count)
    double mad = distributed_kth( comm, abs_res, n_global/2 );
    if( n_global % 2 == 0 ) { mad = 0.5 * ( mad + distributed_kth( comm, abs_res, n_global/2 - 1 ) ); }

    const double cutoff = 1.5;
    double threshold = 0.0;
    for(size_t i=0;i<n;i++)
    {
      const double residual = y[i] - ( gradient*x[i] + intercept );
      if( residual < cutoff*mad ) { threshold = std::max( threshold, y[i] ); }
    }
    MPI_Allreduce( MPI_IN_PLACE, &threshold, 1, MPI_DOUBLE, MPI_MAX, comm );
    return threshold;
  }

  std::vector<Vec3d> grain_colors( long n_grains, unsigned int seed )
  {
    std::default_random_engine rng( seed );
    std::uniform_real_distribution<double> uniform(0.,1.);
    std::vector<Vec3d> colors( std::max(0L,n_grains) );
    for( long g=0; g<n_grains; g++ )
    {
      const double hue = uniform(rng), sat = 1.0-0.8*uniform(rng), val = 1.0-0.5*uniform(rng);
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
      colors[g] = Vec3d{r,gc,b};
    }
    return colors;
  }
}
