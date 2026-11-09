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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/mpi/ghosts_comm_scheme.h>
#include <exanb/mpi/flat_buffer_ghosts.h>

#include <exaStamp/analysis_particle/grains/grain_segmentation_algo.h>
#include <exaStamp/analysis_particle/grains/grain_quat_math.h>

#include <ptm_constants.h>

#include <mpi.h>
#include <vector>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <fstream>
#include <iomanip>

// Grain segmentation, pass (ii): Node-Pair-Sampling clustering of crystalline atoms (edges: bonds to PTM
// template neighbors with the same structure, disorientation under the threshold of pass (i)), then
// orphan adoption. Any number of MPI ranks:
//  1. each rank clusters its owned atoms; node weights include bonds to other ranks' atoms
//  2. collective automatic threshold, local cut into fragments of at most fragment_max_atoms atoms
//  3. fragment graph (bonds within and across ranks) clustered on rank 0, threshold selected again,
//     grains sent back
//  4. orphans take the grain of the closest grain atom along bonds (distributed Dijkstra rounds), then
//     grains are numbered by decreasing final size (ties: smallest atom id)
//  5. per-grain table: atom count, structure, mean orientation (clustered atoms, collective average)
// Not implemented: coherent twin/stacking-fault interfaces.
namespace exaStamp
{
  using namespace exanb;

  struct GrainIdFieldFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const double * const __restrict__ m_atom_grain_id = nullptr; // 0-based, -1 = none
    const Vec3d * const __restrict__ m_grain_color = nullptr;
    inline void operator () ( size_t cell, unsigned int part, double& id_out, Vec3d& color_out ) const
    {
      const double g = m_atom_grain_id[ m_cell_particle_offset[cell] + part ];
      id_out = g;
      color_out = ( g >= 0.0 ) ? m_grain_color[ static_cast<size_t>(g) ] : Vec3d{0.8,0.8,0.8};
    }
  };

  namespace
  {
    static constexpr uint64_t GRAIN_BOND_EMPTY = std::numeric_limits<uint64_t>::max();
    static constexpr uint64_t NO_KEY = std::numeric_limits<uint64_t>::max();

    inline int ptm_template_neighbors( int structure_type )
    {
      switch( structure_type )
      {
        case PTM_MATCH_FCC: return PTM_NUM_NBRS_FCC;
        case PTM_MATCH_HCP: return PTM_NUM_NBRS_HCP;
        case PTM_MATCH_BCC: return PTM_NUM_NBRS_BCC;
        case PTM_MATCH_ICO: return PTM_NUM_NBRS_ICO;
        case PTM_MATCH_SC:  return PTM_NUM_NBRS_SC;
        default: return 0;
      }
    }

    // global fragment id, exact in a double for rank < 2^20
    inline uint64_t fragment_gid( int rank, size_t local_idx ) { return ( uint64_t(rank) << 32 ) | uint64_t(local_idx); }

    struct LocalEdge { int k1, k2; double w; };   // owned compact indices, k1 < k2
    struct RemoteEdge { int k; size_t ghost; uint64_t gid; double w; }; // owned atom k -> ghost of particle gid
    struct LocalBond { int k; double d; };
    struct RemoteBond { size_t ghost; double d; };

    template<class T> MPI_Datatype mpi_type();
    template<> MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_type<uint64_t>() { return MPI_UINT64_T; }

    // variable-size gather on rank 0
    template<class T>
    std::vector<T> gather_on_root( MPI_Comm comm, const std::vector<T>& local, std::vector<int>& counts, std::vector<int>& displs )
    {
      int rank=0, nprocs=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&nprocs);
      const int n = static_cast<int>( local.size() );
      counts.assign( nprocs, 0 ); displs.assign( nprocs, 0 );
      MPI_Gather( &n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm );
      std::vector<T> all;
      if( rank == 0 )
      {
        for(int p=1;p<nprocs;p++) displs[p] = displs[p-1] + counts[p-1];
        all.resize( displs[nprocs-1] + counts[nprocs-1] );
      }
      MPI_Gatherv( local.data(), n, mpi_type<T>(), all.data(), counts.data(), displs.data(), mpi_type<T>(), 0, comm );
      return all;
    }
  }

  template<class GridT>
  class ComputeGrainClusters : public OperatorNode
  {
    std::string m_grains_file_started; // grains_file already started by this run
    ADD_SLOT( MPI_Comm , mpi , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT   , grid , INPUT_OUTPUT );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme , INPUT , REQUIRED , DocString{"Scheme of the last ghost update"} );
    ADD_SLOT( std::string , struct_field , INPUT , std::string("ptm_type") , DocString{"PTM structure-type field (ptm_fields)"} );
    ADD_SLOT( std::string , orient_field , INPUT , std::string("ptm_orientation") , DocString{"PTM orientation field (ptm_fields), for the mean grain orientations. Missing field: identity orientations"} );
    ADD_SLOT( onika::memory::CudaMMVector<uint64_t> , grain_bond_id            , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double>   , grain_bond_distance      , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double>   , grain_bond_disorientation , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<int>      , grain_bond_count         , INPUT , REQUIRED );
    ADD_SLOT( bool   , auto_threshold  , INPUT , true , DocString{"Select the merge threshold automatically, otherwise use manual_threshold_log"} );
    ADD_SLOT( double , manual_threshold_log , INPUT , 0.0 , DocString{"Merge threshold, log of the NPS distance (not degrees)"} );
    ADD_SLOT( long   , fragment_max_atoms , INPUT , 1000 , DocString{"Largest fragment of the per-rank clustering. Smaller: closer to the serial result at rank boundaries; larger: fewer fragments gathered on rank 0"} );
    ADD_SLOT( long   , min_grain_atom_count , INPUT , 100 , DocString{"Smaller clusters are not grains"} );
    ADD_SLOT( bool   , orphan_adoption , INPUT , true , DocString{"Give atoms of no grain the grain of the closest grain atom along bonds"} );
    ADD_SLOT( unsigned int , color_seed , INPUT , 1u , DocString{"Seed of the random grain colors"} );
    ADD_SLOT( std::string , grain_id_field , INPUT , std::string("grain_id") , DocString{"Output grain-id field (0-based, -1 = none)"} );
    ADD_SLOT( std::string , grain_color_field , INPUT , std::string("grain_color") , DocString{"Output RGB grain-color field"} );
    ADD_SLOT( std::string , grains_file , INPUT , std::string("") , DocString{"CSV file with one line per grain and per call (timestep, id, atoms, structure, mean orientation quaternion, color), written by rank 0. Truncated at the first call of the run, appended to afterwards. Empty: no file"} );
    ADD_SLOT( long , timestep , INPUT , 0 , DocString{"Timestep written in the grains_file lines"} );
    ADD_SLOT( long , n_grains , OUTPUT );
    ADD_SLOT( long , n_orphans_adopted , OUTPUT );
    ADD_SLOT( std::vector<long>   , grain_atom_count  , OUTPUT , DocString{"Atoms of each grain, adopted orphans included, indexed by grain id"} );
    ADD_SLOT( std::vector<long>   , grain_structure   , OUTPUT , DocString{"PTM structure type of each grain"} );
    ADD_SLOT( std::vector<double> , grain_orientation , OUTPUT , DocString{"Mean lattice orientation of each grain, quaternion (w,x,y,z) in the fundamental zone, 4 values per grain"} );

  public:
    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      int rank=0, nprocs=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&nprocs);

      if( ! grid->has_allocated_field( field::mk_generic_real( *struct_field ) ) )
      {
        fatal_error() << "compute_grain_clusters: input field '" << *struct_field << "' does not exist (run compute_ptm + ptm_fields first)" << std::endl;
      }

      const size_t n_flat = grid->number_of_particles();
      const size_t n_cells = grid->number_of_cells();
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      const int max_neighbors = static_cast<int>( grain_bond_count->size()>0 ? (grain_bond_id->size()/grain_bond_count->size()) : 0 );
      const double NaN = std::numeric_limits<double>::quiet_NaN();
      const double inf = std::numeric_limits<double>::infinity();

      // ---------- owned atoms (compact index k) and id lookups ----------
      std::vector<size_t> owned;           // k -> flat index
      std::vector<uint64_t> owned_gid;     // k -> particle id
      std::vector<int> owned_type;         // k -> PTM structure type
      std::vector<double> owned_quat;      // k -> orientation quaternion, 4 values
      const bool has_orient = grid->has_allocated_field( field::mk_generic_mat3( *orient_field ) );
      std::unordered_map<uint64_t,int> id_to_owned;
      std::unordered_map<uint64_t,size_t> id_to_ghost; // not owned here -> a ghost row
      {
        auto cells = grid->cells_accessor();
        auto type_acc = grid->field_accessor( field::mk_generic_real( *struct_field ) );
        auto orient_acc = grid->field_accessor( field::mk_generic_mat3( *orient_field ) );
        for(size_t c=0;c<n_cells;c++)
        {
          if( grid->is_ghost_cell(c) ) continue;
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            id_to_owned[ cells[c][field::id][p] ] = static_cast<int>( owned.size() );
            owned.push_back( cell_particle_offset[c] + p );
            owned_gid.push_back( cells[c][field::id][p] );
            owned_type.push_back( static_cast<int>( cells[c][type_acc][p] ) );
            double q[4] = { 1.0, 0.0, 0.0, 0.0 };
            if( has_orient ) { gb_matrix_to_quat( cells[c][orient_acc][p], q ); }
            owned_quat.insert( owned_quat.end(), q, q+4 );
          }
        }
        for(size_t c=0;c<n_cells;c++)
        {
          if( ! grid->is_ghost_cell(c) ) continue;
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            const uint64_t id = cells[c][field::id][p];
            if( id_to_owned.find(id) == id_to_owned.end() ) id_to_ghost.emplace( id, cell_particle_offset[c] + p );
          }
        }
      }
      const size_t n_owned = owned.size();

      // ---------- graph edges and bonds (orphan adoption) ----------
      std::vector<LocalEdge> local_edges;
      std::vector<RemoteEdge> remote_edges;
      std::vector<std::vector<LocalBond>> local_bonds( n_owned );
      std::vector<std::vector<RemoteBond>> remote_bonds( n_owned );
      std::vector<double> node_weight( n_owned, 0.0 );
      std::vector<char> has_candidate( n_owned, 0 );
      {
        std::vector<int> order;
        for(size_t k=0;k<n_owned;k++)
        {
          const size_t base = owned[k] * static_cast<size_t>(max_neighbors);
          const int cnt = (*grain_bond_count)[ owned[k] ];
          order.resize( cnt );
          for(int s=0;s<cnt;s++) order[s] = s;
          std::sort( order.begin(), order.end(), [&](int a, int b){ return (*grain_bond_distance)[base+a] < (*grain_bond_distance)[base+b]; } );
          const int n_template = ptm_template_neighbors( owned_type[k] );
          for(int r=0;r<cnt;r++)
          {
            const size_t slot = base + order[r];
            const uint64_t nid = (*grain_bond_id)[slot];
            if( nid == GRAIN_BOND_EMPTY ) continue;
            const double dist = (*grain_bond_distance)[slot];
            const double deg = (*grain_bond_disorientation)[slot];
            const bool candidate = ( r < n_template && deg >= 0.0 );
            auto it = id_to_owned.find( nid );
            if( it != id_to_owned.end() )
            {
              const int k2 = it->second;
              if( k2 == static_cast<int>(k) ) continue;
              local_bonds[k].push_back( { k2, dist } );
              local_bonds[k2].push_back( { static_cast<int>(k), dist } );
              if( candidate ) local_edges.push_back( { std::min<int>(k,k2), std::max<int>(k,k2), grain_bond_weight(deg) } );
            }
            else
            {
              auto ig = id_to_ghost.find( nid );
              if( ig == id_to_ghost.end() ) continue;
              remote_bonds[k].push_back( { ig->second, dist } );
              if( candidate ) remote_edges.push_back( { static_cast<int>(k), ig->second, nid, grain_bond_weight(deg) } );
            }
          }
        }
        // an edge seen from both atoms counts once
        std::sort( local_edges.begin(), local_edges.end(), [](const LocalEdge& a, const LocalEdge& b){ return a.k1<b.k1 || (a.k1==b.k1 && a.k2<b.k2); } );
        local_edges.erase( std::unique( local_edges.begin(), local_edges.end(), [](const LocalEdge& a, const LocalEdge& b){ return a.k1==b.k1 && a.k2==b.k2; } ), local_edges.end() );
        for( const auto& e : local_edges ) { node_weight[e.k1] += e.w; node_weight[e.k2] += e.w; has_candidate[e.k1] = has_candidate[e.k2] = 1; }
        for( const auto& e : remote_edges ) { node_weight[e.k] += e.w; has_candidate[e.k] = 1; }
      }

      // ---------- 1. local clustering ----------
      std::vector<NPSEdge> nps_edges; nps_edges.reserve( local_edges.size() );
      for( const auto& e : local_edges ) nps_edges.push_back( { e.k1, e.k2, e.w } );
      const std::vector<double> ones( n_owned, 1.0 );
      const std::vector<NPSMerge> dendro1 = nps_dendrogram( node_weight, ones, nps_edges );

      // ---------- 2. threshold, fragments ----------
      double thr1 = *manual_threshold_log;
      if( *auto_threshold )
      {
        std::vector<double> ms(dendro1.size()), ds(dendro1.size());
        for(size_t i=0;i<dendro1.size();i++) { ms[i]=dendro1[i].merge_size; ds[i]=dendro1[i].distance; }
        thr1 = nps_auto_threshold( comm, ms, ds );
      }
      // larger clusters are built at the fragment level, which sees the bonds across rank boundaries
      std::vector<int> uf( n_owned );
      std::vector<long> uf_size( n_owned, 1 );
      for(size_t k=0;k<n_owned;k++) uf[k]=static_cast<int>(k);
      auto find = [&](int a){ while(uf[a]!=a){ uf[a]=uf[uf[a]]; a=uf[a]; } return a; };
      std::vector<double> kept_ms, kept_ds; // applied local merges
      for( const auto& e : dendro1 )
      {
        if( std::log(e.distance) > thr1 ) break;
        const int a=find(e.a), b=find(e.b);
        if( a == b || uf_size[a] + uf_size[b] > *fragment_max_atoms ) continue;
        uf[b]=a; uf_size[a]+=uf_size[b];
        kept_ms.push_back( e.merge_size ); kept_ds.push_back( e.distance );
      }
      std::vector<int> frag_of( n_owned, -1 );
      std::vector<double> frag_weight, frag_size;
      {
        std::unordered_map<int,int> root_to_frag;
        for(size_t k=0;k<n_owned;k++)
        {
          if( !has_candidate[k] ) continue;
          const int r = find( static_cast<int>(k) );
          auto [it,inserted] = root_to_frag.emplace( r, static_cast<int>(frag_weight.size()) );
          if( inserted ) { frag_weight.push_back(0.0); frag_size.push_back(0.0); }
          frag_of[k] = it->second;
          frag_weight[it->second] += node_weight[k];
          frag_size[it->second] += 1.0;
        }
      }

      // fragment ids of ghost atoms
      std::vector<double> flat_buf( n_flat, -1.0 );
      for(size_t k=0;k<n_owned;k++) if( frag_of[k] >= 0 ) flat_buf[ owned[k] ] = static_cast<double>( fragment_gid( rank, frag_of[k] ) );
      {
        GhostFlatBuffer gb { flat_buf.data(), 1, 1, n_flat, GhostFlatBuffer::ValueType::Double };
        flat_buffer_ghost_exchange<true>( comm, *ghost_comm_scheme, cell_particle_offset, gb, 0x6a1 );
      }

      // ---------- 3. fragment graph on rank 0 ----------
      // per edge: fragment a, fragment b, atom pair key (min id, max id) or NO_KEY, weight
      std::vector<uint64_t> edge_ids; std::vector<double> edge_w;
      {
        std::unordered_map<uint64_t,double> intra; // same-rank fragment pairs, summed
        const uint64_t nf = frag_weight.size();
        for( const auto& e : local_edges )
        {
          const int fa = frag_of[e.k1], fb = frag_of[e.k2];
          if( fa == fb ) continue;
          intra[ uint64_t(std::min(fa,fb))*nf + uint64_t(std::max(fa,fb)) ] += e.w;
        }
        for( const auto& [key,w] : intra )
        {
          edge_ids.insert( edge_ids.end(), { fragment_gid(rank,key/nf), fragment_gid(rank,key%nf), NO_KEY, NO_KEY } );
          edge_w.push_back( w );
        }
        // a pair across ranks may come from both sides: deduplicated on rank 0
        for( const auto& e : remote_edges )
        {
          const double fg = flat_buf[ e.ghost ];
          if( fg < 0.0 ) continue;
          const uint64_t a = owned_gid[e.k];
          edge_ids.insert( edge_ids.end(), { fragment_gid(rank,frag_of[e.k]), static_cast<uint64_t>(fg), std::min(a,e.gid), std::max(a,e.gid) } );
          edge_w.push_back( e.w );
        }
      }
      std::vector<int> frag_counts, frag_displs, cnt_tmp, dsp_tmp;
      const std::vector<double> all_frag_weight = gather_on_root( comm, frag_weight, frag_counts, frag_displs );
      const std::vector<double> all_frag_size = gather_on_root( comm, frag_size, cnt_tmp, dsp_tmp );
      const std::vector<uint64_t> all_edge_ids = gather_on_root( comm, edge_ids, cnt_tmp, dsp_tmp );
      const std::vector<double> all_edge_w = gather_on_root( comm, edge_w, cnt_tmp, dsp_tmp );

      std::vector<NPSMerge> dendro2;
      const size_t n_frag_total = all_frag_weight.size(); // 0 except on rank 0
      if( rank == 0 )
      {
        auto to_node = [&]( uint64_t fgid ) { return static_cast<size_t>( frag_displs[ fgid >> 32 ] ) + static_cast<size_t>( fgid & 0xffffffffull ); };
        const size_t ne = all_edge_w.size();
        std::vector<size_t> remote_idx;
        for(size_t i=0;i<ne;i++) if( all_edge_ids[4*i+2] != NO_KEY ) remote_idx.push_back(i);
        std::sort( remote_idx.begin(), remote_idx.end(), [&](size_t a, size_t b){
          const uint64_t* x=&all_edge_ids[4*a]; const uint64_t* y=&all_edge_ids[4*b];
          return x[2]<y[2] || (x[2]==y[2] && x[3]<y[3]); } );
        std::vector<char> use( ne, 1 );
        for(size_t r=1;r<remote_idx.size();r++)
        {
          const uint64_t* x=&all_edge_ids[4*remote_idx[r-1]]; const uint64_t* y=&all_edge_ids[4*remote_idx[r]];
          if( x[2]==y[2] && x[3]==y[3] ) use[ remote_idx[r] ] = 0;
        }
        std::unordered_map<uint64_t,double> pair_w;
        for(size_t i=0;i<ne;i++)
        {
          if( !use[i] ) continue;
          const size_t u = to_node( all_edge_ids[4*i] ), v = to_node( all_edge_ids[4*i+1] );
          if( u == v ) continue;
          pair_w[ uint64_t(std::min(u,v))*n_frag_total + uint64_t(std::max(u,v)) ] += all_edge_w[i];
        }
        std::vector<NPSEdge> fedges; fedges.reserve( pair_w.size() );
        for( const auto& [key,w] : pair_w ) fedges.push_back( { static_cast<int>(key/n_frag_total), static_cast<int>(key%n_frag_total), w } );
        std::sort( fedges.begin(), fedges.end(), [](const NPSEdge& a, const NPSEdge& b){ return a.u<b.u || (a.u==b.u && a.v<b.v); } ); // deterministic
        dendro2 = nps_dendrogram( all_frag_weight, all_frag_size, fedges );
      }

      // final threshold: applied local merges + fragment merges
      double thr2 = *manual_threshold_log;
      if( *auto_threshold )
      {
        for( const auto& e : dendro2 ) { kept_ms.push_back( e.merge_size ); kept_ds.push_back( e.distance ); }
        thr2 = nps_auto_threshold( comm, kept_ms, kept_ds );
      }

      // cut the fragment dendrogram, dissolve small clusters, number grains largest first
      std::vector<long> all_frag_grain; // 1-based, 0 = none
      long ngrains = 0;
      if( rank == 0 )
      {
        std::vector<int> fuf( n_frag_total ); std::vector<double> fsz( all_frag_size );
        for(size_t i=0;i<n_frag_total;i++) fuf[i]=static_cast<int>(i);
        auto ffind = [&](int a){ while(fuf[a]!=a){ fuf[a]=fuf[fuf[a]]; a=fuf[a]; } return a; };
        for( const auto& e : dendro2 )
        {
          if( std::log(e.distance) > thr2 ) break;
          int a=ffind(e.a), b=ffind(e.b);
          if( a==b ) continue;
          if( fsz[a] < fsz[b] ) std::swap(a,b);
          fuf[b]=a; fsz[a]+=fsz[b];
        }
        struct Root { int root; double size; size_t first; };
        std::vector<Root> roots;
        std::unordered_map<int,size_t> root_idx;
        for(size_t i=0;i<n_frag_total;i++)
        {
          const int r = ffind( static_cast<int>(i) );
          if( root_idx.emplace( r, roots.size() ).second ) roots.push_back( { r, fsz[r], i } );
        }
        std::sort( roots.begin(), roots.end(), [](const Root& a, const Root& b){ return a.size>b.size || (a.size==b.size && a.first<b.first); } );
        std::unordered_map<int,long> root_grain;
        for( const auto& r : roots ) { if( r.size >= static_cast<double>(*min_grain_atom_count) ) root_grain[r.root] = ++ngrains; }
        all_frag_grain.assign( n_frag_total, 0 );
        for(size_t i=0;i<n_frag_total;i++) { auto it = root_grain.find( ffind(static_cast<int>(i)) ); if( it != root_grain.end() ) all_frag_grain[i] = it->second; }
      }
      MPI_Bcast( &ngrains, 1, MPI_LONG, 0, comm );
      std::vector<long> frag_grain( frag_weight.size(), 0 );
      MPI_Scatterv( all_frag_grain.data(), frag_counts.data(), frag_displs.data(), MPI_LONG,
                    frag_grain.data(), static_cast<int>(frag_grain.size()), MPI_LONG, 0, comm );

      // ---------- 4. orphan adoption ----------
      std::vector<double> dist( n_owned, inf ), lab( n_owned, 0.0 );
      for(size_t k=0;k<n_owned;k++)
      {
        const long g = ( frag_of[k] >= 0 ) ? frag_grain[ frag_of[k] ] : 0;
        if( g > 0 ) { dist[k] = 0.0; lab[k] = static_cast<double>(g); }
      }
      std::vector<double> lab_before = lab;
      long rounds = 0;
      if( *orphan_adoption && ngrains > 0 )
      {
        std::vector<double> rows( 2*n_flat, NaN );
        using QE = std::pair<double,int>;
        int any_change = 1;
        while( any_change )
        {
          for(size_t k=0;k<n_owned;k++) { rows[2*owned[k]] = dist[k]; rows[2*owned[k]+1] = lab[k]; }
          GhostFlatBuffer gb { rows.data(), 2, 2, n_flat, GhostFlatBuffer::ValueType::Double };
          flat_buffer_ghost_exchange<true>( comm, *ghost_comm_scheme, cell_particle_offset, gb, 0x6a2 );
          bool changed = false;
          std::priority_queue< QE, std::vector<QE>, std::greater<QE> > pq;
          for(size_t k=0;k<n_owned;k++) if( dist[k] < inf ) pq.push( { dist[k], static_cast<int>(k) } );
          for(size_t k=0;k<n_owned;k++)
          {
            for( const auto& rb : remote_bonds[k] )
            {
              const double gd = rows[2*rb.ghost];
              if( !( gd < inf ) ) continue; // unlabeled or NaN
              const double cand = gd + rb.d;
              if( cand < dist[k] ) { dist[k] = cand; lab[k] = rows[2*rb.ghost+1]; pq.push( { cand, static_cast<int>(k) } ); changed = true; }
            }
          }
          while( !pq.empty() )
          {
            const auto [d,k] = pq.top(); pq.pop();
            if( d > dist[k] ) continue;
            for( const auto& lb : local_bonds[k] )
            {
              const double cand = d + lb.d;
              if( cand < dist[lb.k] ) { dist[lb.k] = cand; lab[lb.k] = lab[k]; pq.push( { cand, lb.k } ); changed = true; }
            }
          }
          any_change = changed ? 1 : 0;
          MPI_Allreduce( MPI_IN_PLACE, &any_change, 1, MPI_INT, MPI_MAX, comm );
          ++ rounds;
        }
      }
      long adopted = 0;
      for(size_t k=0;k<n_owned;k++) if( lab_before[k] == 0.0 && lab[k] > 0.0 ) ++adopted;
      MPI_Allreduce( MPI_IN_PLACE, &adopted, 1, MPI_LONG, MPI_SUM, comm );

      // ---------- renumber grains by final size (adopted orphans included), ties by smallest atom id ----------
      // independent of the number of ranks, unlike the cluster sizes before adoption
      std::vector<long> gsize( ngrains, 0 );
      {
        std::vector<unsigned long long> gmin( ngrains, std::numeric_limits<unsigned long long>::max() );
        for(size_t k=0;k<n_owned;k++) if( lab[k] > 0.0 )
        {
          const size_t g = static_cast<size_t>(lab[k]) - 1;
          ++gsize[g]; gmin[g] = std::min<unsigned long long>( gmin[g], owned_gid[k] );
        }
        MPI_Allreduce( MPI_IN_PLACE, gsize.data(), static_cast<int>(ngrains), MPI_LONG, MPI_SUM, comm );
        MPI_Allreduce( MPI_IN_PLACE, gmin.data(), static_cast<int>(ngrains), MPI_UNSIGNED_LONG_LONG, MPI_MIN, comm );
        std::vector<long> order( ngrains );
        for(long g=0; g<ngrains; g++) order[g] = g;
        std::sort( order.begin(), order.end(), [&](long a, long b){ return gsize[a]>gsize[b] || (gsize[a]==gsize[b] && gmin[a]<gmin[b]); } );
        std::vector<double> new_lab( ngrains+1, 0.0 ); // old 1-based label -> new 1-based label
        std::vector<long> new_size( ngrains );
        for(long r=0; r<ngrains; r++) { new_lab[ order[r]+1 ] = static_cast<double>(r+1); new_size[r] = gsize[ order[r] ]; }
        for(size_t k=0;k<n_owned;k++)
        {
          lab[k] = new_lab[ static_cast<size_t>(lab[k]) ];
          lab_before[k] = new_lab[ static_cast<size_t>(lab_before[k]) ];
        }
        gsize.swap( new_size );
      }

      // ---------- output fields, ghosts included ----------
      std::fill( flat_buf.begin(), flat_buf.end(), -1.0 );
      for(size_t k=0;k<n_owned;k++) flat_buf[ owned[k] ] = lab[k] - 1.0;
      {
        GhostFlatBuffer gb { flat_buf.data(), 1, 1, n_flat, GhostFlatBuffer::ValueType::Double };
        flat_buffer_ghost_exchange<true>( comm, *ghost_comm_scheme, cell_particle_offset, gb, 0x6a3 );
      }
      const std::vector<Vec3d> colors = grain_colors( ngrains, *color_seed );
      auto id_acc = grid->field_accessor( field::mk_generic_real( *grain_id_field ) );
      auto color_acc = grid->field_accessor( field::mk_generic_vec3( *grain_color_field ) );
      GrainIdFieldFunctor func = { cell_particle_offset, flat_buf.data(), colors.data() };
      compute_cell_particles( *grid, true, func, onika::make_flat_tuple( id_acc, color_acc ), parallel_execution_context() );

      // ---------- per-grain table ----------
      std::vector<long> mgrain; std::vector<uint64_t> mid; std::vector<int> mtype; std::vector<double> mquat;
      for(size_t k=0;k<n_owned;k++)
      {
        if( !( lab_before[k] > 0.0 ) ) continue; // clustered atoms only, not adopted orphans
        mgrain.push_back( static_cast<long>(lab_before[k]) - 1 );
        mid.push_back( owned_gid[k] );
        mtype.push_back( owned_type[k] );
        mquat.insert( mquat.end(), owned_quat.begin()+4*k, owned_quat.begin()+4*k+4 );
      }
      std::vector<int> gstruct;
      const std::vector<double> gquat = grain_mean_orientations( comm, ngrains, mgrain, mid, mtype, mquat, gstruct );
      if( rank == 0 && ! grains_file->empty() )
      {
        // first call of the run (or new file name): start a new file, append afterwards
        const bool new_file = ( *grains_file != m_grains_file_started );
        std::ofstream out( *grains_file, new_file ? std::ios::trunc : std::ios::app );
        if( ! out ) { fatal_error() << "compute_grain_clusters: cannot write " << *grains_file << std::endl; }
        if( new_file ) { out << "timestep,grain_id,atoms,structure,qw,qx,qy,qz,r,g,b\n"; m_grains_file_started = *grains_file; }
        out << std::setprecision(8);
        for(long g=0; g<ngrains; g++)
        {
          out << *timestep << ',' << g << ',' << gsize[g] << ',' << gstruct[g];
          for(int c=0;c<4;c++) out << ',' << gquat[4*g+c];
          out << ',' << colors[g].x << ',' << colors[g].y << ',' << colors[g].z << '\n';
        }
      }
      grain_atom_count->assign( gsize.begin(), gsize.end() );
      grain_structure->assign( gstruct.begin(), gstruct.end() );
      *grain_orientation = gquat;

      *n_grains = ngrains;
      *n_orphans_adopted = adopted;
      unsigned long long n_atoms = n_owned;
      MPI_Allreduce( MPI_IN_PLACE, &n_atoms, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm );
      lout << "compute_grain_clusters: " << ngrains << " grains found (" << *min_grain_atom_count << "+ atoms each) in "
           << n_atoms << " atoms, " << adopted << " orphan atoms adopted (" << rounds << " rounds), "
           << n_frag_total << " fragments on " << nprocs << " ranks, merge threshold (log) = " << thr2
           << " (fragments cut at " << thr1 << ")" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Grain segmentation, pass (ii): Node-Pair-Sampling clustering of crystalline atoms by
misorientation with an automatic threshold, then orphan adoption. Any number of MPI ranks. Writes a
grain_id field (0-based, -1 = none, ghosts included) and a grain_color field, and a per-grain table
(atom count, structure, mean orientation) in the grain_* output slots and, if grains_file is set, in
a CSV file. Needs a ghost communication scheme (run after a ghost update).

Usage example:

compute_ptm: { rcut: 3.6 ang }
ptm_fields: { struct_field: ptm_type, orient_field: ptm_orientation }
ghost_update_opt: { opt_fields: [ ptm_type, ptm_orientation ] }
compute_grain_bond_misorientation: { rcut: 3.6 ang }
compute_grain_clusters: { min_grain_atom_count: 100, orphan_adoption: true }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_grain_clusters)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_grain_clusters", make_grid_variant_operator< ComputeGrainClusters > );
  }

}
