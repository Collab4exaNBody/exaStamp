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
#include <ptm_quat.h>

#include <mpi.h>
#include <vector>
#include <unordered_map>
#include <map>
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
//  1-3. exact mode (up to serial_max_atoms atoms): the candidate-bond graph is gathered and clustered on
//     rank 0, the result does not depend on the number of ranks. Two-level mode (larger systems):
//     1. each rank clusters its owned atoms; node weights include bonds to other ranks' atoms
//     2. collective automatic threshold, local cut into fragments of at most fragment_max_atoms atoms
//     3. fragment graph (bonds within and across ranks) clustered on rank 0, threshold selected again,
//        grains sent back
//  4. coherent interfaces (handle_coherent_interfaces): an HCP grain (stacking fault, twin plane in FCC)
//     joins the adjacent FCC grains it has coherent bonds with; FCC grains on both sides of it are merged
//     only if their orientations match (stacking fault), otherwise (twin) it joins the side with the most
//     coherent bonds and the twin boundary stays a grain boundary
//  5. orphans take the grain of the closest grain atom along bonds (distributed Dijkstra rounds; bonds
//     symmetrized across ranks, lengths in integer micro-angstroms, equal lengths: smallest grain), then
//     grains are numbered by decreasing final size (ties: smallest atom id)
//  6. per-grain table: atom count, structure, mean orientation (clustered atoms, collective average)
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
    // adoption bonds, lengths in integer micro-angstroms (exact sums and ties)
    struct LocalBond { int k; double d; };
    struct RemoteBond { size_t ghost; double d; };
    inline double quantized_length( double d ) { return std::round( d * 1.0e6 ); }

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

    // candidate bond seen from atom a (owned) to atom b
    struct SideEdge { uint64_t a, b; double w; };

    // Exact clustering, independent of the number of ranks: the candidate-bond graph is gathered on rank 0
    // and clustered there as a whole. Returns the 1-based grain of each owned atom (0 = none), numbered by
    // decreasing cluster size (ties: smallest atom id). An edge seen from both atoms takes the weight
    // computed on the side of the smaller id.
    std::vector<long> serial_grain_clustering( MPI_Comm comm, const std::vector<uint64_t>& owned_gid, const std::vector<SideEdge>& side_edges,
                                               bool auto_thr, double manual_thr, long min_size, long& ngrains, double& thr, size_t& n_nodes_out )
    {
      int rank=0; MPI_Comm_rank(comm,&rank);
      std::vector<uint64_t> eids; std::vector<double> ews;
      eids.reserve( 2*side_edges.size() ); ews.reserve( side_edges.size() );
      for( const auto& e : side_edges ) { eids.push_back( e.a ); eids.push_back( e.b ); ews.push_back( e.w ); }
      std::vector<int> c1, d1, id_counts, id_displs;
      const std::vector<uint64_t> all_eids = gather_on_root( comm, eids, c1, d1 );
      const std::vector<double> all_ews = gather_on_root( comm, ews, c1, d1 );
      const std::vector<uint64_t> all_ids = gather_on_root( comm, owned_gid, id_counts, id_displs );

      std::vector<long> all_grain; // per gathered owned atom
      ngrains = 0; thr = 0.0; n_nodes_out = 0;
      if( rank == 0 )
      {
        // canonical edges: (min id, max id), weight from the smaller-id side if it saw the bond
        struct E { uint64_t a, b; int from_b; double w; };
        std::vector<E> edges; edges.reserve( all_ews.size() );
        for(size_t i=0;i<all_ews.size();i++)
        {
          const uint64_t x = all_eids[2*i], y = all_eids[2*i+1];
          if( x == y ) continue;
          edges.push_back( { std::min(x,y), std::max(x,y), x < y ? 0 : 1, all_ews[i] } );
        }
        std::sort( edges.begin(), edges.end(), [](const E& p, const E& q){
          return p.a<q.a || (p.a==q.a && (p.b<q.b || (p.b==q.b && (p.from_b<q.from_b || (p.from_b==q.from_b && p.w<q.w))))); } );
        edges.erase( std::unique( edges.begin(), edges.end(), [](const E& p, const E& q){ return p.a==q.a && p.b==q.b; } ), edges.end() );
        std::vector<uint64_t> nodes;
        for( const auto& e : edges ) { nodes.push_back( e.a ); nodes.push_back( e.b ); }
        std::sort( nodes.begin(), nodes.end() );
        nodes.erase( std::unique( nodes.begin(), nodes.end() ), nodes.end() );
        const size_t nn = nodes.size();
        n_nodes_out = nn;
        auto node_of = [&]( uint64_t id ) { return static_cast<int>( std::lower_bound( nodes.begin(), nodes.end(), id ) - nodes.begin() ); };
        std::vector<double> wnode( nn, 0.0 );
        std::vector<NPSEdge> nps_edges; nps_edges.reserve( edges.size() );
        for( const auto& e : edges )
        {
          const int u = node_of(e.a), v = node_of(e.b);
          wnode[u] += e.w; wnode[v] += e.w;
          nps_edges.push_back( { u, v, e.w } );
        }
        const std::vector<NPSMerge> dendro = nps_dendrogram( wnode, std::vector<double>( nn, 1.0 ), nps_edges );
        thr = manual_thr;
        if( auto_thr )
        {
          std::vector<double> ms( dendro.size() ), ds( dendro.size() );
          for(size_t i=0;i<dendro.size();i++) { ms[i]=dendro[i].merge_size; ds[i]=dendro[i].distance; }
          thr = nps_auto_threshold( MPI_COMM_SELF, ms, ds );
        }
        std::vector<int> uf( nn ); std::vector<long> sz( nn, 1 );
        for(size_t i=0;i<nn;i++) uf[i]=static_cast<int>(i);
        auto find = [&](int a){ while(uf[a]!=a){ uf[a]=uf[uf[a]]; a=uf[a]; } return a; };
        for( const auto& e : dendro )
        {
          if( std::log(e.distance) > thr ) break;
          int a=find(e.a), b=find(e.b);
          if( a==b ) continue;
          if( sz[a] < sz[b] ) std::swap(a,b);
          uf[b]=a; sz[a]+=sz[b];
        }
        // roots by decreasing size, ties: smallest atom id (nodes are sorted by id)
        std::vector<int> first( nn, -1 ); std::vector<int> roots;
        for(size_t i=0;i<nn;i++) { const int r=find(static_cast<int>(i)); if( first[r] < 0 ) { first[r]=static_cast<int>(i); roots.push_back(r); } }
        std::sort( roots.begin(), roots.end(), [&](int a, int b){ return sz[a]>sz[b] || (sz[a]==sz[b] && first[a]<first[b]); } );
        std::vector<long> root_grain( nn, 0 );
        for( int r : roots ) { if( sz[r] >= min_size ) root_grain[r] = ++ngrains; }
        all_grain.assign( all_ids.size(), 0 );
        for(size_t i=0;i<all_ids.size();i++)
        {
          auto it = std::lower_bound( nodes.begin(), nodes.end(), all_ids[i] );
          if( it != nodes.end() && *it == all_ids[i] ) all_grain[i] = root_grain[ find( static_cast<int>( it - nodes.begin() ) ) ];
        }
      }
      MPI_Bcast( &ngrains, 1, MPI_LONG, 0, comm );
      MPI_Bcast( &thr, 1, MPI_DOUBLE, 0, comm );
      std::vector<long> grain( owned_gid.size(), 0 );
      MPI_Scatterv( all_grain.data(), id_counts.data(), id_displs.data(), MPI_LONG,
                    grain.data(), static_cast<int>(grain.size()), MPI_LONG, 0, comm );
      return grain;
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
    ADD_SLOT( long   , serial_max_atoms , INPUT , 2000000 , DocString{"Up to this number of atoms, the whole bond graph is clustered on rank 0 (exact, same result for any number of ranks); above it, two-level clustering (per-rank fragments, then fragment graph on rank 0)"} );
    ADD_SLOT( long   , fragment_max_atoms , INPUT , 1000 , DocString{"Largest fragment of the per-rank clustering. Smaller: closer to the serial result at rank boundaries; larger: fewer fragments gathered on rank 0"} );
    ADD_SLOT( long   , min_grain_atom_count , INPUT , 100 , DocString{"Smaller clusters are not grains"} );
    ADD_SLOT( bool   , orphan_adoption , INPUT , true , DocString{"Give atoms of no grain the grain of the closest grain atom along bonds"} );
    ADD_SLOT( bool   , handle_coherent_interfaces , INPUT , true , DocString{"Merge HCP grains at coherent interfaces (stacking faults, twin planes) into the adjacent FCC grains"} );
    ADD_SLOT( double , misorientation_threshold , INPUT , 4.0 , DocString{"Largest disorientation (degrees) of a coherent FCC-HCP bond, and of two FCC grains merged through a stacking fault"} );
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
      std::vector<int> flat_type( n_flat, PTM_MATCH_NONE ); // owned and ghost rows (coherent interfaces)
      std::vector<double> flat_quat( 4*n_flat, 0.0 );
      std::vector<uint64_t> flat_id( n_flat, GRAIN_BOND_EMPTY );
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
            flat_type[ owned.back() ] = owned_type.back();
            flat_id[ owned.back() ] = owned_gid.back();
            for(int i=0;i<4;i++) flat_quat[ 4*owned.back()+i ] = q[i];
          }
        }
        for(size_t c=0;c<n_cells;c++)
        {
          if( ! grid->is_ghost_cell(c) ) continue;
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            const uint64_t id = cells[c][field::id][p];
            const size_t f = cell_particle_offset[c] + p;
            if( id_to_owned.find(id) == id_to_owned.end() ) id_to_ghost.emplace( id, f );
            flat_type[f] = static_cast<int>( cells[c][type_acc][p] );
            flat_id[f] = id;
            double q[4] = { 1.0, 0.0, 0.0, 0.0 };
            if( has_orient ) { gb_matrix_to_quat( cells[c][orient_acc][p], q ); }
            for(int i=0;i<4;i++) flat_quat[4*f+i] = q[i];
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
      // coherent FCC-HCP bonds: owned HCP atom k, FCC neighbor (owned or ghost row), from the HCP side only
      struct InterfaceBond { int k; size_t nbh; };
      std::vector<SideEdge> side_edges; // every candidate bond, from the owned side (exact clustering)
      std::vector<InterfaceBond> interface_bonds;
      const bool coherent = *handle_coherent_interfaces && has_orient;
      const double max_interface_rad = (*misorientation_threshold) * M_PI / 180.0;
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
            const double deg = (*grain_bond_disorientation)[slot];
            const bool candidate = ( r < n_template && deg >= 0.0 );
            auto it = id_to_owned.find( nid );
            if( coherent && r < n_template && owned_type[k] == PTM_MATCH_HCP )
            {
              size_t f = std::numeric_limits<size_t>::max();
              if( it != id_to_owned.end() ) { f = owned[ it->second ]; }
              else { auto ig = id_to_ghost.find( nid ); if( ig != id_to_ghost.end() ) f = ig->second; }
              if( f != std::numeric_limits<size_t>::max() && flat_type[f] == PTM_MATCH_FCC )
              {
                double qfcc[4], qhcp[4];
                for(int i=0;i<4;i++) { qfcc[i] = flat_quat[4*f+i]; qhcp[i] = owned_quat[4*k+i]; }
                if( ptm::quat_disorientation_hexagonal_to_cubic( qfcc, qhcp ) < max_interface_rad ) interface_bonds.push_back( { static_cast<int>(k), f } );
              }
            }
            if( it != id_to_owned.end() )
            {
              const int k2 = it->second;
              if( k2 == static_cast<int>(k) ) continue;
              if( candidate ) local_edges.push_back( { std::min<int>(k,k2), std::max<int>(k,k2), grain_bond_weight(deg) } );
              if( candidate ) side_edges.push_back( { owned_gid[k], nid, grain_bond_weight(deg) } );
            }
            else
            {
              auto ig = id_to_ghost.find( nid );
              if( ig == id_to_ghost.end() ) continue;
              if( candidate ) side_edges.push_back( { owned_gid[k], nid, grain_bond_weight(deg) } );
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

      // adoption bonds, independent of the decomposition: two atoms are bonded if either lists the other (the
      // bond lists of the ghosts are exchanged), with the length recorded by the atom of smaller id
      if( *orphan_adoption )
      {
        const size_t W = static_cast<size_t>( max_neighbors );
        std::vector<double> nb( 2*W*n_flat, -1.0 ); // per row: neighbor id, length, ...
        for(size_t k=0;k<n_owned;k++)
        {
          const size_t f = owned[k];
          const int cnt = (*grain_bond_count)[f];
          for(int s=0;s<cnt;s++)
          {
            const uint64_t nid = (*grain_bond_id)[f*W+s];
            if( nid == GRAIN_BOND_EMPTY ) continue;
            nb[2*W*f+2*s] = static_cast<double>( nid ); nb[2*W*f+2*s+1] = (*grain_bond_distance)[f*W+s];
          }
        }
        {
          GhostFlatBuffer gb { nb.data(), 2*W, 2*W, n_flat, GhostFlatBuffer::ValueType::Double };
          flat_buffer_ghost_exchange<true>( comm, *ghost_comm_scheme, cell_particle_offset, gb, 0x6a5 );
        }
        // (owned atom, neighbor id, length, length recorded by the smaller id)
        struct AB { int k; uint64_t nid; double d; int from_max; }; // from_max: recorded by the larger id
        std::vector<AB> ab;
        for(size_t f=0;f<n_flat;f++)
        {
          const uint64_t fid = flat_id[f];
          if( fid == GRAIN_BOND_EMPTY ) continue;
          auto io = id_to_owned.find( fid );
          const bool own_row = ( io != id_to_owned.end() ) && owned[io->second] == f;
          if( !own_row )
          {
            // one ghost row per particle, not the periodic images of owned particles
            auto ig = id_to_ghost.find( fid );
            if( ig == id_to_ghost.end() || ig->second != f ) continue;
          }
          for(size_t s=0;s<W;s++)
          {
            const double v = nb[2*W*f+2*s];
            if( v < 0.0 ) continue;
            const uint64_t nid = static_cast<uint64_t>( v );
            if( nid == fid ) continue;
            const double d = nb[2*W*f+2*s+1];
            if( own_row ) ab.push_back( { io->second, nid, d, fid < nid ? 0 : 1 } );   // listed by the owned atom
            auto it = id_to_owned.find( nid );
            if( it != id_to_owned.end() ) ab.push_back( { it->second, fid, d, fid < nid ? 0 : 1 } ); // listed by the other atom
          }
        }
        std::sort( ab.begin(), ab.end(), [](const AB& a, const AB& b){
          return a.k<b.k || (a.k==b.k && (a.nid<b.nid || (a.nid==b.nid && (a.from_max<b.from_max || (a.from_max==b.from_max && a.d<b.d))))); } );
        for(size_t i=0;i<ab.size();i++)
        {
          if( i>0 && ab[i].k==ab[i-1].k && ab[i].nid==ab[i-1].nid ) continue;
          const double d = quantized_length( ab[i].d );
          auto it = id_to_owned.find( ab[i].nid );
          if( it != id_to_owned.end() ) { local_bonds[ab[i].k].push_back( { it->second, d } ); continue; }
          auto ig = id_to_ghost.find( ab[i].nid );
          if( ig != id_to_ghost.end() ) remote_bonds[ab[i].k].push_back( { ig->second, d } );
        }
      }

      unsigned long long n_atoms = n_owned;
      MPI_Allreduce( MPI_IN_PLACE, &n_atoms, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm );
      const bool serial = static_cast<long long>(n_atoms) <= *serial_max_atoms;
      std::vector<double> flat_buf( n_flat, -1.0 );
      std::vector<long> atom_grain( n_owned, 0 ); // 1-based grain of each owned atom after clustering, 0 = none
      long ngrains = 0;
      double thr1 = *manual_threshold_log, thr2 = *manual_threshold_log;
      size_t n_frag_total = 0;
      if( serial )
      {
        atom_grain = serial_grain_clustering( comm, owned_gid, side_edges, *auto_threshold, *manual_threshold_log, *min_grain_atom_count, ngrains, thr2, n_frag_total );
        thr1 = thr2;
      }
      else
      {
      // ---------- 1. local clustering ----------
      std::vector<NPSEdge> nps_edges; nps_edges.reserve( local_edges.size() );
      for( const auto& e : local_edges ) nps_edges.push_back( { e.k1, e.k2, e.w } );
      const std::vector<double> ones( n_owned, 1.0 );
      const std::vector<NPSMerge> dendro1 = nps_dendrogram( node_weight, ones, nps_edges );

      // ---------- 2. threshold, fragments ----------
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
      n_frag_total = all_frag_weight.size(); // 0 except on rank 0
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
      if( *auto_threshold )
      {
        for( const auto& e : dendro2 ) { kept_ms.push_back( e.merge_size ); kept_ds.push_back( e.distance ); }
        thr2 = nps_auto_threshold( comm, kept_ms, kept_ds );
      }

      // cut the fragment dendrogram, dissolve small clusters, number grains largest first
      std::vector<long> all_frag_grain; // 1-based, 0 = none
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
      for(size_t k=0;k<n_owned;k++) atom_grain[k] = ( frag_of[k] >= 0 ) ? frag_grain[ frag_of[k] ] : 0;
      } // two-level clustering

      std::vector<double> dist( n_owned, inf ), lab( n_owned, 0.0 );
      for(size_t k=0;k<n_owned;k++)
      {
        const long g = atom_grain[k];
        if( g > 0 ) { dist[k] = 0.0; lab[k] = static_cast<double>(g); }
      }

      // ---------- 4. coherent interfaces ----------
      long n_interface_grains = 0;
      if( coherent && ngrains > 0 )
      {
        const size_t ng = static_cast<size_t>( ngrains );
        // grain labels of ghost atoms
        std::fill( flat_buf.begin(), flat_buf.end(), 0.0 );
        for(size_t k=0;k<n_owned;k++) flat_buf[ owned[k] ] = lab[k];
        {
          GhostFlatBuffer gb { flat_buf.data(), 1, 1, n_flat, GhostFlatBuffer::ValueType::Double };
          flat_buffer_ghost_exchange<true>( comm, *ghost_comm_scheme, cell_particle_offset, gb, 0x6a4 );
        }
        // grains before the merge: structure, mean orientation, smallest atom id (tie breaks)
        std::vector<long> mgrain; std::vector<uint64_t> mid; std::vector<int> mtype; std::vector<double> mquat;
        std::vector<unsigned long long> gmin( ng, std::numeric_limits<unsigned long long>::max() );
        for(size_t k=0;k<n_owned;k++)
        {
          if( !( lab[k] > 0.0 ) ) continue;
          const long g = static_cast<long>(lab[k]) - 1;
          mgrain.push_back( g ); mid.push_back( owned_gid[k] ); mtype.push_back( owned_type[k] );
          mquat.insert( mquat.end(), owned_quat.begin()+4*k, owned_quat.begin()+4*k+4 );
          gmin[g] = std::min<unsigned long long>( gmin[g], owned_gid[k] );
        }
        MPI_Allreduce( MPI_IN_PLACE, gmin.data(), static_cast<int>(ng), MPI_UNSIGNED_LONG_LONG, MPI_MIN, comm );
        std::vector<int> gstruct0;
        const std::vector<double> gquat0 = grain_mean_orientations( comm, ngrains, mgrain, mid, mtype, mquat, gstruct0 );

        // coherent bonds per (HCP grain, FCC grain) pair
        std::unordered_map<uint64_t,uint64_t> pair_count;
        for( const auto& ib : interface_bonds )
        {
          const long h = static_cast<long>( lab[ib.k] ) - 1, c = static_cast<long>( flat_buf[ib.nbh] ) - 1;
          if( h < 0 || c < 0 || h == c ) continue;
          if( gstruct0[h] != PTM_MATCH_HCP || gstruct0[c] != PTM_MATCH_FCC ) continue;
          ++ pair_count[ uint64_t(h)*ng + uint64_t(c) ];
        }
        std::vector<uint64_t> local_pairs;
        for( const auto& [key,n] : pair_count ) { local_pairs.push_back( key ); local_pairs.push_back( n ); }
        std::vector<int> pc_counts, pc_displs;
        const std::vector<uint64_t> all_pairs = gather_on_root( comm, local_pairs, pc_counts, pc_displs );

        std::vector<long> new_label( ng+1, 0 ); // old 1-based label -> new 1-based label
        long ng_new = 0;
        if( rank == 0 )
        {
          std::map< long, std::map<long,uint64_t> > candidates; // HCP grain -> FCC grain -> coherent bonds
          for(size_t i=0;i+1<all_pairs.size();i+=2) candidates[ all_pairs[i]/ng ][ all_pairs[i]%ng ] += all_pairs[i+1];
          std::vector<long> guf( ng );
          for(size_t g=0;g<ng;g++) guf[g] = static_cast<long>(g);
          auto gfind = [&](long a){ while(guf[a]!=a){ guf[a]=guf[guf[a]]; a=guf[a]; } return a; };
          for( const auto& [h,cands] : candidates )
          {
            // FCC neighbors with matching orientations form one group (both sides of a stacking fault);
            // different groups are twins of each other
            std::vector<long> cg; std::vector<uint64_t> cn;
            for( const auto& [c,n] : cands ) { cg.push_back(c); cn.push_back(n); }
            const size_t m = cg.size();
            std::vector<size_t> grp( m );
            for(size_t i=0;i<m;i++) grp[i] = i;
            for(size_t i=0;i<m;i++) for(size_t j=0;j<i;j++)
            {
              if( grp[i] != i ) break;
              if( gb_disorientation_deg( PTM_MATCH_FCC, gquat0.data()+4*cg[i], gquat0.data()+4*cg[j] ) < *misorientation_threshold ) grp[i] = grp[j];
            }
            // the group with the most coherent bonds (ties: smallest atom id)
            std::map<size_t,std::pair<uint64_t,unsigned long long>> score; // group -> (bonds, smallest atom id)
            for(size_t i=0;i<m;i++)
            {
              auto [it,ins] = score.emplace( grp[i], std::make_pair( uint64_t(0), std::numeric_limits<unsigned long long>::max() ) );
              it->second.first += cn[i];
              it->second.second = std::min( it->second.second, gmin[ cg[i] ] );
            }
            size_t best = score.begin()->first;
            for( const auto& [g,sc] : score )
            {
              const auto& b = score[best];
              if( sc.first > b.first || ( sc.first == b.first && sc.second < b.second ) ) best = g;
            }
            for(size_t i=0;i<m;i++) if( grp[i] == best ) { const long a = gfind(h), b = gfind(cg[i]); if( a != b ) guf[a] = b; }
            ++ n_interface_grains;
          }
          std::vector<long> root_label( ng, 0 );
          for(size_t g=0;g<ng;g++)
          {
            const long r = gfind( static_cast<long>(g) );
            if( root_label[r] == 0 ) root_label[r] = ++ng_new;
            new_label[g+1] = root_label[r];
          }
        }
        MPI_Bcast( &ng_new, 1, MPI_LONG, 0, comm );
        MPI_Bcast( &n_interface_grains, 1, MPI_LONG, 0, comm );
        MPI_Bcast( new_label.data(), static_cast<int>(ng+1), MPI_LONG, 0, comm );
        for(size_t k=0;k<n_owned;k++) lab[k] = static_cast<double>( new_label[ static_cast<size_t>(lab[k]) ] );
        ngrains = ng_new;
      }

      // ---------- 5. orphan adoption ----------
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
              const double gl = rows[2*rb.ghost+1];
              if( cand < dist[k] || ( cand == dist[k] && gl < lab[k] ) ) { dist[k] = cand; lab[k] = gl; pq.push( { cand, static_cast<int>(k) } ); changed = true; }
            }
          }
          while( !pq.empty() )
          {
            const auto [d,k] = pq.top(); pq.pop();
            if( d > dist[k] ) continue;
            for( const auto& lb : local_bonds[k] )
            {
              const double cand = d + lb.d;
              if( cand < dist[lb.k] || ( cand == dist[lb.k] && lab[k] < lab[lb.k] ) ) { dist[lb.k] = cand; lab[lb.k] = lab[k]; pq.push( { cand, lb.k } ); changed = true; }
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
      // structure of a grain: the smallest PTM type of its members (FCC over the HCP atoms of coherent interfaces)
      std::vector<int> gtype( ngrains, std::numeric_limits<int>::max() );
      for(size_t k=0;k<n_owned;k++) if( lab_before[k] > 0.0 )
      {
        int& t = gtype[ static_cast<size_t>(lab_before[k]) - 1 ];
        t = std::min( t, owned_type[k] );
      }
      MPI_Allreduce( MPI_IN_PLACE, gtype.data(), static_cast<int>(ngrains), MPI_INT, MPI_MIN, comm );
      std::vector<long> mgrain; std::vector<uint64_t> mid; std::vector<int> mtype; std::vector<double> mquat;
      for(size_t k=0;k<n_owned;k++)
      {
        if( !( lab_before[k] > 0.0 ) ) continue; // clustered atoms only, not adopted orphans
        if( owned_type[k] != gtype[ static_cast<size_t>(lab_before[k]) - 1 ] ) continue;
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
      lout << "compute_grain_clusters: " << ngrains << " grains found (" << *min_grain_atom_count << "+ atoms each) in "
           << n_atoms << " atoms, " << n_interface_grains << " coherent-interface clusters merged, "
           << adopted << " orphan atoms adopted (" << rounds << " rounds), "
           << ( serial ? "exact clustering of " : "two-level clustering, " ) << n_frag_total << ( serial ? " atoms" : " fragments" )
           << " on " << nprocs << " ranks, merge threshold (log) = " << thr2;
      if( !serial ) lout << " (fragments cut at " << thr1 << ")";
      lout << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Grain segmentation, pass (ii): Node-Pair-Sampling clustering of crystalline atoms by
misorientation with an automatic threshold, then orphan adoption. Any number of MPI ranks. Writes a
grain_id field (0-based, -1 = none, ghosts included) and a grain_color field, and a per-grain table
(atom count, structure, mean orientation) in the grain_* output slots and, if grains_file is set, in
a CSV file. Needs a ghost communication scheme (run after a ghost update).

With handle_coherent_interfaces (default), HCP grains with coherent bonds to FCC grains (stacking
faults, twin planes) are merged into them: FCC grains on both sides of a stacking fault become one
grain, a twin plane joins one of the two twins (the side with more coherent bonds) and the twin
boundary stays a grain boundary. HCP crystals are left as they are (hexagonal into cubic only).

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
