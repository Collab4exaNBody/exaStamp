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

#include <exaStamp/delaunay/delaunay_tessellation.h>
#include <exaStamp/delaunay/interface_mesh.h>
#include <exaStamp/delaunay/dxa_crystal_path.h>
#include <exaStamp/delaunay/dxa_core_ownership.h>

#include <algorithm>
#include <array>
#include <deque>
#include <unordered_map>

// Extends compute_dxa_circuit_sweep's own per-triangle boundary claims (DXATriangleOwnership,
// derived from that operator's own internal facet_owner) into a full per-tetrahedron "core atom"
// marking, matching OVITO's own markCoreAtoms *result* (which tets/atoms belong to which
// dislocation's own core) without porting OVITO's own MECHANISM (a continuous 3D spatial-query +
// tetrahedron-triangle intersection test run incrementally as each circuit's own cap sweeps
// forward, see DislocationTracer.cpp's own traceSegment/appendLinePoint and
// DislocationAnalysisEngine.cpp's own assignCoreAtomDislocationIDs).
//
// Mechanism here instead: a discrete, combinatorial multi-source breadth-first flood through
// bad-tet-to-bad-tet face adjacency, seeded from every "bad" tet directly behind a triangle
// compute_dxa_circuit_sweep's own circuit already claimed (DXATriangleOwnership -- itself entirely
// unchanged/reused, no new geometry needed there). Every seed's own front expands ONE HOP AT A TIME,
// simultaneously across all dislocations at once (not sequentially, dislocation by dislocation) --
// an UNCLAIMED bad tet is assigned to whichever dislocation's own front reaches it in the FEWEST
// hops; a tet reached by two (or more) different dislocations' own fronts in the exact SAME round
// is a genuine tie, broken deterministically by processing seeds in ascending dislocation_id order
// (so a lower dislocation_id always wins a same-round tie) -- reproducible, not queue-order-
// accidental.
//
// Why this needs to be a genuine multi-source COMPETITION, not N independent single-source flood-
// fills run one after another: the "bad" (defective) tet region is one single CONNECTED 3D blob
// wherever two dislocations' own cores are close enough to physically touch (most obviously at a
// real multi-way junction, but also just two lines passing near each other) -- flooding outward
// from ONE dislocation's own claimed boundary with no competitor present would run straight through
// into a NEIGHBORING dislocation's own territory too, since there is no wall between them in the
// bad-tet interior (the only real boundary that exists at all is the 2D interface mesh between good
// and bad tets, which the circuit sweep already walks -- that boundary says nothing about where one
// dislocation's own INTERIOR territory ends and another's begins). Running every dislocation's own
// front outward AT THE SAME TIME instead gives each contested tet to whichever front is physically
// closer (fewest tet-hops) -- the discrete, purely combinatorial analog of OVITO's own continuous
// swept-cap claiming, without needing to replicate its underlying 3D geometry kernel at all. An
// isolated dislocation far from any other gets its own entire tube uncontested, same as OVITO's own
// result; only the exact seam between two adjacent dislocations' own cores can differ from OVITO's
// own (itself sweep-order-dependent, not objectively "more correct") choice of exactly where to draw
// that internal line.
//
// Output stays entirely LOCAL to this rank, in this operator's own per-rank dislocation_id numbering
// (matching compute_dxa_circuit_sweep's own DXADislocationLines::dislocation_id at this rank, BEFORE
// MPI stitching) -- compute_dxa_mpi_stitch_lines' own new local-to-final remap (see that file's own
// header comment) is what a consumer (dxa_mark_core_atoms) needs to translate this into the FINAL,
// post-stitch dislocation id before actually marking atoms.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXACoreAtoms : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation                 , delaunay_tessellation                  , INPUT , REQUIRED );
    ADD_SLOT( DXAElasticMappingTetClassification   , dxa_elastic_mapping_tet_classification , INPUT , REQUIRED );
    ADD_SLOT( InterfaceMesh                        , interface_mesh                         , INPUT , REQUIRED );
    ADD_SLOT( DXATriangleOwnership                 , dxa_triangle_ownership                 , INPUT , REQUIRED );
    ADD_SLOT( DXACoreTetOwnership                  , dxa_core_tet_ownership                 , OUTPUT );
    ADD_SLOT( long                                 , n_core_tets                            , OUTPUT , DocString{"Number of tetrahedra assigned to some dislocation's own core (out of the total bad-tet count)"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXAElasticMappingTetClassification& cls = *dxa_elastic_mapping_tet_classification;
      const InterfaceMesh& iface = *interface_mesh;
      const DXATriangleOwnership& tri_own = *dxa_triangle_ownership;
      const size_t n_tets = mesh.tetrahedra.size();

      // Bad-tet-to-bad-tet face adjacency: same face-hash technique compute_dxa_elastic_interface_
      // mesh.cpp already uses to find good/bad tet PAIRS sharing a face -- here instead keeping
      // only BAD/BAD pairs (an interior face of the defective region, never part of the interface
      // mesh at all, which only ever records good/bad boundary faces).
      //
      // A bad tet may only participate in this adjacency if it isn't a ghost-halo-limit artifact --
      // same exact criterion compute_dxa_elastic_interface_mesh.cpp already uses when deciding what
      // to mesh (unresolved AND touching ghost territory), NOT a blanket "unresolved" exclusion:
      // this session's own earlier investigation found "unresolved" also fires legitimately deep in
      // the interior near genuine dislocation cores (a real physical signal, not just a ghost-data
      // gap), and blanket-excluding it there caused a real, measured regression (quadrupole
      // dislocation count 9->17 at np=1, where there's no ghost boundary at all to speak of).
      // Missing this ghost-touching-specific exclusion here (i.e. treating ALL "bad" tets as valid
      // adjacency, unresolved or not, ghost-touching or not) let the flood leak through the much
      // WIDER unresolved halo surrounding the true core (44936 unresolved tets on the quadrupole
      // test vs only 712 genuinely-bad-and-resolved ones) instead of staying confined to the real
      // defect region -- caught via the screw-dipole test (2 physically symmetric, equal-length
      // dislocations should get equal core-atom counts; they didn't, because the flood was leaking
      // into the asymmetric ghost-halo-limit noise around each one, not the real, symmetric core).
      static constexpr int face_lv[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };
      struct FaceKeyHash
      {
        inline size_t operator()( const std::array<uint32_t,3>& k ) const noexcept
        {
          size_t h = 1469598103934665603ull;
          for(int i=0;i<3;i++) { h ^= k[i]; h *= 1099511628211ull; }
          return h;
        }
      };
      auto touches_ghost = [&]( size_t t ) -> bool
      {
        if( mesh.vertex_is_owned.empty() ) { return false; }
        for( uint32_t v : mesh.tetrahedra[t] ) { if( !mesh.vertex_is_owned[v] ) { return true; } }
        return false;
      };
      auto is_valid_defect_tet = [&]( size_t t ) -> bool
      {
        if( cls.good[t] != 0.0 ) { return false; }
        if( cls.unresolved[t] && touches_ghost(t) ) { return false; }
        return true;
      };
      std::unordered_map<std::array<uint32_t,3>, std::array<int64_t,2>, FaceKeyHash> face_tets;
      face_tets.reserve( n_tets * 2 );
      for(size_t t=0;t<n_tets;t++)
      {
        if( !is_valid_defect_tet(t) ) { continue; } // only bad tets that aren't a ghost-halo-limit artifact contribute faces to this adjacency
        const auto& tet = mesh.tetrahedra[t];
        for(int f=0; f<4; f++)
        {
          std::array<uint32_t,3> key = { tet[face_lv[f][0]], tet[face_lv[f][1]], tet[face_lv[f][2]] };
          std::sort( key.begin(), key.end() );
          auto it = face_tets.find(key);
          if( it == face_tets.end() ) { face_tets[key] = { static_cast<int64_t>(t), -1 }; }
          else { it->second[1] = static_cast<int64_t>(t); }
        }
      }
      std::vector<std::vector<uint32_t>> bad_adj( n_tets );
      for( const auto& [key, tt] : face_tets )
      {
        if( tt[0] < 0 || tt[1] < 0 ) { continue; } // boundary face (mesh edge or good/bad interface) -- no bad-bad adjacency here
        bad_adj[ static_cast<uint32_t>(tt[0]) ].push_back( static_cast<uint32_t>(tt[1]) );
        bad_adj[ static_cast<uint32_t>(tt[1]) ].push_back( static_cast<uint32_t>(tt[0]) );
      }

      // Seeds: every (bad_tet, dislocation_id) pair a claimed triangle points at, deduplicated and
      // sorted by ascending dislocation_id so same-round ties in the BFS below are always broken
      // the same, reproducible way (lower dislocation_id wins).
      std::vector<std::pair<int32_t,uint32_t>> seeds; // (dislocation_id, tet)
      for(size_t t=0;t<tri_own.triangle_dislocation_id.size();t++)
      {
        const int32_t did = tri_own.triangle_dislocation_id[t];
        if( did < 0 ) { continue; }
        seeds.push_back( { did, iface.bad_tet[t] } );
      }
      std::sort( seeds.begin(), seeds.end() );
      seeds.erase( std::unique( seeds.begin(), seeds.end() ), seeds.end() );

      DXACoreTetOwnership& result = *dxa_core_tet_ownership;
      result.tet_dislocation_id.assign( n_tets, -1 );
      std::deque<uint32_t> queue;
      for( const auto& [did, t] : seeds )
      {
        if( result.tet_dislocation_id[t] != -1 ) { continue; } // already claimed by a lower dislocation_id's own seed
        result.tet_dislocation_id[t] = did;
        queue.push_back(t);
      }
      while( !queue.empty() )
      {
        const uint32_t t = queue.front(); queue.pop_front();
        const int32_t did = result.tet_dislocation_id[t];
        for( uint32_t nb : bad_adj[t] )
        {
          if( result.tet_dislocation_id[nb] != -1 ) { continue; }
          result.tet_dislocation_id[nb] = did;
          queue.push_back(nb);
        }
      }

      long n_claimed = 0;
      for( int32_t d : result.tet_dislocation_id ) { if( d != -1 ) { ++n_claimed; } }
      *n_core_tets = n_claimed;
      lout << "compute_dxa_core_atoms: " << n_claimed << " / " << n_tets << " tetrahedra assigned to a dislocation's own core ("
           << seeds.size() << " distinct boundary seeds)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Extends compute_dxa_circuit_sweep's own per-triangle boundary claims (DXATriangleOwnership) into a
full per-tetrahedron core-atom marking via a multi-source breadth-first flood through bad-tet-to-
bad-tet face adjacency -- see this file's own header comment for the full mechanism and how it
differs from (while matching the physical meaning of) OVITO's own markCoreAtoms. Output
(DXACoreTetOwnership) stays in this operator's own LOCAL, per-rank dislocation_id numbering --
consumed by dxa_mark_core_atoms after compute_dxa_mpi_stitch_lines' own local-to-final remap.

Usage example:

compute_dxa_elastic_mapping_tet_classification: {}
compute_dxa_elastic_interface_mesh: {}
compute_dxa_circuit_sweep: {}
compute_dxa_core_atoms: {}

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_dxa_core_atoms)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_core_atoms", make_simple_operator< ComputeDXACoreAtoms > );
  }

}
