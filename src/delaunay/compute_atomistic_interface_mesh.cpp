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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exaStamp/delaunay/delaunay_tessellation.h>
#include <exaStamp/delaunay/dxa_edge_vectors.h>
#include <exaStamp/delaunay/interface_mesh.h>

#include <algorithm>
#include <cmath>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>

// Alternative DXA step (v): builds the interface mesh the way Stukowski's own reference DXA
// implementation actually does it (verified against
// /home/lafourcadep/CODES/ANALYSE/DXA_SOURCE/DXA1.3.6/src/analysis/interfacemesh/InterfaceMesh.cpp
// and src/lattice/LatticeTypeBCC.cpp), which is a genuinely different construction from
// compute_interface_mesh.cpp's own (kept separate/untouched -- both are useful, see that file's
// own header comment and src/delaunay/README.md):
//
//   - compute_interface_mesh.cpp: tessellates ALL atoms with compute_delaunay, classifies every
//     TETRAHEDRON as good/bad, and the mesh is the boundary between good and bad tets.
//   - This operator: mesh VERTICES are non-crystalline atoms themselves (real atom positions, not
//     tet-derived points) -- specifically the ones bordering a crystalline atom. Mesh FACETS are
//     built purely from each *crystalline* atom's own fixed local lattice geometry (which of its
//     neighbor slots are non-crystalline), never from an independent tessellation of the whole
//     system. Reference source, verbatim comment: "For each non-crystalline atom that has at
//     least one crystalline neighbor (the so-called interface atoms) a node is created for the
//     interface mesh." This is why, in OVITO, deleting all BCC atoms with the interface mesh still
//     shown makes the remaining (non-BCC) atoms coincide exactly with the mesh's triangle corners.
//
// BCC only for now (the real dislocation-quadrupole test case is BCC Ta) -- FCC/HCP need their own
// "8 Thompson tetrahedra" template (a different, triangle-based decomposition of the 12-neighbor
// shell) instead of BCC's 6-quad template; same technique, not implemented here.
// ponytail: BCC-only; add FCC/HCP Thompson-tetrahedra templates (also already extracted from
// DXA1.3.6, see src/delaunay/README.md) if/when a non-BCC target_structure is needed here.
//
// Algorithm, per crystalline atom A (vertex_matches_target[A], from compute_dxa_edge_vectors):
//   1. Resolve up to 14 of A's real Delaunay neighbors against DXA1.3.6's own canonical BCC
//      template directions (8 first-shell <111>-type + 6 second-shell <100>-type, in A's own
//      ptm_orientation frame) -- same one-sided nearest-ideal-direction snap technique as
//      compute_dxa_edge_vectors.cpp, just done independently for BOTH endpoints of every edge
//      here (that operator only resolves the lower-indexed endpoint's own view, by design; this
//      operator needs every crystalline atom's own full local neighbor map, not just half of
//      them), reusing its already-deduplicated edge list and vertex_matches_target flags.
//   2. For each of the 6 canonical BCC quads (4 first-shell slots + 1 second-shell slot): if all 4
//      first-shell neighbors are themselves non-crystalline, emit one quad facet (as 2 triangles)
//      connecting them; else, for each adjacent pair of first-shell slots that are BOTH
//      non-crystalline, fan a triangle through the second-shell neighbor (if it's also
//      non-crystalline) -- exactly DXAInterfaceMesh::createBCCMeshFacets's own two branches.
//   3. Triangles are deduplicated (a sorted-vertex-tuple set) and oriented so their normal points
//      toward the generating crystalline atom A (i.e. outward from the non-crystalline region),
//      matching InterfaceMesh's own established sign convention.
//
// ponytail: hole-closing only handles clean simple loops (see the pass itself, below) -- not
// DXA1.3.6's full bounded backtracking search (`closeFacetHoles`/`removeUnnecessaryFacets`/
// `duplicateSharedMeshNodes`/`fixMeshEdges`, a real half-edge-mesh cleanup machinery). Upgrade path
// is that full search if n_interface_open_edges shows real data needs it.
//
// CRITICAL requirement, found empirically: dxa_edge_vectors' struct_field (hence
// vertex_matches_target here) MUST come from compute_cna, not compute_ptm/ptm_fields. On the real
// quadrupole test case, compute_ptm's continuous RMSD fit (even at a seemingly tight
// rmsd_cutoff=0.2) accepts a MUCH wider strain tolerance than a real crystal-structure classifier
// should: it found only 281 non-crystalline (owned) atoms, vs. compute_cna's 1024 and OVITO's own
// real CNA analysis on the identical file: 1026 (see data/regression_new/delaunay/ovitodata/
// output_cna_ovito.xyz) -- effectively an exact match once compute_cna's own signature-counting bug
// was fixed (see src/cna/compute_cna.cu, src/delaunay/README.md). Feeding that correctly-sized
// non-crystalline population into this operator gave 2020 triangles vs. OVITO's 1648 (interface
// mesh, output_dxa_ovito.vtk) -- within 22%, down from measured 91-93 triangles when
// dxa_edge_vectors was built against ptm_type instead (the atom population was simply ~5.5x too
// small for the mesh to have anything to enclose). ptm_orientation is still used for the per-atom
// ideal-direction snap below -- only the crystalline/non-crystalline *classification* needs to be
// compute_cna's, not compute_ptm's; PTM's own fitted orientation remains a reasonable proxy for
// "which way is this atom's own local lattice frame rotated" regardless of which classifier judged
// its structure type (same reasoning ptm_shrink_disorder and the earlier CNA struct_field
// experiments already relied on).
namespace exaStamp
{
  using namespace exanb;

  // DXA1.3.6's own canonical BCC lattice template (src/lattice/LatticeTypeBCC.cpp), transcribed
  // verbatim (up to overall scale, irrelevant here since only directions are ever compared):
  // 8 first-shell <111>-type neighbors (indices 0-7) then 6 second-shell <100>-type (indices 8-13).
  static constexpr double bcc_ideal_raw[14][3] = {
    {-0.5,-0.5,-0.5}, { 0.5,-0.5,-0.5}, { 0.5, 0.5,-0.5}, {-0.5, 0.5,-0.5},
    {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5},
    { 0.0, 0.0,-1.0}, { 0.0, 0.0, 1.0}, { 0.0,-1.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}
  };

  // The 6 quads formed by the nearest neighbors of a BCC atom (same source): each is 4 first-shell
  // slot indices (cyclic order around the square face) plus the second-shell slot whose direction
  // passes through that square's center.
  struct BCCQuad { int idx[4]; int second; };
  static constexpr BCCQuad bcc_quads[6] = {
    { {0,1,2,3}, 8  },
    { {0,4,5,1}, 10 },
    { {1,5,6,2}, 13 },
    { {2,6,7,3}, 11 },
    { {3,7,4,0}, 12 },
    { {7,6,5,4}, 9  },
  };

  template<class GridT>
  class ComputeAtomisticInterfaceMesh : public OperatorNode
  {
    ADD_SLOT( GridT                , grid                   , INPUT , REQUIRED );
    ADD_SLOT( Domain               , domain                 , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation  , INPUT , REQUIRED );
    ADD_SLOT( DXAEdgeVectors       , dxa_edge_vectors       , INPUT , REQUIRED );
    ADD_SLOT( std::string          , target_structure       , INPUT , std::string("BCC") , DocString{"Crystal structure to build the atomistic interface mesh against -- only BCC implemented so far (see file header comment)"} );
    ADD_SLOT( std::string          , orient_field           , INPUT , std::string("ptm_orientation") , DocString{"Name of the per-particle lattice-orientation rotation-tensor grid field (written by ptm_fields)"} );
    ADD_SLOT( double               , angle_tolerance        , INPUT , 40.0 , DocString{"Sanity bound (degrees) for snapping a real neighbor direction onto one of the 14 canonical BCC template directions, same role/default as compute_dxa_edge_vectors' own angle_tolerance"} );
    ADD_SLOT( InterfaceMesh        , interface_mesh         , OUTPUT );
    ADD_SLOT( long                 , n_interface_triangles  , OUTPUT , DocString{"Number of atomistic interface mesh triangles generated"} );
    ADD_SLOT( long                 , n_interface_open_edges , OUTPUT , DocString{"Number of mesh edges used by exactly one triangle instead of two -- either a domain-decomposition cutoff or a real gap this operator's construction left unclosed (see file header comment, no hole-closing pass yet)"} );

  public:
    inline void execute () override final
    {
      if( *target_structure != "BCC" )
      {
        fatal_error() << "compute_atomistic_interface_mesh: only target_structure=BCC is implemented so far, got '"
                       << *target_structure << "'" << std::endl;
      }

      // flatten orient_field into a per-particle array, same convention as
      // compute_dxa_edge_vectors.cpp
      auto cells = grid->cells_accessor();
      auto orient_acc = grid->field_const_accessor( field::mk_generic_mat3( *orient_field ) );
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells_grid = grid->number_of_cells();
      const size_t n_particles = grid->number_of_particles();

      std::vector<Mat3d> flat_orient( n_particles, make_identity_matrix() );
      for(size_t c=0;c<n_cells_grid;c++)
      {
        const size_t np = cells[c].size();
        for(size_t p=0;p<np;p++)
        {
          flat_orient[ cell_particle_offset[c] + p ] = cells[c][orient_acc][p];
        }
      }

      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXAEdgeVectors& ev = *dxa_edge_vectors;
      const size_t n_vertices = mesh.vertices.size();

      // unit-normalize the canonical template directions once
      Vec3d bcc_ideal[14];
      for(int k=0;k<14;k++)
      {
        Vec3d v { bcc_ideal_raw[k][0], bcc_ideal_raw[k][1], bcc_ideal_raw[k][2] };
        bcc_ideal[k] = v / norm(v);
      }

      // per-vertex 14-slot neighbor map: slot_neighbor[14*v+k] = vertex index resolved into A=v's
      // own k-th canonical BCC direction, or -1 if none/unresolved. Filled from both ends of every
      // edge independently (see file header comment -- unlike compute_dxa_edge_vectors, which only
      // resolves the lower-indexed endpoint).
      std::vector<int32_t> slot_neighbor( 14*n_vertices, -1 );
      const double angle_tol = *angle_tolerance;
      const double cos_tol = std::cos( angle_tol * (M_PI/180.0) );

      auto resolve_slot = [&]( uint32_t center, uint32_t other, const Vec3d& dir_hat, const Mat3d& orient )
      {
        const Vec3d d_ideal = transpose(orient) * dir_hat;
        int best_k = -1; double best_dot = -2.0;
        for(int k=0;k<14;k++)
        {
          const double dot = d_ideal.x*bcc_ideal[k].x + d_ideal.y*bcc_ideal[k].y + d_ideal.z*bcc_ideal[k].z;
          if( dot > best_dot ) { best_dot = dot; best_k = k; }
        }
        if( best_k >= 0 && best_dot >= cos_tol ) { slot_neighbor[ 14*size_t(center) + best_k ] = static_cast<int32_t>(other); }
      };

      for(size_t e=0; e<ev.edges.size(); e++)
      {
        const uint32_t vu = ev.edges[e][0];
        const uint32_t vv = ev.edges[e][1];
        const Vec3d d = domain->xform() * ( mesh.vertices[vv] - mesh.vertices[vu] );
        const double dlen = norm(d);
        if( dlen <= 0.0 ) { continue; }
        const Vec3d d_hat = d / dlen;

        if( ev.vertex_matches_target[vu] )
        {
          const uint32_t pu = mesh.vertex_particle_index[vu];
          resolve_slot( vu, vv, d_hat, flat_orient[pu] );
        }
        if( ev.vertex_matches_target[vv] )
        {
          const uint32_t pv = mesh.vertex_particle_index[vv];
          resolve_slot( vv, vu, Vec3d{-d_hat.x,-d_hat.y,-d_hat.z}, flat_orient[pv] );
        }
      }

      InterfaceMesh& result = *interface_mesh;
      result.triangles.clear();
      result.good_tet.clear();
      result.bad_tet.clear();
      result.edge_triangles.clear();

      std::set<std::array<uint32_t,3>> seen_triangles;

      auto emit_triangle = [&]( uint32_t A, uint32_t a, uint32_t b, uint32_t c )
      {
        // orient so the normal points toward A (the generating crystalline atom), i.e. outward
        // from the non-crystalline region -- same sign convention as compute_interface_mesh.cpp.
        const Vec3d pa = mesh.vertices[a], pb = mesh.vertices[b], pc = mesh.vertices[c];
        const Vec3d centroid = (pa+pb+pc) / 3.0;
        const Vec3d normal = cross( pb-pa, pc-pa );
        const Vec3d toward_A = mesh.vertices[A] - centroid;
        uint32_t v0=a, v1=b, v2=c;
        if( normal.x*toward_A.x + normal.y*toward_A.y + normal.z*toward_A.z < 0.0 ) { std::swap(v1,v2); }

        std::array<uint32_t,3> key { v0, v1, v2 };
        std::array<uint32_t,3> sorted_key = key;
        std::sort( sorted_key.begin(), sorted_key.end() );
        if( !seen_triangles.insert( sorted_key ).second ) { return; } // already emitted (e.g. from a symmetric quad configuration)

        const uint32_t tri_index = static_cast<uint32_t>( result.triangles.size() );
        result.triangles.push_back( key );
        result.edge_triangles[ InterfaceMesh::edge_key(v0,v1) ].push_back( tri_index );
        result.edge_triangles[ InterfaceMesh::edge_key(v1,v2) ].push_back( tri_index );
        result.edge_triangles[ InterfaceMesh::edge_key(v2,v0) ].push_back( tri_index );
      };

      for(uint32_t A=0; A<n_vertices; A++)
      {
        if( !ev.vertex_matches_target[A] ) { continue; }
        const int32_t* slots = &slot_neighbor[ 14*size_t(A) ];

        for(const auto& quad : bcc_quads)
        {
          const int32_t second_v = slots[ quad.second ];
          if( second_v < 0 ) { continue; }
          const bool second_is_hole = !ev.vertex_matches_target[second_v];

          int32_t vtx[4]; bool present[4];
          for(int v=0; v<4; v++)
          {
            const int32_t nb = slots[ quad.idx[v] ];
            present[v] = ( nb >= 0 ) && !ev.vertex_matches_target[nb];
            vtx[v] = present[v] ? nb : -1;
          }

          if( present[0] && present[1] && present[2] && present[3] )
          {
            emit_triangle( A, vtx[0], vtx[1], vtx[2] );
            emit_triangle( A, vtx[0], vtx[2], vtx[3] );
          }
          else if( second_is_hole )
          {
            for(int v1=0; v1<4; v1++)
            {
              const int v2 = (v1+1)%4;
              if( present[v1] && present[v2] ) { emit_triangle( A, vtx[v1], vtx[v2], static_cast<uint32_t>(second_v) ); }
            }
          }
        }
      }

      // Hole closing (a bounded analog of DXA1.3.6's own closeFacetHoles/constructFacetRecursive):
      // the per-atom quad construction above only fires where a whole quad-face's worth of
      // neighbors are holes at once, which under-covers a *thin* (1-2 atom radius) dislocation
      // core -- exactly the common case here. Real boundary of an oriented facet is itself a
      // consistently-oriented 1-manifold, so: for every open (used-by-exactly-1-triangle) edge,
      // the closing facet must traverse it in the OPPOSITE direction from however the existing
      // triangle already winds it -- collect that required direction for every open edge, then
      // walk each vertex's unique required outgoing edge until back at the start, and fan-
      // triangulate the loop in its own traversal order (preserves the induced orientation).
      // ponytail: only closes clean simple cycles (every loop vertex has exactly one required
      // outgoing direction) -- a branch point (3+ open edges) or a genuine dangling edge
      // (domain-decomposition cutoff) makes every loop touching it unwalkable and is left open;
      // upgrade path is DXA1.3.6's full bounded backtracking search if that turns out to matter.
      long n_loops_closed = 0;
      {
        std::unordered_map<uint32_t,uint32_t> outgoing; // from -> to, only where unambiguous (see below)
        std::unordered_map<uint32_t,int> out_count;
        for(const auto& tri : result.triangles)
        {
          for(int e=0;e<3;e++)
          {
            const uint32_t a = tri[e], b = tri[(e+1)%3];
            if( result.edge_triangles[ InterfaceMesh::edge_key(a,b) ].size() != 1 ) { continue; }
            // this triangle winds a->b; the closing facet must use b->a
            outgoing[b] = a;
            ++out_count[b];
          }
        }
        // an ambiguous vertex (more than one candidate outgoing edge, i.e. a branch point) can't
        // be walked unambiguously -- drop it from 'outgoing' so any loop reaching it stops there.
        for(const auto& kv : out_count) { if( kv.second != 1 ) { outgoing.erase(kv.first); } }

        std::set<uint32_t> consumed;
        for(const auto& start_kv : outgoing)
        {
          const uint32_t start = start_kv.first;
          if( consumed.count(start) ) { continue; }

          std::vector<uint32_t> loop_verts;
          uint32_t cur = start;
          bool ok = false;
          for(int step=0; step<64; step++)
          {
            auto it = outgoing.find(cur);
            if( it == outgoing.end() || consumed.count(cur) ) { break; }
            loop_verts.push_back(cur);
            cur = it->second;
            if( cur == start ) { ok = true; break; }
          }
          if( ok && loop_verts.size() >= 3 )
          {
            for(uint32_t v : loop_verts) { consumed.insert(v); }
            for(size_t i=1; i+1<loop_verts.size(); i++)
            {
              const std::array<uint32_t,3> tri { loop_verts[0], loop_verts[i], loop_verts[i+1] };
              const uint32_t tri_index = static_cast<uint32_t>( result.triangles.size() );
              result.triangles.push_back( tri );
              result.edge_triangles[ InterfaceMesh::edge_key(tri[0],tri[1]) ].push_back( tri_index );
              result.edge_triangles[ InterfaceMesh::edge_key(tri[1],tri[2]) ].push_back( tri_index );
              result.edge_triangles[ InterfaceMesh::edge_key(tri[2],tri[0]) ].push_back( tri_index );
            }
            ++n_loops_closed;
          }
        }
      }

      long n_open = 0;
      for(const auto& kv : result.edge_triangles) { if( kv.second.size() == 1 ) { ++n_open; } }

      *n_interface_triangles = static_cast<long>( result.triangles.size() );
      *n_interface_open_edges = n_open;
      lout << "compute_atomistic_interface_mesh: " << result.triangles.size() << " interface triangles ("
           << n_loops_closed << " small holes closed), " << n_open << " open edges remaining" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA interface mesh, built the way Stukowski's own reference DXA implementation does it (verified
against DXA1.3.6 source), as an alternative to compute_interface_mesh's tet-classification-boundary
approach: mesh vertices are non-crystalline atoms themselves (real positions), and facets are built
from each crystalline atom's own local BCC lattice template (6 quads of first/second-shell
neighbors) wherever a template slot lands on a non-crystalline neighbor -- see this file's own
header comment for the full rationale and DXA1.3.6 source references. BCC only for now.

Produces an InterfaceMesh, directly usable by write_interface_mesh (same struct as
compute_interface_mesh's own output).

IMPORTANT: dxa_edge_vectors' struct_field must be compute_cna's output, not compute_ptm's -- PTM's
continuous RMSD fit is far more strain-tolerant than a real crystal-structure classifier, and
under-populates the non-crystalline atom set this operator needs (see file header comment for the
measured numbers). ptm_orientation is still needed for per-atom lattice orientation (compute_cna
doesn't produce one).

Usage example:

compute_ptm: { rcut: 5.0 ang }
ptm_fields: {}
compute_cna: { rcut: 5.0 ang }
cna_fields: { struct_field: cna_type }
compute_delaunay: {}
compute_dxa_edge_vectors: { struct_field: cna_type, target_structure: BCC, angle_tolerance: 40.0 }
compute_atomistic_interface_mesh: { target_structure: BCC, angle_tolerance: 40.0 }
write_interface_mesh: { filename: "paraview/interface_atomistic" }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_atomistic_interface_mesh)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_atomistic_interface_mesh", make_grid_variant_operator< ComputeAtomisticInterfaceMesh > );
  }

}
