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
#include <exaStamp/delaunay/dxa_crystal_path.h>
#include <exaStamp/delaunay/interface_mesh.h>

#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>

// DXA pipeline step (v), the real mechanism: builds the tet-boundary interface mesh (same
// triangle-extraction algorithm as compute_interface_mesh.cpp -- one triangle per tet face shared
// by exactly one good and one bad tet) using the new elastic-mapping tet classification
// (compute_dxa_elastic_mapping_tet_classification, stage 3) instead of the old edge-resolution-
// count heuristic. Produces the SAME InterfaceMesh struct compute_dxa_circuit_sweep already
// consumes, so the sweep itself needs no changes at all -- it only ever depended on
// DelaunayTessellation + InterfaceMesh, never on which classifier built the mesh.
//
// Every interface-mesh edge here IS a real Delaunay tessellation edge (unlike the atomistic mesh,
// whose vertices are disordered atoms with no lattice correspondence of their own and needed a
// separate "generating crystalline atom's own template slots" derivation) -- so
// InterfaceMesh::edge_ideal_vector is populated directly from DXACrystalPathEdgeVectors, no
// re-derivation needed.
//
// Known simplification (ponytail): edge_ideal_vector is copied directly in each edge's own v0
// (lower-indexed vertex)'s cluster frame, with no cluster-transition-aware transport when
// compute_dxa_circuit_sweep later sums vectors around a circuit -- correct for every test case so
// far (a single dislocation network within one crystal, even when buildClusters happens to split it
// into several disconnected same-orientation "cluster islands" with no actual misorientation
// between them), but not a general multi-grain Burgers-vector transport. Revisit if a genuine
// polycrystal test case is ever needed here.
namespace exaStamp
{
  using namespace exanb;

  struct FaceKeyHash
  {
    inline size_t operator()( const std::array<uint32_t,3>& k ) const noexcept
    {
      size_t h = 1469598103934665603ull; // FNV-1a
      for(int i=0;i<3;i++) { h ^= k[i]; h *= 1099511628211ull; }
      return h;
    }
  };

  class ComputeDXAElasticInterfaceMesh : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation                 , delaunay_tessellation                  , INPUT , REQUIRED );
    ADD_SLOT( DXAElasticMappingTetClassification   , dxa_elastic_mapping_tet_classification , INPUT , REQUIRED );
    ADD_SLOT( DXACrystalPathEdgeVectors            , dxa_crystal_path_edge_vectors           , INPUT , REQUIRED );
    ADD_SLOT( InterfaceMesh                        , interface_mesh                         , OUTPUT );
    ADD_SLOT( long                                 , n_interface_triangles                  , OUTPUT , DocString{"Number of interface mesh triangles (good/bad tetrahedron face pairs)"} );
    ADD_SLOT( long                                 , n_interface_cutoff_edges                , OUTPUT , DocString{"Number of interface mesh edges shared by only 1 triangle -- domain-decomposition boundary, not a real defect edge"} );

  public:
    inline void execute () override final
    {
      static constexpr int face_lv[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };

      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXAElasticMappingTetClassification& cls = *dxa_elastic_mapping_tet_classification;
      const DXACrystalPathEdgeVectors& ev = *dxa_crystal_path_edge_vectors;
      const size_t n_tets = mesh.tetrahedra.size();

      std::unordered_map<std::array<uint32_t,3>, std::vector<std::pair<uint32_t,int>>, FaceKeyHash> face_tets;
      face_tets.reserve( n_tets * 2 );
      for(uint32_t t=0; t<n_tets; t++)
      {
        const auto& tet = mesh.tetrahedra[t];
        for(int f=0; f<4; f++)
        {
          std::array<uint32_t,3> key = { tet[face_lv[f][0]], tet[face_lv[f][1]], tet[face_lv[f][2]] };
          std::sort( key.begin(), key.end() );
          face_tets[key].emplace_back( t, f );
        }
      }

      InterfaceMesh& result = *interface_mesh;
      result.triangles.clear();
      result.good_tet.clear();
      result.bad_tet.clear();
      result.edge_triangles.clear();
      result.edge_ideal_vector.clear();

      for(const auto& [key, refs] : face_tets)
      {
        if( refs.size() != 2 ) { continue; } // this rank's own tessellation boundary, far side unknown

        const auto [ta, fa] = refs[0];
        const auto [tb, fb] = refs[1];
        const bool good_a = cls.good[ta] != 0.0;
        const bool good_b = cls.good[tb] != 0.0;
        if( good_a == good_b ) { continue; } // both good or both bad -- not an interface

        const uint32_t good_t = good_a ? ta : tb;
        const uint32_t bad_t  = good_a ? tb : ta;
        const int bad_f = good_a ? fb : fa;
        const auto& bad_tet_verts = mesh.tetrahedra[bad_t];

        std::array<uint32_t,3> tri = { bad_tet_verts[face_lv[bad_f][0]], bad_tet_verts[face_lv[bad_f][1]], bad_tet_verts[face_lv[bad_f][2]] };

        // orient consistently: right-hand-rule normal points away from the bad tet's own 4th
        // vertex (outward, into the good region) -- same convention as compute_interface_mesh.cpp.
        const uint32_t opposite_v = bad_tet_verts[bad_f];
        const Vec3d& p0 = mesh.vertices[tri[0]];
        const Vec3d& p1 = mesh.vertices[tri[1]];
        const Vec3d& p2 = mesh.vertices[tri[2]];
        const Vec3d& pd = mesh.vertices[opposite_v];
        if( dot( p1-p0, cross( p2-p0, pd-p0 ) ) > 0.0 ) { std::swap( tri[1], tri[2] ); }

        result.triangles.push_back( tri );
        result.good_tet.push_back( good_t );
        result.bad_tet.push_back( bad_t );
      }

      for(uint32_t ti=0; ti<result.triangles.size(); ti++)
      {
        const auto& tri = result.triangles[ti];
        for(int e=0;e<3;e++)
        {
          const uint32_t va = tri[e];
          const uint32_t vb = tri[(e+1)%3];
          result.edge_triangles[ InterfaceMesh::edge_key( va, vb ) ].push_back(ti);

          // Populate edge_ideal_vector directly from the already-resolved tessellation edge --
          // canonical v0<v1 direction, same storage convention InterfaceMesh itself documents.
          const auto ekey = InterfaceMesh::edge_key( va, vb );
          if( result.edge_ideal_vector.find(ekey) == result.edge_ideal_vector.end() )
          {
            const auto it = ev.edge_index.find( DXACrystalPathEdgeVectors::key(va,vb) );
            if( it != ev.edge_index.end() && ev.resolved[it->second] )
            {
              result.edge_ideal_vector[ekey] = ev.ideal_vector[it->second];
            }
          }
        }
      }

      size_t n_cutoff = 0;
      for(const auto& [key, tris] : result.edge_triangles) { if( tris.size() == 1 ) { ++n_cutoff; } }

      *n_interface_triangles = static_cast<long>( result.triangles.size() );
      *n_interface_cutoff_edges = static_cast<long>( n_cutoff );
      lout << "compute_dxa_elastic_interface_mesh: " << result.triangles.size() << " interface triangles, "
           << n_cutoff << " cutoff edges (domain-decomposition boundary)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline step (v), the real mechanism: builds the tet-boundary interface mesh using the new
elastic-mapping tet classification (Burgers-circuit-closure + Frank-rotation test, compute_dxa_
elastic_mapping_tet_classification) instead of the old edge-resolution-count heuristic. Produces the
same InterfaceMesh struct compute_dxa_circuit_sweep already consumes -- see this file's own header
comment and src/delaunay/README.md "Remaining work" item 0.

Usage example:

compute_delaunay: {}
compute_dxa_lattice_correspondence: { rcut: 6.0 ang }
compute_dxa_lattice_clusters: {}
compute_dxa_crystal_path_edge_vectors: { crystal_path_steps: 4 }
compute_dxa_elastic_mapping_tet_classification: {}
compute_dxa_elastic_interface_mesh: {}
compute_dxa_circuit_sweep: {}

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_dxa_elastic_interface_mesh)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_elastic_interface_mesh", make_simple_operator< ComputeDXAElasticInterfaceMesh > );
  }

}
