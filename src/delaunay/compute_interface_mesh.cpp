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
#include <exaStamp/delaunay/dxa_edge_vectors.h>
#include <exaStamp/delaunay/interface_mesh.h>

#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>

// DXA pipeline step (v): build the interface mesh -- the 2D surface separating "good" from "bad"
// tetrahedra (compute_dxa_tet_classification). One triangle per tetrahedron face shared by
// exactly one good and one bad tet. Faces are found by hashing every tet's 4 faces (sorted
// vertex-triple key) and grouping -- a face referenced by exactly 2 tets is interior (the normal
// case); referenced by only 1 means the neighboring tet wasn't kept in this rank's own
// DelaunayTessellation (a domain-decomposition/ghost-fringe artifact, same trust concern
// compute_delaunay.cpp already reasons about for tets themselves), so its far side's good/bad
// status isn't actually known -- deliberately excluded rather than guessed at.
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

  class ComputeInterfaceMesh : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation , INPUT , REQUIRED );
    ADD_SLOT( DXATetClassification , dxa_tet_classification , INPUT , REQUIRED );
    ADD_SLOT( InterfaceMesh        , interface_mesh         , OUTPUT );
    ADD_SLOT( long                 , n_interface_triangles  , OUTPUT , DocString{"Number of interface mesh triangles (good/bad tetrahedron face pairs)"} );
    ADD_SLOT( long                 , n_interface_cutoff_edges , OUTPUT , DocString{"Number of interface mesh edges shared by only 1 triangle -- where this rank's own kept-tetrahedron data cuts off the surface (domain-decomposition boundary), not a real edge of a defect"} );

  public:
    inline void execute () override final
    {
      static constexpr int face_lv[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };

      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXATetClassification& cls = *dxa_tet_classification;
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

        // orient consistently: right-hand-rule normal must point away from the bad tet's own 4th
        // vertex (outward, into the good region), for every triangle in the mesh -- otherwise
        // winding is whatever the bad tet's arbitrary vertex order happens to give, and ParaView's
        // backface culling/lighting can make faces look "missing" depending on view angle even
        // though the mesh is topologically closed.
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
        for(int e=0;e<3;e++) { result.edge_triangles[ InterfaceMesh::edge_key( tri[e], tri[(e+1)%3] ) ].push_back(ti); }
      }

      size_t n_cutoff = 0;
      for(const auto& [key, tris] : result.edge_triangles) { if( tris.size() == 1 ) { ++n_cutoff; } }

      *n_interface_triangles = static_cast<long>( result.triangles.size() );
      *n_interface_cutoff_edges = static_cast<long>( n_cutoff );
      lout << "compute_interface_mesh: " << result.triangles.size() << " interface triangles, "
           << n_cutoff << " cutoff edges (domain-decomposition boundary)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline step (v): builds the interface mesh, the 2D surface separating "good" from "bad"
tetrahedra (compute_dxa_tet_classification) -- one triangle per tetrahedron face shared by exactly
one good and one bad tet. A face whose neighboring tet wasn't kept in this rank's own
DelaunayTessellation (domain-decomposition/ghost-fringe boundary) is excluded rather than guessed
at, since its far side's status isn't known. See write_interface_mesh to export the result.

Usage example:

compute_delaunay: {}
compute_dxa_edge_vectors: { target_structure: BCC }
compute_dxa_tet_classification: {}
compute_interface_mesh: {}
write_interface_mesh: { filename: "paraview/interface" }

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_interface_mesh)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_interface_mesh", make_simple_operator< ComputeInterfaceMesh > );
  }

}
