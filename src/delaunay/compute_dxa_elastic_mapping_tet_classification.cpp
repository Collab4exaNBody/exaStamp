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
#include <exaStamp/delaunay/dxa_crystal_path.h>

// DXA pipeline step (iv), the real mechanism: applies dxa_is_elastic_mapping_compatible() (see
// dxa_elastic_mapping_compatible.cpp) to every tetrahedron of the Delaunay tessellation -- port of
// how OVITO's real InterfaceMesh::createMesh() decides good/bad (via
// ElasticMapping::isElasticMappingCompatible), replacing the old compute_dxa_tet_classification's
// edge-resolution-count heuristic. Not yet wired into compute_interface_mesh (stage 4) -- see
// src/delaunay/README.md "Remaining work" item 0.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXAElasticMappingTetClassification : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation           , delaunay_tessellation          , INPUT , REQUIRED );
    ADD_SLOT( DXACrystalPathEdgeVectors      , dxa_crystal_path_edge_vectors  , INPUT , REQUIRED );
    ADD_SLOT( DXALatticeClusters             , dxa_lattice_clusters           , INPUT_OUTPUT , REQUIRED ); // graph may gain cached transitions
    ADD_SLOT( DXAElasticMappingTetClassification , dxa_elastic_mapping_tet_classification , OUTPUT );
    ADD_SLOT( long                           , n_tets_good                    , OUTPUT , DocString{"Number of tetrahedra classified good (out of delaunay_tessellation->tetrahedra.size())"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXACrystalPathEdgeVectors& ev = *dxa_crystal_path_edge_vectors;
      DXALatticeClusters& clusters = *dxa_lattice_clusters;
      const size_t n_tets = mesh.tetrahedra.size();

      DXAElasticMappingTetClassification& result = *dxa_elastic_mapping_tet_classification;
      result.good.assign( n_tets, 0.0 );

      static constexpr int edge_lv[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };

      size_t n_good = 0;
      long n_bad_missing_edge = 0, n_bad_all_resolved = 0;
      for(size_t t=0;t<n_tets;t++)
      {
        const bool ok = dxa_is_elastic_mapping_compatible( ev, clusters, mesh.tetrahedra[t] );
        if( ok )
        {
          result.good[t] = 1.0;
          ++n_good;
        }
        else
        {
          bool all_resolved = true;
          const auto& tet = mesh.tetrahedra[t];
          for(int i=0;i<6 && all_resolved;i++)
          {
            const auto it = ev.edge_index.find( DXACrystalPathEdgeVectors::key( tet[edge_lv[i][0]], tet[edge_lv[i][1]] ) );
            all_resolved = ( it != ev.edge_index.end() ) && ev.resolved[it->second];
          }
          if( all_resolved ) { ++n_bad_all_resolved; } else { ++n_bad_missing_edge; }
        }
      }

      *n_tets_good = static_cast<long>( n_good );
      lout << "compute_dxa_elastic_mapping_tet_classification: " << n_good << " / " << n_tets << " tetrahedra classified good"
           << " (bad breakdown: " << n_bad_missing_edge << " missing an edge, " << n_bad_all_resolved
           << " all 6 edges resolved but Burgers/Frank test itself failed)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline step (iv), the real mechanism: classifies each Delaunay tetrahedron good/bad via a
genuine per-tetrahedron Burgers-circuit-closure + Frank-rotation test (dxa_is_elastic_mapping_
compatible, port of OVITO's ElasticMapping::isElasticMappingCompatible), not a per-vertex/edge-count
heuristic. See this file's own header comment and src/delaunay/README.md "Remaining work" item 0.

Usage example:

compute_delaunay: {}
compute_dxa_lattice_correspondence: { rcut: 6.0 ang }
compute_dxa_lattice_clusters: {}
compute_dxa_crystal_path_edge_vectors: { crystal_path_steps: 4 }
compute_dxa_elastic_mapping_tet_classification: {}

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_dxa_elastic_mapping_tet_classification)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_elastic_mapping_tet_classification", make_simple_operator< ComputeDXAElasticMappingTetClassification > );
  }

}
