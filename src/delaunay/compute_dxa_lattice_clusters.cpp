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

#include <exaStamp/delaunay/dxa_lattice_correspondence.h>

#include <vector>

// Thin grid wrapper around dxa_build_lattice_clusters() (dxa_lattice_clusters_algo.cpp) -- see that
// file's own header comment for the actual buildClusters()/connectClusters() port. This operator
// only flattens real (already ghost-duplicated, unwrapped) atom positions, same convention as
// compute_delaunay.cpp's own point flattening, and hands off to the grid-independent algorithm.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDXALatticeClusters : public OperatorNode
  {
    ADD_SLOT( GridT                     , grid                        , INPUT , REQUIRED );
    ADD_SLOT( Domain                    , domain                      , INPUT , REQUIRED );
    ADD_SLOT( DXALatticeCorrespondence  , dxa_lattice_correspondence  , INPUT_OUTPUT ); // neighbor_atom rows of unclassified atoms get extended in place
    ADD_SLOT( DXALatticeClusters        , dxa_lattice_clusters        , OUTPUT );
    ADD_SLOT( long                      , n_clusters                  , OUTPUT , DocString{"Number of clusters formed (excluding the reserved null cluster 0)"} );
    ADD_SLOT( long                      , n_transitions               , OUTPUT , DocString{"Number of distinct cluster-to-cluster transitions found (connectClusters)"} );

  public:
    inline void execute () override final
    {
      const size_t n_particles = dxa_lattice_correspondence->structure_type.size();

      std::vector<Vec3d> pos( n_particles );
      std::vector<uint64_t> global_id( n_particles );
      {
        auto cells = grid->cells();
        const size_t * const cpo = grid->cell_particle_offset_data();
        const size_t n_cells = grid->number_of_cells();
        for(size_t c=0;c<n_cells;c++)
        {
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++)
          {
            const size_t i = cpo[c] + p;
            pos[i] = domain->xform() * Vec3d{ cells[c][field::rx][p], cells[c][field::ry][p], cells[c][field::rz][p] };
            global_id[i] = cells[c][field::id][p];
          }
        }
      }

      dxa_build_lattice_clusters( *dxa_lattice_correspondence, pos, global_id, *dxa_lattice_clusters );

      *n_clusters = static_cast<long>( dxa_lattice_clusters->graph.clusters.size() ) - 1;
      long n_conn = 0;
      for( const auto& t : dxa_lattice_clusters->graph.transitions ) { if( !t.is_self_transition() && t.distance == 1 ) { ++n_conn; } }
      *n_transitions = n_conn / 2; // each non-self transition is stored in both directions

      lout << "compute_dxa_lattice_clusters: " << *n_clusters << " clusters, " << *n_transitions << " transitions"
           << " (" << dxa_lattice_clusters->n_appends_dropped_diagnostic << " classified neighbors dropped -- unclassified atom's own row was full)" << std::endl;
      for( size_t ci=1; ci<dxa_lattice_clusters->graph.clusters.size(); ci++ )
      {
        lout << "compute_dxa_lattice_clusters:   cluster " << ci << ": " << dxa_lattice_clusters->graph.clusters[ci].atom_count << " atoms" << std::endl;
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA elastic-mapping step: resolves compute_dxa_lattice_correspondence's per-atom labeling ambiguity
into globally-consistent clusters (buildClusters) and finds the transitions between adjacent
clusters (connectClusters) -- port of OVITO's real StructureAnalysis::buildClusters()/
connectClusters(), see dxa_lattice_clusters_algo.cpp's own file header comment.

Usage example:

compute_dxa_lattice_correspondence: { rcut: 5.0 ang }
compute_dxa_lattice_clusters: {}

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_dxa_lattice_clusters)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_lattice_clusters", make_grid_variant_operator< ComputeDXALatticeClusters > );
  }

}
