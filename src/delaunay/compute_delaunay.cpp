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

#include <geogram/basic/common.h>
#include <geogram/delaunay/periodic_delaunay_3d.h>

#include <vector>

// Builds a 3D Delaunay tessellation of this rank's owned+ghost particles using Geogram's
// PeriodicDelaunay3d, constructed in NON-periodic mode -- exaStamp's own ghost particles already
// are literal duplicated position copies for both periodic images and cross-rank neighbors (the
// same ghost layer every other operator in this codebase relies on), so Geogram never needs to
// handle periodicity itself. MPI parallelism therefore falls out of exaStamp's existing domain
// decomposition (one independent tessellation per rank); OMP parallelism falls out of
// PeriodicDelaunay3d's own internal multithreading (confirmed OpenMP-backed in this build:
// libgeogram.so links libgomp and calls omp_get_max_threads/etc).
//
// Ownership/trust filtering: a tetrahedron has no fixed radius of influence the way a physical
// interaction cutoff does, so tets near the outer edge of the ghost halo aren't guaranteed
// correct. Rather than track a separate correctness margin, ownership and correctness are handled
// by the SAME single rule: a tetrahedron is kept by this rank only if its centroid falls inside
// one of this rank's own OWNED (non-ghost) cells. This guarantees no gaps and no duplicates once
// every rank's output is combined (every point in space belongs to exactly one rank's owned
// region), and it also discards exactly the tets that could be unreliable, since those can only
// occur out near the fringe of the ghost halo -- far from any owned cell.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDelaunay : public OperatorNode
  {
    ADD_SLOT( GridT               , grid                 , INPUT );
    ADD_SLOT( Domain              , domain               , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation, delaunay_tessellation, OUTPUT );
    ADD_SLOT( bool                , keep_ghost_tets      , INPUT , false , DocString{"false (default): keep only tets whose centroid falls in an owned cell -- the original no-gaps/no-duplicates convention every other consumer (write_delaunay_vtk, the 96000-tet MPI/ghost regression test, ...) relies on. true: ALSO keep tets whose centroid falls in a ghost cell (still discarding genuinely out-of-block ones) -- gives each rank a full local view extending into its own ghost halo, deliberately overlapping with neighboring ranks' own equally-extended view of the same shared region, so a consumer like compute_dxa_circuit_sweep can fully grow a circuit even when its real core sits at or straddles this rank's own owned/ghost seam, rather than being cut off right there. Populates the new DelaunayTessellation::vertex_is_owned so such a consumer can still tell owned from ghost territory and clamp/stitch accordingly."} );

  public:
    inline void execute () override final
    {
      GEO::initialize(); // idempotent (function-local static singleton), safe to call every step

      DelaunayTessellation& result = *delaunay_tessellation;
      result.vertices.clear();
      result.tetrahedra.clear();
      result.vertex_particle_index.clear();
      result.vertex_global_id.clear();
      result.vertex_is_owned.clear();

      const size_t n_cells = grid->number_of_cells();
      auto cells = grid->cells();

      // point_flat_index[pv] is pv's flat particle index in grid->cell_particle_offset_data()'s
      // own convention -- read explicitly here rather than assumed equal to pv, so any consumer
      // of DelaunayTessellation::vertex_particle_index can look up other per-particle data
      // (PTM's flat output buffers, compute_bispectrum's, ...) without relying on iteration-order
      // coincidence.
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      // `points` stays in the SAME raw/reduced (pre-xform) frame as field::rx/ry/rz -- this is the
      // frame domain_periodic_location()/cell partitioning itself operates in (see the classification
      // loop below, and its own comment), so it's the right one for CLASSIFYING which cell a point
      // falls in. It is NOT, in general, real Euclidean space -- Geogram's own Delaunay computation
      // (nearest-neighbor/circumsphere tests) is only geometrically meaningful in real space, and a
      // non-trivial xform (shear, non-uniform scale) can change which triangulation is actually
      // Delaunay between the two frames. `points_real` (== `points` whenever the domain's own xform
      // is the identity, which every test case run so far uses -- this bug was real but silently
      // masked until now) is the one actually fed to Geogram, and the one whose values become
      // DelaunayTessellation::vertices, so every downstream consumer (interface mesh, circuit sweep's
      // own recorded positions/Burgers vectors, compute_dxa_mpi_stitch_lines' own atom-position
      // matching) gets real physical coordinates.
      std::vector<double> points;
      std::vector<double> points_real;
      std::vector<uint32_t> point_flat_index;
      std::vector<uint64_t> point_global_id;
      std::vector<uint8_t> point_is_owned;
      points.reserve( grid->number_of_particles() * 3 );
      points_real.reserve( grid->number_of_particles() * 3 );
      point_flat_index.reserve( grid->number_of_particles() );
      point_global_id.reserve( grid->number_of_particles() );
      point_is_owned.reserve( grid->number_of_particles() );
      const Mat3d xform = domain->xform();
      for(size_t c=0;c<n_cells;c++)
      {
        const size_t np = cells[c].size();
        const uint8_t cell_owned = grid->is_ghost_cell(c) ? 0 : 1;
        for(size_t p=0;p<np;p++)
        {
          const Vec3d r { cells[c][field::rx][p], cells[c][field::ry][p], cells[c][field::rz][p] };
          const Vec3d r_real = xform * r;
          points.push_back( r.x ); points.push_back( r.y ); points.push_back( r.z );
          points_real.push_back( r_real.x ); points_real.push_back( r_real.y ); points_real.push_back( r_real.z );
          point_flat_index.push_back( static_cast<uint32_t>( cell_particle_offset[c] + p ) );
          point_global_id.push_back( cells[c][field::id][p] );
          point_is_owned.push_back( cell_owned );
        }
      }
      const size_t nb_points = points.size() / 3;
      if( nb_points < 4 )
      {
        lout << "compute_delaunay: not enough particles ("<<nb_points<<") to tessellate" << std::endl;
        return;
      }

      GEO::PeriodicDelaunay3d delaunay( false );
      delaunay.set_vertices( static_cast<GEO::index_t>(nb_points), points_real.data() );
      delaunay.compute();

      // classify + compact: keep only tets whose centroid falls in one of this rank's owned cells
      const IJK local_dims = grid->dimension();
      const IJK block_start = grid->block().start;
      const IJK domain_dims = domain->grid_dimension();

      // Ranks whose ghost-inclusive local block straddles a periodic domain edge have a
      // block_start that can be negative (or block_start+local_dims > domain_dims) -- e.g.
      // rank 0 on a periodic axis owns global cells [0,5) and ghost_layers=2 gives block_start=-2.
      // domain_periodic_location() always wraps into the canonical [0,domain_dims) range, which
      // is the WRONG representation to subtract block_start from in that case (it would place
      // the cell far outside this rank's local extent even though it's really a nearby ghost
      // cell, just expressed on "the other side" of the wrap) -- so resolve each axis directly
      // against this rank's own local_dims/block_start instead of going through the domain's
      // canonical wrapped frame.
      auto resolve_local_axis = [&]( ssize_t raw, ssize_t start, ssize_t local_n, ssize_t dom_n, bool periodic ) -> ssize_t
      {
        ssize_t local = raw - start;
        if( periodic )
        {
          while( local < 0 )        { local += dom_n; }
          while( local >= local_n ) { local -= dom_n; }
        }
        return local;
      };

      std::vector<int32_t> vertex_remap( nb_points, -1 );
      size_t n_kept = 0;
      size_t n_oob = 0;
      size_t n_ghost = 0;
      const size_t nb_cells_delaunay = delaunay.nb_cells(); // finite cells only (keeps_infinite()==false by default)
      for(size_t t=0;t<nb_cells_delaunay;t++)
      {
        uint32_t v[4];
        Vec3d centroid { 0., 0., 0. };
        for(int lv=0;lv<4;lv++)
        {
          v[lv] = static_cast<uint32_t>( delaunay.cell_vertex( static_cast<GEO::index_t>(t), lv ) );
          const Vec3d p { points[3*v[lv]+0], points[3*v[lv]+1], points[3*v[lv]+2] };
          centroid += p;
        }
        centroid = centroid * 0.25;

        // raw (unwrapped) global cell coordinates -- same "reduced" position frame as field::rx/ry/rz,
        // no xform needed (cell partitioning operates in that frame, same as domain_periodic_location).
        const Vec3d rel_pos = centroid - domain->origin();
        const IJK raw_ijk = make_ijk( rel_pos / domain->cell_size() );

        IJK local_ijk;
        local_ijk.i = resolve_local_axis( raw_ijk.i, block_start.i, local_dims.i, domain_dims.i, domain->periodic_boundary_x() );
        local_ijk.j = resolve_local_axis( raw_ijk.j, block_start.j, local_dims.j, domain_dims.j, domain->periodic_boundary_y() );
        local_ijk.k = resolve_local_axis( raw_ijk.k, block_start.k, local_dims.k, domain_dims.k, domain->periodic_boundary_z() );

        const bool inside_local_grid =
             local_ijk.i>=0 && local_ijk.i<local_dims.i
          && local_ijk.j>=0 && local_ijk.j<local_dims.j
          && local_ijk.k>=0 && local_ijk.k<local_dims.k;
        if( !inside_local_grid ) { n_oob++; continue; }
        const bool is_ghost_tet = grid->is_ghost_cell(local_ijk);
        if( is_ghost_tet )
        {
          n_ghost++;
          if( !*keep_ghost_tets ) { continue; }
        }

        std::array<uint32_t,4> local_v;
        for(int lv=0;lv<4;lv++)
        {
          if( vertex_remap[v[lv]] < 0 )
          {
            vertex_remap[v[lv]] = static_cast<int32_t>( result.vertices.size() );
            result.vertices.push_back( { points_real[3*v[lv]+0], points_real[3*v[lv]+1], points_real[3*v[lv]+2] } );
            result.vertex_particle_index.push_back( point_flat_index[v[lv]] );
            result.vertex_global_id.push_back( point_global_id[v[lv]] );
            result.vertex_is_owned.push_back( point_is_owned[v[lv]] );
          }
          local_v[lv] = static_cast<uint32_t>( vertex_remap[v[lv]] );
        }
        result.tetrahedra.push_back( local_v );
        n_kept++;
      }

      ldbg << "compute_delaunay: " << nb_points << " points (owned+ghost) -> "
           << nb_cells_delaunay << " tetrahedra, " << n_kept << " kept by this rank ("
           << result.vertices.size() << " vertices), " << n_oob << " out-of-local-grid, "
           << n_ghost << " in ghost cells (" << ( *keep_ghost_tets ? "kept" : "discarded" ) << ")" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Builds a 3D Delaunay tessellation of this rank's owned+ghost particles, using Geogram's
PeriodicDelaunay3d in non-periodic mode (exaStamp's own ghost layer already encodes
periodicity/domain decomposition, so Geogram doesn't need to). By default (keep_ghost_tets: false),
only tetrahedra whose centroid falls in one of this rank's own owned (non-ghost) cells are kept, in
a DelaunayTessellation OUTPUT slot (compacted local vertex/tetrahedra numbering) -- this guarantees
no gaps or duplicates when every rank's output is combined, and discards exactly the tets that could
be unreliable near the ghost fringe. See write_delaunay_vtk to export the result.

With keep_ghost_tets: true, ghost-centroid tets are ALSO kept (only genuinely out-of-block ones are
discarded), giving each rank a full local view extending into its own ghost halo -- deliberately
overlapping with neighboring ranks' own equally-extended view of the same shared region, so a
consumer needing to grow something (a circuit, a region) across what would otherwise be a hard
owned/ghost cutoff has room to do so. DelaunayTessellation::vertex_is_owned distinguishes owned from
ghost vertices in this mode (trivially all-owned in the default mode).

Usage example:

compute_delaunay: {}
write_delaunay_vtk: { filename: "delaunay" }

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_delaunay)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_delaunay", make_grid_variant_operator< ComputeDelaunay > );
  }

}
