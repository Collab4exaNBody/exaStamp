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

// Delaunay tessellation (Geogram, non-periodic) of this rank's owned and ghost particles. A tetrahedron
// is kept if its centroid lies in an owned cell, so the ranks' outputs combine without gaps or
// duplicates and the unreliable tetrahedra near the ghost fringe are dropped.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDelaunay : public OperatorNode
  {
    ADD_SLOT( GridT               , grid                 , INPUT );
    ADD_SLOT( Domain              , domain               , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation, delaunay_tessellation, OUTPUT );
    ADD_SLOT( bool                , keep_ghost_tets      , INPUT , false , DocString{"Also keep tetrahedra whose centroid lies in a ghost cell (overlapping with neighbor ranks); see DelaunayTessellation::vertex_is_owned"} );

  public:
    inline void execute () override final
    {
      GEO::initialize(); // idempotent

      DelaunayTessellation& result = *delaunay_tessellation;
      result.vertices.clear();
      result.tetrahedra.clear();
      result.vertex_particle_index.clear();
      result.vertex_global_id.clear();
      result.vertex_is_owned.clear();

      const size_t n_cells = grid->number_of_cells();
      auto cells = grid->cells();

      // point_flat_index[pv]: flat particle index (cell_particle_offset[cell] + particle)
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      // `points`: reduced (pre-xform) frame, used for cell classification; `points_real`: real coordinates,
      // given to Geogram (Delaunay is only meaningful in real space) and output as vertices
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

      // keep tetrahedra whose centroid lies in an owned cell
      const IJK local_dims = grid->dimension();
      const IJK block_start = grid->block().start;
      const IJK domain_dims = domain->grid_dimension();

      // near a periodic edge block_start can be negative: resolve each axis against the local block, not the
      // wrapped domain frame
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
      const size_t nb_cells_delaunay = delaunay.nb_cells(); // finite cells only
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

        // unwrapped global cell coordinates, reduced frame
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

Delaunay tessellation (Geogram) of each rank's owned and ghost particles. Only tetrahedra whose
centroid lies in an owned cell are kept (no gaps or duplicates across ranks). With keep_ghost_tets,
tetrahedra with a ghost-cell centroid are kept too, and DelaunayTessellation::vertex_is_owned tells
owned from ghost vertices. Export with write_delaunay_vtk.

Usage example:

compute_delaunay: {}
write_delaunay_vtk: { filename: "delaunay" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_delaunay)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_delaunay", make_grid_variant_operator< ComputeDelaunay > );
  }

}
