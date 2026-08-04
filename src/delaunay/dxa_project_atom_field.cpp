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

#include <exaStamp/delaunay/dxa_dislocation_lines.h>

#include <mpi.h>
#include <string>
#include <vector>
#include <cmath>

// User-requested operator: projects an arbitrary per-atom scalar field (any exanb generic_real grid
// field -- CNA/PTM type, a cluster id, a per-atom stress/strain component, ...) onto every
// dislocation-line node, for ParaView visualization (colour the extracted dislocation lines by a
// real physical/structural quantity instead of only Burgers vector or dislocation id).
//
// DXADislocationLines only has real content on rank 0 after compute_dxa_mpi_stitch_lines (every
// other rank's own copy is cleared), while atoms are distributed across every rank -- so this
// gathers every rank's own OWNED atom (position, chosen field value) to rank 0 first (same
// MPI_Gatherv pattern write_dxa_ca_file/write_ovito_interface_mesh/compute_dxa_mpi_stitch_lines
// already use), then does the actual per-node projection entirely on rank 0. For each line node,
// averages every gathered atom within `cutoff` (periodic-image-aware, same minimum-image wrap
// convention as compute_dxa_mpi_stitch_lines' own wrap_to_reference); if literally no atom falls
// within `cutoff` of a given node (shouldn't normally happen for a reasonable cutoff, but not
// assumed impossible), falls back to the single nearest atom regardless of distance, so a node is
// never left without some real, physically meaningful value.
//
// ponytail: brute-force O(n_atoms * n_nodes) nearest/within-cutoff search, no spatial index -- fine
// at the scale of this pipeline's own test systems (10^5 atoms, 10^2-10^3 line nodes, a one-shot
// analysis step, not a hot loop). Add a spatial index (e.g. a uniform grid bucket, same idea as
// DelaunayTessellationSpatialQuery in OVITO's own source) if this ever needs to scale to a much
// larger system.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class DXAProjectAtomField : public OperatorNode
  {
    ADD_SLOT( MPI_Comm               , mpi                        , INPUT , REQUIRED );
    ADD_SLOT( GridT                  , grid                       , INPUT , REQUIRED );
    ADD_SLOT( Domain                 , domain                     , INPUT , REQUIRED , DocString{"Only used for periodic-image-aware distance/averaging near a domain boundary"} );
    ADD_SLOT( DXADislocationLines    , dxa_dislocation_lines      , INPUT , REQUIRED );
    ADD_SLOT( std::string            , atom_field                 , INPUT , REQUIRED , DocString{"Name of the per-particle field (field::mk_generic_real, e.g. a field written by dxa_lattice_cluster_fields/dxa_mark_core_atoms/ptm_fields/cna_fields, or any other generic_real field) to project onto dislocation-line nodes"} );
    ADD_SLOT( double                 , cutoff                     , INPUT , 5.0 , DocString{"Max distance (same units as line_positions) an atom may be from a line node to contribute to that node's own averaged value -- kept generous relative to ordinary lattice spacing so every node normally finds several nearby atoms"} );
    ADD_SLOT( DXALineNodeFieldValues , dxa_line_node_field_values , OUTPUT );

  public:
    inline void execute () override final
    {
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      // Flatten this rank's own owned atoms: (x,y,z,field_value) per atom.
      std::vector<double> local_buf;
      {
        auto cells = grid->cells_accessor();
        auto field_acc = grid->field_const_accessor( field::mk_generic_real( *atom_field ) );
        const Mat3d xform = domain->xform();
        const size_t n_cells = grid->number_of_cells();
        for(size_t c=0;c<n_cells;c++)
        {
          if( grid->is_ghost_cell(c) ) { continue; }
          const size_t np_cell = cells[c].size();
          for(size_t p=0;p<np_cell;p++)
          {
            const Vec3d real_pos = xform * Vec3d{ cells[c][field::rx][p], cells[c][field::ry][p], cells[c][field::rz][p] };
            local_buf.push_back( real_pos.x ); local_buf.push_back( real_pos.y ); local_buf.push_back( real_pos.z );
            local_buf.push_back( cells[c][field_acc][p] );
          }
        }
      }
      const int local_count = static_cast<int>( local_buf.size() );

      std::vector<int> recvcounts, displs;
      if( rank == 0 ) { recvcounts.resize(np); }
      MPI_Gather( &local_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, *mpi );

      std::vector<double> all_buf;
      if( rank == 0 )
      {
        displs.resize(np);
        int off = 0;
        for(int r=0;r<np;r++) { displs[r] = off; off += recvcounts[r]; }
        all_buf.resize(off);
      }
      MPI_Gatherv( local_buf.data(), local_count, MPI_DOUBLE,
                   rank==0 ? all_buf.data() : nullptr,
                   rank==0 ? recvcounts.data() : nullptr,
                   rank==0 ? displs.data() : nullptr,
                   MPI_DOUBLE, 0, *mpi );

      DXALineNodeFieldValues& result = *dxa_line_node_field_values;
      result.field_name = *atom_field;
      result.value.clear();
      if( rank != 0 ) { return; } // only rank 0 has the real DXADislocationLines geometry to project onto

      const size_t n_atoms = all_buf.size() / 4;
      const Mat3d xform = domain->xform();
      const Vec3d reduced_size = domain->bounds_size();
      const Vec3d box_size { norm( xform * Vec3d{reduced_size.x,0.,0.} ),
                             norm( xform * Vec3d{0.,reduced_size.y,0.} ),
                             norm( xform * Vec3d{0.,0.,reduced_size.z} ) };
      const bool periodic_x = domain->periodic_boundary_x();
      const bool periodic_y = domain->periodic_boundary_y();
      const bool periodic_z = domain->periodic_boundary_z();
      auto wrap_axis = [&]( double d, double box, bool periodic ) -> double
      {
        if( !periodic || box <= 0. ) { return d; }
        while( d >  0.5*box ) { d -= box; }
        while( d < -0.5*box ) { d += box; }
        return d;
      };
      auto min_image_dist = [&]( const Vec3d& p, const Vec3d& q ) -> double
      {
        const Vec3d d{ wrap_axis( p.x-q.x, box_size.x, periodic_x ),
                       wrap_axis( p.y-q.y, box_size.y, periodic_y ),
                       wrap_axis( p.z-q.z, box_size.z, periodic_z ) };
        return norm(d);
      };

      const double cut = *cutoff;
      const DXADislocationLines& dl = *dxa_dislocation_lines;
      result.value.resize( dl.line_positions.size() );
      long n_fallback_nearest = 0;
      for(size_t li=0; li<dl.line_positions.size(); li++)
      {
        result.value[li].resize( dl.line_positions[li].size() );
        for(size_t k=0; k<dl.line_positions[li].size(); k++)
        {
          const Vec3d& node_pos = dl.line_positions[li][k];
          double sum = 0.0; long n_in_cutoff = 0;
          double nearest_d = 1e300, nearest_v = 0.0;
          for(size_t a=0;a<n_atoms;a++)
          {
            const Vec3d atom_pos{ all_buf[4*a], all_buf[4*a+1], all_buf[4*a+2] };
            const double d = min_image_dist( node_pos, atom_pos );
            if( d < nearest_d ) { nearest_d = d; nearest_v = all_buf[4*a+3]; }
            if( d < cut ) { sum += all_buf[4*a+3]; ++n_in_cutoff; }
          }
          if( n_in_cutoff > 0 ) { result.value[li][k] = sum / static_cast<double>(n_in_cutoff); }
          else { result.value[li][k] = nearest_v; ++n_fallback_nearest; } // no atom within cutoff -- fall back to the single nearest one rather than leave the node without a real value
        }
      }

      long n_nodes = 0; for( const auto& v : result.value ) { n_nodes += static_cast<long>(v.size()); }
      lout << "dxa_project_atom_field: projected '" << *atom_field << "' onto " << n_nodes << " dislocation-line nodes across "
           << dl.line_positions.size() << " lines (" << n_atoms << " atoms gathered, " << n_fallback_nearest
           << " nodes fell back to nearest-atom, no atom within " << cut << " Ang)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Projects an arbitrary per-atom scalar field onto every dislocation-line node, for ParaView
visualization -- see this file's own header comment for the full mechanism. Run after
compute_dxa_mpi_stitch_lines (needs the FINAL, post-stitch line geometry); output
(DXALineNodeFieldValues) lives on rank 0 only, same convention as DXADislocationLines itself.

Usage example:

compute_dxa_mpi_stitch_lines: {}
dxa_project_atom_field: { atom_field: dxa_cluster_id, cutoff: 5.0 ang }
write_dxa_dislocation_lines: { filename: "paraview/dxa_lines" }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(dxa_project_atom_field)
  {
    OperatorNodeFactory::instance()->register_factory( "dxa_project_atom_field", make_grid_variant_operator< DXAProjectAtomField > );
  }

}
