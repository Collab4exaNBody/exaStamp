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

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/domain.h>

#include <exaStamp/delaunay/dxa_dislocation_lines.h>
#include <exaStamp/delaunay/dxa_lattice_correspondence.h>

#include <mpi.h>
#include <fstream>
#include <string>
#include <vector>

// Writes this operator's DXADislocationLines to OVITO's own "Crystal Analysis" (.ca) file format,
// so the result can be loaded directly into OVITO for visual comparison against its own DXA output
// -- reverse-engineered from OVITO's real (non-public) exporter/importer source, obtained directly
// by the user (CAExporter.cpp/CAImporter.cpp), not from public documentation.
//
// Deliberately minimal/simplified relative to a real OVITO-written file:
//   - Only one STRUCTURE_TYPE is declared ("bcc"), but its own BURGERS_VECTOR_FAMILY table is now a
//     faithful copy of OVITO's real one (cross-checked directly against a real OVITO-exported .ca
//     reference file's own BCC section) -- **this matters for more than cosmetics**: OVITO's CA
//     importer classifies each dislocation's "type" purely by matching its own Burgers vector against
//     the declared families (via the structure's symmetry group), there's no separate per-dislocation
//     type field in the DISLOCATIONS section itself. A first version of this operator declared only
//     the catch-all "Other" family (reference vector (0,0,0)) -- meaning EVERY dislocation, regardless
//     of its real Burgers vector, necessarily matched only that one family and showed up as "Other" in
//     OVITO, a real bug the user caught by inspecting OVITO's own classification panel. Fixed by
//     declaring the same 4 families OVITO's own BCC section has: Other, 1/2<111> (Perfect), <100>
//     (Hirth-like), <110>.
//   - Per-point trailing value (OVITO's own "core size"/local dislocation-core width) is now the real
//     per-point `DXADislocationLines::core_size` value (same data `smooth_dxa_dislocation_lines`
//     already consumes for its own coarsening weighting) -- a first version hardcoded this to a
//     constant 0 for every point, which isn't wrong for import (OVITO tolerates it) but throws away
//     real, already-computed data this writer has no good reason not to use.
//   - DISLOCATION_JUNCTIONS: this operator doesn't reconstruct real multi-way junction connectivity
//     (see compute_dxa_circuit_sweep.cpp's own "ponytail" note), so every dislocation is written as
//     a trivial self-referential 2-cycle (both of its own ends point back to each other) rather than
//     real cross-dislocation junction rings. OVITO will render every line as an independent,
//     unconnected segment -- exactly what this operator actually knows, no more.
//   - Works with any number of MPI ranks: the .ca format has no native multi-piece convention (unlike
//     write_dxa_dislocation_lines.cpp's per-rank-piece + .pvtu convention), so every rank's own
//     dislocation lines are serialized into a flat double buffer (self-delimiting: burgers.x/y/z,
//     npoints, then npoints*4 doubles) and gathered to rank 0 via MPI_Gatherv; only rank 0 writes
//     the file, walking each rank's own segment of the concatenated buffer in turn and assigning a
//     fresh, globally-sequential index to each line as it's written (the original per-rank index
//     doesn't need to survive the gather -- nothing else in this file cross-references it).
//   - No DEFECT_MESH section (OVITO's own interface-mesh triangulation) -- write_interface_mesh.cpp
//     already covers that separately, in VTK format.
//   - CLUSTER_SIZE: **real bug, found and fixed** -- the user noticed every dislocation showed up
//     as "Other" in OVITO regardless of its real type, and suspected BURGERS_VECTOR_FAMILY. That
//     table was indeed incomplete (fixed above, now a faithful copy of OVITO's real 4-family BCC
//     table) but turned out NOT to be the actual cause: bisected directly against a real OVITO-
//     exported .ca reference by splicing/patching individual fields (confirmed via `ovitos`,
//     `data.tables['disloc-lengths']`) that neither the Burgers vector's exact sign/permutation, nor
//     CLUSTER_ORIENTATION, nor the STRUCTURE_TYPE's own numeric ID mattered -- but forcing a known-
//     good reference file's `CLUSTER_SIZE` from its real value down to 0 alone broke classification
//     completely (every dislocation became "Other"). This operator always wrote a hardcoded
//     `CLUSTER_SIZE 0` (never had the real atom count available at all) -- OVITO's importer
//     apparently treats a zero-size cluster as invalid for Burgers-vector-family matching. Fixed by
//     taking `DXALatticeClusters` as a new input and writing cluster 1's own real `atom_count`
//     (summed across MPI ranks via `MPI_Reduce`, since `compute_dxa_lattice_clusters` is per-rank
//     local) instead. `CLUSTER_ORIENTATION` is also now cluster 1's own real least-squares
//     orientation fit (rank 0's own local value; this field is display/diagnostic only in OVITO, not
//     used by the elastic mapping itself, so no cross-rank reconciliation is needed) instead of a
//     hardcoded identity matrix.
namespace exaStamp
{
  using namespace exanb;

  class WriteDXACAFile : public OperatorNode
  {
    ADD_SLOT( MPI_Comm             , mpi                    , INPUT , REQUIRED );
    ADD_SLOT( Domain               , domain                 , INPUT , REQUIRED );
    ADD_SLOT( DXADislocationLines  , dxa_dislocation_lines  , INPUT , REQUIRED );
    ADD_SLOT( DXALatticeClusters   , dxa_lattice_clusters   , INPUT , REQUIRED , DocString{"Only cluster 1's own atom_count/orientation are used (CLUSTER_SIZE/CLUSTER_ORIENTATION) -- see this file's own header comment for why CLUSTER_SIZE matters for OVITO's own Burgers-vector-family classification"} );
    ADD_SLOT( std::string          , filename               , INPUT , std::string("dislocation_lines.ca") );

  public:
    inline void execute () override final
    {
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      const DXADislocationLines& dl = *dxa_dislocation_lines;
      const size_t n_local = dl.line_positions.empty() ? dl.lines.size() : dl.line_positions.size();

      // Serialize this rank's own lines into one flat, self-delimiting double buffer: per line,
      // burgers.x/y/z then npoints then npoints*4 doubles (x,y,z,core_size) -- "arcs" isn't included
      // since this writer always emits exactly 1 (see below), so there's nothing to transmit for it.
      std::vector<double> local_buf;
      for(size_t i=0;i<n_local;i++)
      {
        const Vec3d& b = dl.burgers_vector[i];
        local_buf.push_back(b.x); local_buf.push_back(b.y); local_buf.push_back(b.z);
        if( !dl.line_positions.empty() )
        {
          const auto& pts = dl.line_positions[i];
          const bool has_core_size = !dl.core_size.empty() && dl.core_size[i].size() == pts.size();
          local_buf.push_back( static_cast<double>( pts.size() ) );
          for(size_t k=0;k<pts.size();k++)
          {
            const Vec3d& p = pts[k];
            local_buf.push_back(p.x); local_buf.push_back(p.y); local_buf.push_back(p.z);
            local_buf.push_back( has_core_size ? static_cast<double>( dl.core_size[i][k] ) : 0. );
          }
        }
        else
        {
          local_buf.push_back(0.);
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

      // CLUSTER_SIZE must be the REAL total atom count, not 0 -- see this file's own header comment
      // (OVITO's importer treats a zero-size cluster as invalid for Burgers-vector-family matching).
      // compute_dxa_lattice_clusters is per-rank local, so sum every rank's own cluster 1 atom_count.
      const auto& graph = dxa_lattice_clusters->graph;
      const int local_cluster_size = ( graph.clusters.size() > 1 ) ? graph.clusters[1].atom_count : 0;
      int total_cluster_size = 0;
      MPI_Reduce( &local_cluster_size, &total_cluster_size, 1, MPI_INT, MPI_SUM, 0, *mpi );

      if( rank != 0 ) { return; }

      // Display/diagnostic only in OVITO (not used by the elastic mapping itself), so rank 0's own
      // local orientation fit is a perfectly valid representative -- no cross-rank reconciliation
      // needed, unlike CLUSTER_SIZE.
      const Mat3d cluster_orientation = ( graph.clusters.size() > 1 ) ? graph.clusters[1].orientation : onika::math::make_identity_matrix();

      // Walk the concatenated buffer once, per rank segment (recvcounts/displs boundaries), to
      // decode every line's own (burgers, points) record -- self-delimiting via each record's own
      // npoints field, so rank boundaries don't need special handling beyond bounding the walk.
      struct DecodedLine { Vec3d burgers; std::vector<Vec3d> pts; std::vector<int32_t> core_size; };
      std::vector<DecodedLine> all_lines;
      for(int r=0;r<np;r++)
      {
        size_t pos = static_cast<size_t>( displs[r] );
        const size_t end = pos + static_cast<size_t>( recvcounts[r] );
        while( pos < end )
        {
          DecodedLine line;
          line.burgers = Vec3d{ all_buf[pos], all_buf[pos+1], all_buf[pos+2] };
          pos += 3;
          const size_t npoints = static_cast<size_t>( all_buf[pos] );
          pos += 1;
          line.pts.reserve(npoints);
          line.core_size.reserve(npoints);
          for(size_t k=0;k<npoints;k++)
          {
            line.pts.push_back( Vec3d{ all_buf[pos], all_buf[pos+1], all_buf[pos+2] } );
            line.core_size.push_back( static_cast<int32_t>( all_buf[pos+3] ) );
            pos += 4;
          }
          all_lines.push_back( std::move(line) );
        }
      }
      const size_t n = all_lines.size();

      std::ofstream f( *filename );
      f << "CA_FILE_VERSION 6\n";
      f << "CA_LIB_VERSION 0.0.0\n";
      // OVITO itself always declares all 5 of its own built-in structure types in a .ca file, even
      // when a given DXA run only ever analyzes one of them (confirmed by the user against a real
      // OVITO-exported reference file) -- reproduced verbatim here (same ids/names/reference vectors/
      // colors) for full format fidelity, not just the single structure this pipeline actually uses.
      // This also gives `CLUSTER_STRUCTURE` below its real, OVITO-matching id (bcc=3) instead of an
      // arbitrary 1 -- confirmed harmless either way for classification itself (see this file's own
      // header comment: the real fix for the "Other" bug was CLUSTER_SIZE, not this id), but matching
      // OVITO's own numbering is more faithful and costs nothing.
      f << "STRUCTURE_TYPES 5\n";
      f << "STRUCTURE_TYPE 1\n";
      f << "NAME fcc\n";
      f << "FULL_NAME FCC\n";
      f << "COLOR 0.4 1.0 0.4\n";
      f << "TYPE LATTICE\n";
      f << "BURGERS_VECTOR_FAMILIES 6\n";
      f << "BURGERS_VECTOR_FAMILY ID 0\nOther\n0.0 0.0 0.0\n0.9 0.2 0.2\n";
      f << "BURGERS_VECTOR_FAMILY ID 1\n1/2<110> (Perfect)\n0.5 0.5 0.0\n0.2 0.2 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 2\n1/6<112> (Shockley)\n0.1666666667 0.1666666667 0.3333333433\n0.0 1.0 0.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 3\n1/6<110> (Stair-rod)\n0.1666666667 0.1666666667 0.0\n1.0 0.0 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 4\n1/3<100> (Hirth)\n0.3333333333 0.0 0.0\n1.0 1.0 0.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 5\n1/3<111> (Frank)\n0.3333333333 0.3333333333 0.3333333333\n0.0 1.0 1.0\n";
      f << "END_STRUCTURE_TYPE\n";
      f << "STRUCTURE_TYPE 2\n";
      f << "NAME hcp\n";
      f << "FULL_NAME HCP\n";
      f << "COLOR 1.0 0.4 0.4\n";
      f << "TYPE LATTICE\n";
      f << "BURGERS_VECTOR_FAMILIES 6\n";
      f << "BURGERS_VECTOR_FAMILY ID 0\nOther\n0.0 0.0 0.0\n0.9 0.2 0.2\n";
      f << "BURGERS_VECTOR_FAMILY ID 1\n1/3<1-210>\n0.7071067812 0.0 0.0\n0.0 1.0 0.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 2\n<0001>\n0.0 0.0 1.1547005384\n0.2 0.2 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 3\n<1-100>\n0.0 1.2247448714 0.0\n1.0 0.0 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 4\n1/3<1-100>\n0.0 0.4082482905 0.0\n1.0 0.5 0.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 5\n1/3<1-213>\n0.7071067812 0.0 1.1547005384\n1.0 1.0 0.0\n";
      f << "END_STRUCTURE_TYPE\n";
      f << "STRUCTURE_TYPE 3\n";
      f << "NAME bcc\n";
      f << "FULL_NAME BCC\n";
      f << "COLOR 0.4 0.4 1.0\n";
      f << "TYPE LATTICE\n";
      f << "BURGERS_VECTOR_FAMILIES 4\n";
      f << "BURGERS_VECTOR_FAMILY ID 0\nOther\n0.0 0.0 0.0\n0.9 0.2 0.2\n";
      f << "BURGERS_VECTOR_FAMILY ID 1\n1/2<111>\n0.5 0.5 0.5\n0.0 1.0 0.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 2\n<100>\n1.0 0.0 0.0\n1.0 0.3 0.8\n";
      f << "BURGERS_VECTOR_FAMILY ID 3\n<110>\n1.0 1.0 0.0\n0.2 0.5 1.0\n";
      f << "END_STRUCTURE_TYPE\n";
      f << "STRUCTURE_TYPE 4\n";
      f << "NAME diamond\n";
      f << "FULL_NAME Cubic diamond\n";
      f << "COLOR 0.0745098069 0.6274510026 0.9960784316\n";
      f << "TYPE LATTICE\n";
      f << "BURGERS_VECTOR_FAMILIES 5\n";
      f << "BURGERS_VECTOR_FAMILY ID 0\nOther\n0.0 0.0 0.0\n0.9 0.2 0.2\n";
      f << "BURGERS_VECTOR_FAMILY ID 1\n1/2<110>\n0.5 0.5 0.0\n0.2 0.2 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 2\n1/6<112>\n0.1666666667 0.1666666667 0.3333333333\n0.0 1.0 0.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 3\n1/6<110>\n0.1666666667 0.1666666667 0.0\n1.0 0.0 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 4\n1/3<111>\n0.3333333333 0.3333333333 0.3333333333\n0.0 1.0 1.0\n";
      f << "END_STRUCTURE_TYPE\n";
      f << "STRUCTURE_TYPE 5\n";
      f << "NAME hex_diamond\n";
      f << "FULL_NAME Hexagonal diamond\n";
      f << "COLOR 0.9960784316 0.5372549295 0.0\n";
      f << "TYPE LATTICE\n";
      f << "BURGERS_VECTOR_FAMILIES 5\n";
      f << "BURGERS_VECTOR_FAMILY ID 0\nOther\n0.0 0.0 0.0\n0.9 0.2 0.2\n";
      f << "BURGERS_VECTOR_FAMILY ID 1\n1/3<1-210>\n0.7071067812 0.0 0.0\n0.0 1.0 0.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 2\n<0001>\n0.0 0.0 1.1547005384\n0.2 0.2 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 3\n<1-100>\n0.0 1.2247448714 0.0\n1.0 0.0 1.0\n";
      f << "BURGERS_VECTOR_FAMILY ID 4\n1/3<1-100>\n0.0 0.4082482905 0.0\n1.0 0.5 0.0\n";
      f << "END_STRUCTURE_TYPE\n";

      const Vec3d origin = domain->origin();
      const Vec3d ex = domain->xform() * Vec3d{ domain->bounds_size().x, 0., 0. };
      const Vec3d ey = domain->xform() * Vec3d{ 0., domain->bounds_size().y, 0. };
      const Vec3d ez = domain->xform() * Vec3d{ 0., 0., domain->bounds_size().z };
      f << "SIMULATION_CELL_ORIGIN " << origin.x << " " << origin.y << " " << origin.z << "\n";
      f << "SIMULATION_CELL_MATRIX\n";
      f << ex.x << " " << ey.x << " " << ez.x << "\n";
      f << ex.y << " " << ey.y << " " << ez.y << "\n";
      f << ex.z << " " << ey.z << " " << ez.z << "\n";
      f << "PBC_FLAGS " << (int)domain->periodic_boundary_x() << " " << (int)domain->periodic_boundary_y() << " " << (int)domain->periodic_boundary_z() << "\n";

      f << "CLUSTERS 1\n";
      f << "CLUSTER 1\n";
      f << "CLUSTER_STRUCTURE 3\n"; // bcc, matching OVITO's own numbering (see STRUCTURE_TYPE 3 above)
      f << "CLUSTER_ORIENTATION\n";
      f << cluster_orientation.m11 << " " << cluster_orientation.m12 << " " << cluster_orientation.m13 << "\n";
      f << cluster_orientation.m21 << " " << cluster_orientation.m22 << " " << cluster_orientation.m23 << "\n";
      f << cluster_orientation.m31 << " " << cluster_orientation.m32 << " " << cluster_orientation.m33 << "\n";
      f << "CLUSTER_COLOR 1.0 1.0 1.0\n";
      f << "CLUSTER_SIZE " << total_cluster_size << "\n";
      f << "END_CLUSTER\n";
      f << "CLUSTER_TRANSITIONS 0\n";

      f << "DISLOCATIONS " << n << "\n";
      for(size_t i=0;i<n;i++)
      {
        const DecodedLine& line = all_lines[i];
        f << i << "\n";
        f << line.burgers.x << " " << line.burgers.y << " " << line.burgers.z << "\n";
        f << "1\n";
        f << line.pts.size() << "\n";
        for(size_t k=0;k<line.pts.size();k++)
        {
          const Vec3d& p = line.pts[k];
          f << p.x << " " << p.y << " " << p.z << " " << line.core_size[k] << "\n";
        }
      }

      // no real junction connectivity reconstructed -- each dislocation's own two ends form a
      // trivial self-referential 2-cycle (see this file's own header comment).
      f << "DISLOCATION_JUNCTIONS\n";
      for(size_t i=0;i<n;i++)
      {
        f << "0 " << i << "\n";
        f << "1 " << i << "\n";
      }

      f << "DEFECT_MESH_VERTICES 0\n";
      f << "DEFECT_MESH_FACETS 0\n";
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes a DXADislocationLines to OVITO's own "Crystal Analysis" (.ca) file format, so the result can
be loaded directly into OVITO (File > Load Trajectory, "CA File" importer) for visual comparison
against OVITO's own DXA output. Works with any number of MPI ranks (gathers to rank 0, which alone
writes the file). See this file's own header comment for the specific simplifications (minimal
structure-type stub, trivial self-referential junction records).

Usage example:

compute_dxa_circuit_sweep: {}
write_dxa_ca_file: { filename: "paraview/dislocation_lines.ca" }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(write_dxa_ca_file)
  {
    OperatorNodeFactory::instance()->register_factory( "write_dxa_ca_file", make_simple_operator< WriteDXACAFile > );
  }

}
