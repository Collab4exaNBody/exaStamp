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

#include <exaStamp/delaunay/delaunay_tessellation.h>
#include <exaStamp/delaunay/interface_mesh.h>

#include <mpi.h>
#include <fstream>
#include <vector>
#include <array>
#include <string>

// Writes this rank's InterfaceMesh in the SAME legacy VTK format OVITO Pro itself writes for its
// own interface mesh export (`export_file(..., "vtk/trimesh", ...)`, cross-checked byte-for-byte
// structure against a real OVITO-exported reference file: `# vtk DataFile Version 3.0` header,
// ASCII DATASET UNSTRUCTURED_GRID, POINTS <n> double, CELLS <n> <n*4> ("3 v0 v1 v2" per row,
// VTK_TRIANGLE), CELL_TYPES <n> (all "5"), then CELL_DATA/POINT_DATA each carrying one
// `SCALARS cap unsigned_char` field. Unlike write_interface_mesh.cpp's own per-rank-piece + .pvtu
// XML convention, legacy VTK has no multi-piece convention at all -- it's inherently a single,
// self-contained file -- so this operator actually GATHERS every rank's mesh to rank 0 (real
// MPI_Gatherv, not the per-piece trick), remapping each rank's own local vertex indices by that
// rank's running point-count offset before writing, and only rank 0 touches the filesystem.
//
// Honest caveat on the "cap" field: OVITO's own interface mesh synthesizes extra "cap" triangles to
// close off the surface at a non-periodic domain boundary, so a real OVITO file can have cap==1
// cells (confirmed on the screw-dipole reference: 506/3004 cells capped). This operator's own
// InterfaceMesh doesn't synthesize any such triangles -- an edge at this rank's own tessellation
// boundary is simply left open (see InterfaceMesh::edge_triangles's own doc comment) -- so `cap` is
// always written as 0 here. This writer matches OVITO's FILE FORMAT exactly (loads into OVITO/
// ParaView identically), not its capping semantics; synthesizing real cap triangles would be a
// separate, substantially bigger feature.
namespace exaStamp
{
  using namespace exanb;

  class WriteOvitoInterfaceMesh : public OperatorNode
  {
    ADD_SLOT( MPI_Comm             , mpi                  , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation, INPUT , REQUIRED );
    ADD_SLOT( InterfaceMesh        , interface_mesh       , INPUT , REQUIRED );
    ADD_SLOT( std::string          , filename             , INPUT , std::string("interface_mesh_ovito.vtk") );

  public:
    inline void execute () override final
    {
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const InterfaceMesh& iface = *interface_mesh;
      const size_t nb_cells_local = iface.triangles.size();

      // compact to only the vertices actually referenced by an interface triangle, same as
      // write_interface_mesh.cpp -- LOCAL indices for now, rebased to global once gathered.
      std::vector<int64_t> vertex_remap( mesh.vertices.size(), -1 );
      std::vector<Vec3d> out_vertices;
      std::vector<std::array<uint32_t,3>> out_triangles( nb_cells_local );
      for(size_t c=0;c<nb_cells_local;c++)
      {
        for(int lv=0;lv<3;lv++)
        {
          const uint32_t v = iface.triangles[c][lv];
          if( vertex_remap[v] < 0 )
          {
            vertex_remap[v] = static_cast<int64_t>( out_vertices.size() );
            out_vertices.push_back( mesh.vertices[v] );
          }
          out_triangles[c][lv] = static_cast<uint32_t>( vertex_remap[v] );
        }
      }
      const int nb_points_local = static_cast<int>( out_vertices.size() );
      const int nb_cells_local_i = static_cast<int>( nb_cells_local );

      std::vector<int> points_per_rank, cells_per_rank;
      if( rank == 0 ) { points_per_rank.resize(np); cells_per_rank.resize(np); }
      MPI_Gather( &nb_points_local, 1, MPI_INT, points_per_rank.data(), 1, MPI_INT, 0, *mpi );
      MPI_Gather( &nb_cells_local_i, 1, MPI_INT, cells_per_rank.data(), 1, MPI_INT, 0, *mpi );

      // flatten this rank's own points (3 doubles/point) and triangle connectivity (3 uint32/cell,
      // still LOCAL indices at this point) for the two gather calls below.
      std::vector<double> points_flat( 3 * nb_points_local );
      for(int i=0;i<nb_points_local;i++)
      {
        points_flat[3*i+0] = out_vertices[i].x;
        points_flat[3*i+1] = out_vertices[i].y;
        points_flat[3*i+2] = out_vertices[i].z;
      }
      std::vector<uint32_t> tris_flat( 3 * nb_cells_local );
      for(size_t i=0;i<nb_cells_local;i++)
      {
        tris_flat[3*i+0] = out_triangles[i][0];
        tris_flat[3*i+1] = out_triangles[i][1];
        tris_flat[3*i+2] = out_triangles[i][2];
      }

      std::vector<int> point_recvcounts, point_displs, cell_recvcounts, cell_displs;
      int total_points = 0, total_cells = 0;
      if( rank == 0 )
      {
        point_recvcounts.resize(np); point_displs.resize(np);
        cell_recvcounts.resize(np);  cell_displs.resize(np);
        int poff=0, coff=0;
        for(int r=0;r<np;r++)
        {
          point_recvcounts[r] = 3 * points_per_rank[r];
          point_displs[r] = poff;
          poff += point_recvcounts[r];
          cell_recvcounts[r] = 3 * cells_per_rank[r];
          cell_displs[r] = coff;
          coff += cell_recvcounts[r];
        }
        total_points = poff / 3;
        total_cells = coff / 3;
      }

      std::vector<double> all_points( rank==0 ? 3*static_cast<size_t>(total_points) : 0 );
      std::vector<uint32_t> all_tris_local( rank==0 ? 3*static_cast<size_t>(total_cells) : 0 );
      MPI_Gatherv( points_flat.data(), 3*nb_points_local, MPI_DOUBLE,
                   rank==0 ? all_points.data() : nullptr,
                   rank==0 ? point_recvcounts.data() : nullptr,
                   rank==0 ? point_displs.data() : nullptr,
                   MPI_DOUBLE, 0, *mpi );
      MPI_Gatherv( tris_flat.data(), 3*nb_cells_local_i, MPI_UINT32_T,
                   rank==0 ? all_tris_local.data() : nullptr,
                   rank==0 ? cell_recvcounts.data() : nullptr,
                   rank==0 ? cell_displs.data() : nullptr,
                   MPI_UINT32_T, 0, *mpi );

      if( rank != 0 ) { return; }

      // rebase each rank's own LOCAL vertex indices by that rank's running point-count offset --
      // cell_displs[r]/3 gives that rank's own triangle offset (in triangles, not ints), and the
      // matching point offset is point_displs[r]/3 (in points, not doubles).
      std::vector<uint32_t> all_tris( all_tris_local.size() );
      {
        int r = 0;
        for(size_t i=0;i<all_tris_local.size()/3;i++)
        {
          while( r+1 < np && static_cast<size_t>(cell_displs[r+1]/3) <= i ) { ++r; }
          const uint32_t point_offset = static_cast<uint32_t>( point_displs[r] / 3 );
          all_tris[3*i+0] = all_tris_local[3*i+0] + point_offset;
          all_tris[3*i+1] = all_tris_local[3*i+1] + point_offset;
          all_tris[3*i+2] = all_tris_local[3*i+2] + point_offset;
        }
      }

      std::ofstream f( *filename );
      f << "# vtk DataFile Version 3.0\n";
      f << "Triangle surface mesh written by exaStamp write_ovito_interface_mesh\n";
      f << "ASCII\n";
      f << "DATASET UNSTRUCTURED_GRID\n";
      f << "POINTS " << total_points << " double\n";
      for(int i=0;i<total_points;i++)
      {
        f << all_points[3*i+0] << " " << all_points[3*i+1] << " " << all_points[3*i+2] << "\n";
      }
      f << "CELLS " << total_cells << " " << (4*total_cells) << "\n";
      for(int i=0;i<total_cells;i++)
      {
        f << "3 " << all_tris[3*i+0] << " " << all_tris[3*i+1] << " " << all_tris[3*i+2] << "\n";
      }
      f << "CELL_TYPES " << total_cells << "\n";
      for(int i=0;i<total_cells;i++) { f << "5\n"; }
      f << "CELL_DATA " << total_cells << "\n";
      f << "SCALARS cap unsigned_char\n";
      f << "LOOKUP_TABLE default\n";
      for(int i=0;i<total_cells;i++) { f << "0\n"; }
      f << "POINT_DATA " << total_points << "\n";
      f << "SCALARS cap unsigned_char\n";
      f << "LOOKUP_TABLE default\n";
      for(int i=0;i<total_points;i++) { f << "0\n"; }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes an InterfaceMesh (see compute_interface_mesh / compute_dxa_elastic_interface_mesh) to the
same legacy VTK file format OVITO Pro itself writes for its own interface-mesh export -- loads
directly into OVITO or ParaView. Works with any number of MPI ranks: gathers every rank's mesh to
rank 0 (legacy VTK has no multi-piece convention, unlike write_interface_mesh.cpp's XML .pvtu), only
rank 0 writes the file. The per-cell/per-point "cap" scalar is always 0 here (this operator doesn't
synthesize OVITO-style boundary-closing cap triangles) -- see this file's own header comment.

Usage example:

compute_dxa_elastic_interface_mesh: {}
write_ovito_interface_mesh: { filename: "ovitodata/our_interface_mesh.vtk" }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(write_ovito_interface_mesh)
  {
    OperatorNodeFactory::instance()->register_factory( "write_ovito_interface_mesh", make_simple_operator< WriteOvitoInterfaceMesh > );
  }

}
