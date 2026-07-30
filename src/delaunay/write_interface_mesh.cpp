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
#include <sstream>
#include <filesystem>

// Writes this rank's InterfaceMesh (see compute_interface_mesh, DXA pipeline step (v)) as a VTK
// unstructured grid piece (triangles, not tetrahedra), plus (rank 0 only) a .pvtu master file --
// same per-rank-piece + master-file convention as write_delaunay_vtk.cpp. Only vertices actually
// referenced by an interface triangle are emitted (interface triangles are typically a small
// fraction of DelaunayTessellation::vertices, so writing the full vertex list here would be
// mostly-unreferenced points).
namespace exaStamp
{
  using namespace exanb;

  class WriteInterfaceMesh : public OperatorNode
  {
    ADD_SLOT( MPI_Comm             , mpi                  , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation, INPUT , REQUIRED );
    ADD_SLOT( InterfaceMesh        , interface_mesh       , INPUT , REQUIRED );
    ADD_SLOT( std::string          , filename             , INPUT , std::string("interface") );

  public:
    inline void execute () override final
    {
      namespace fs = std::filesystem;

      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      const std::string filepath = *filename;
      std::string dirname, basename = filepath;
      const size_t lastSlash = filepath.rfind('/');
      if( lastSlash != std::string::npos )
      {
        dirname = filepath.substr(0, lastSlash);
        basename = filepath.substr(lastSlash + 1);
      }
      if( basename.rfind('.') != std::string::npos ) { basename = basename.substr(0, basename.rfind('.')); }
      const std::string outputDir = dirname.empty() ? basename : (dirname + "/" + basename);

      if( rank == 0 )
      {
        fs::remove_all(outputDir);
        std::error_code ec;
        fs::create_directories(outputDir, ec);
      }
      MPI_Barrier(*mpi);

      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const InterfaceMesh& iface = *interface_mesh;
      const size_t nb_cells = iface.triangles.size();

      // compact to only the vertices actually referenced by an interface triangle
      std::vector<int64_t> vertex_remap( mesh.vertices.size(), -1 );
      std::vector<Vec3d> out_vertices;
      std::vector<std::array<uint32_t,3>> out_triangles( nb_cells );
      for(size_t c=0;c<nb_cells;c++)
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
      const size_t nb_points = out_vertices.size();

      // per-rank piece file
      {
        std::ostringstream oss; oss << outputDir << "/piece" << rank << ".vtu";
        std::ofstream vtu( oss.str() );
        vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        vtu << "  <UnstructuredGrid>\n";
        vtu << "    <Piece NumberOfPoints=\""<<nb_points<<"\" NumberOfCells=\""<<nb_cells<<"\">\n";
        vtu << "      <Points>\n";
        vtu << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for(const auto& v : out_vertices) { vtu << "          " << v.x << " " << v.y << " " << v.z << "\n"; }
        vtu << "        </DataArray>\n";
        vtu << "      </Points>\n";
        vtu << "      <Cells>\n";
        vtu << "        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
        for(const auto& tri : out_triangles) { vtu << "          " << tri[0] << " " << tri[1] << " " << tri[2] << "\n"; }
        vtu << "        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n          ";
        for(size_t c=0;c<nb_cells;c++) { vtu << (3*(c+1)) << " "; }
        vtu << "\n        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
        for(size_t c=0;c<nb_cells;c++) { vtu << "5 "; } // VTK_TRIANGLE
        vtu << "\n        </DataArray>\n";
        vtu << "      </Cells>\n";
        vtu << "      <CellData Scalars=\"rank\">\n";
        vtu << "        <DataArray type=\"Int32\" Name=\"rank\" format=\"ascii\">\n          ";
        for(size_t c=0;c<nb_cells;c++) { vtu << rank << " "; }
        vtu << "\n        </DataArray>\n";
        vtu << "      </CellData>\n";
        vtu << "    </Piece>\n";
        vtu << "  </UnstructuredGrid>\n";
        vtu << "</VTKFile>\n";
      }

      // rank-0 master file, referencing every rank's piece
      if( rank == 0 )
      {
        std::ofstream pvtu( outputDir + ".pvtu" );
        pvtu << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        pvtu << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
        pvtu << "    <PPoints>\n";
        pvtu << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
        pvtu << "    </PPoints>\n";
        pvtu << "    <PCells>\n";
        pvtu << "      <PDataArray type=\"UInt32\" Name=\"connectivity\"/>\n";
        pvtu << "      <PDataArray type=\"UInt64\" Name=\"offsets\"/>\n";
        pvtu << "      <PDataArray type=\"UInt8\" Name=\"types\"/>\n";
        pvtu << "    </PCells>\n";
        pvtu << "    <PCellData Scalars=\"rank\">\n";
        pvtu << "      <PDataArray type=\"Int32\" Name=\"rank\"/>\n";
        pvtu << "    </PCellData>\n";
        for(int i=0;i<np;i++) { pvtu << "    <Piece Source=\""<<basename<<"/piece"<<i<<".vtu\"/>\n"; }
        pvtu << "  </PUnstructuredGrid>\n";
        pvtu << "</VTKFile>\n";
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes an InterfaceMesh (see compute_interface_mesh) to VTK unstructured grid format (triangles)
for ParaView. Same per-rank-piece + master-.pvtu convention as write_delaunay_vtk. Only vertices
actually referenced by an interface triangle are emitted, so the point count is unrelated to the
full DelaunayTessellation's own vertex count.

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
  ONIKA_AUTORUN_INIT(write_interface_mesh)
  {
    OperatorNodeFactory::instance()->register_factory( "write_interface_mesh", make_simple_operator< WriteInterfaceMesh > );
  }

}
