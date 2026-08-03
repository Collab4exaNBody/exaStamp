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
#include <exaStamp/delaunay/dxa_dislocation_lines.h>

#include <mpi.h>
#include <fstream>
#include <sstream>
#include <filesystem>

// Writes this rank's DXADislocationLines (see compute_dxa_dislocation_lines) as a VTK
// unstructured grid piece -- one VTK_POLY_LINE (type 4) cell per extracted line (open segment or
// closed loop, the latter with its first vertex repeated at the end) plus one VTK_VERTEX (type 1)
// cell per junction atom -- plus (rank 0 only) a .pvtu master file. Same per-rank-piece +
// master-file convention as write_delaunay_vtk/write_interface_mesh. Only vertices actually
// referenced by a line or junction are emitted.
namespace exaStamp
{
  using namespace exanb;

  class WriteDXADislocationLines : public OperatorNode
  {
    ADD_SLOT( MPI_Comm             , mpi                    , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation  , INPUT , REQUIRED );
    ADD_SLOT( DXADislocationLines  , dxa_dislocation_lines  , INPUT , REQUIRED );
    ADD_SLOT( std::string          , filename               , INPUT , std::string("dislocation_lines") );

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
      const DXADislocationLines& dl = *dxa_dislocation_lines;
      const size_t n_lines = dl.lines.size();
      const size_t n_junctions = dl.junction_vertices.size();
      const size_t nb_cells = n_lines + n_junctions;

      // compact to only the vertices actually referenced by a line or junction
      std::vector<int64_t> vertex_remap( mesh.vertices.size(), -1 );
      std::vector<Vec3d> out_points;
      auto remap = [&]( uint32_t v ) -> uint32_t
      {
        if( vertex_remap[v] < 0 ) { vertex_remap[v] = static_cast<int64_t>( out_points.size() ); out_points.push_back( mesh.vertices[v] ); }
        return static_cast<uint32_t>( vertex_remap[v] );
      };

      std::vector<std::vector<uint32_t>> connectivity( nb_cells );
      std::vector<int> cell_kind( nb_cells ); // 0=open line, 1=closed loop, 2=junction
      std::vector<Vec3d> cell_burgers( nb_cells, Vec3d{0.,0.,0.} );
      // which physical dislocation this segment belongs to after merging (compute_dxa_circuit_sweep
      // only); -1 for junction cells and for other extractors that don't populate dislocation_id.
      const bool has_dislocation_id = dl.dislocation_id.size() == n_lines;
      std::vector<int32_t> cell_dislocation_id( nb_cells, -1 );

      // compute_dxa_circuit_sweep's own line vertices are synthetic swept-circuit centroids with
      // no corresponding tessellation vertex -- it populates line_positions directly instead of
      // lines (see DXADislocationLines' own field comments). Prefer that when present.
      const bool use_positions = !dl.line_positions.empty();
      for(size_t li=0; li<n_lines; li++)
      {
        if( use_positions )
        {
          connectivity[li].reserve( dl.line_positions[li].size() );
          for(const auto& p : dl.line_positions[li])
          {
            connectivity[li].push_back( static_cast<uint32_t>( out_points.size() ) );
            out_points.push_back( p );
          }
        }
        else
        {
          connectivity[li].reserve( dl.lines[li].size() );
          for(uint32_t v : dl.lines[li]) { connectivity[li].push_back( remap(v) ); }
        }
        cell_kind[li] = dl.is_loop[li] ? 1 : 0;
        cell_burgers[li] = dl.burgers_vector[li];
        if( has_dislocation_id ) { cell_dislocation_id[li] = dl.dislocation_id[li]; }
      }
      for(size_t ji=0; ji<n_junctions; ji++)
      {
        const size_t c = n_lines + ji;
        connectivity[c] = { remap( dl.junction_vertices[ji] ) };
        cell_kind[c] = 2;
      }

      const size_t nb_points = out_points.size();

      // per-rank piece file
      {
        std::ostringstream oss; oss << outputDir << "/piece" << rank << ".vtu";
        std::ofstream vtu( oss.str() );
        vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        vtu << "  <UnstructuredGrid>\n";
        vtu << "    <Piece NumberOfPoints=\""<<nb_points<<"\" NumberOfCells=\""<<nb_cells<<"\">\n";
        vtu << "      <Points>\n";
        vtu << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for(const auto& p : out_points) { vtu << "          " << p.x << " " << p.y << " " << p.z << "\n"; }
        vtu << "        </DataArray>\n";
        vtu << "      </Points>\n";
        vtu << "      <Cells>\n";
        vtu << "        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
        for(const auto& c : connectivity) { vtu << "          "; for(uint32_t idx : c) { vtu << idx << " "; } vtu << "\n"; }
        vtu << "        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n          ";
        { uint64_t off=0; for(const auto& c : connectivity) { off += c.size(); vtu << off << " "; } }
        vtu << "\n        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
        for(int k : cell_kind) { vtu << ( (k==2) ? 1 : 4 ) << " "; } // 1=VTK_VERTEX, 4=VTK_POLY_LINE
        vtu << "\n        </DataArray>\n";
        vtu << "      </Cells>\n";
        vtu << "      <CellData Scalars=\"kind\">\n";
        vtu << "        <DataArray type=\"Int32\" Name=\"kind\" format=\"ascii\">\n          "; // 0=open line, 1=closed loop, 2=junction
        for(int k : cell_kind) { vtu << k << " "; }
        vtu << "\n        </DataArray>\n";
        vtu << "        <DataArray type=\"Float64\" Name=\"burgers_vector\" NumberOfComponents=\"3\" format=\"ascii\">\n          ";
        for(const auto& b : cell_burgers) { vtu << b.x << " " << b.y << " " << b.z << "  "; }
        vtu << "\n        </DataArray>\n";
        if( has_dislocation_id )
        {
          vtu << "        <DataArray type=\"Int32\" Name=\"dislocation_id\" format=\"ascii\">\n          ";
          for(int32_t d : cell_dislocation_id) { vtu << d << " "; }
          vtu << "\n        </DataArray>\n";
        }
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
        pvtu << "    <PCellData Scalars=\"kind\">\n";
        pvtu << "      <PDataArray type=\"Int32\" Name=\"kind\"/>\n";
        pvtu << "      <PDataArray type=\"Float64\" Name=\"burgers_vector\" NumberOfComponents=\"3\"/>\n";
        if( has_dislocation_id ) { pvtu << "      <PDataArray type=\"Int32\" Name=\"dislocation_id\"/>\n"; }
        pvtu << "    </PCellData>\n";
        for(int i=0;i<np;i++) { pvtu << "    <Piece Source=\""<<basename<<"/piece"<<i<<".vtu\"/>\n"; }
        pvtu << "  </PUnstructuredGrid>\n";
        pvtu << "</VTKFile>\n";
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes a DXADislocationLines (see compute_dxa_dislocation_lines) to VTK unstructured grid format
for ParaView -- one VTK_POLY_LINE cell per extracted line (open segment or closed loop) plus one
VTK_VERTEX cell per junction atom, with per-cell "kind" (0=open line, 1=closed loop, 2=junction) and
"burgers_vector" CellData. When the extractor populates DXADislocationLines::dislocation_id (only
compute_dxa_circuit_sweep does, after merging clean two-way junction chains back into one physical
dislocation), also emits an Int32 "dislocation_id" CellData array so ParaView can color/group merged
segments together. Same per-rank-piece + master-.pvtu convention as write_interface_mesh.

Usage example:

compute_dxa_dislocation_lines: { min_core_depth: 1 }
write_dxa_dislocation_lines: { filename: "paraview/dislocation_lines" }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(write_dxa_dislocation_lines)
  {
    OperatorNodeFactory::instance()->register_factory( "write_dxa_dislocation_lines", make_simple_operator< WriteDXADislocationLines > );
  }

}
