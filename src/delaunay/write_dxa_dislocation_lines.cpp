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
//
// dxa_line_node_field_values (optional, auto-wired from dxa_project_atom_field) additionally emits
// a per-point PointData array, named after DXALineNodeFieldValues::field_name -- only lines (not
// junction points) have a real projected value; junction points get 0. Only rank 0 ever has real
// content (dxa_project_atom_field's own convention), so this array is empty/omitted on every other
// rank's own piece, same as DXADislocationLines itself post-stitch.
namespace exaStamp
{
  using namespace exanb;

  class WriteDXADislocationLines : public OperatorNode
  {
    ADD_SLOT( MPI_Comm               , mpi                        , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation   , delaunay_tessellation      , INPUT , REQUIRED );
    ADD_SLOT( DXADislocationLines    , dxa_dislocation_lines      , INPUT , REQUIRED );
    ADD_SLOT( DXALineNodeFieldValues , dxa_line_node_field_values , INPUT , OPTIONAL , DocString{"Per-line-node projected atom field (dxa_project_atom_field) -- when present, emitted as an extra PointData array named after its own field_name"} );
    ADD_SLOT( std::string            , filename                   , INPUT , std::string("dislocation_lines") );

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
      const size_t n_xrank_junctions = dl.junction_positions.size();
      const size_t nb_cells = n_lines + n_junctions + n_xrank_junctions;

      // compact to only the vertices actually referenced by a line or junction
      std::vector<int64_t> vertex_remap( mesh.vertices.size(), -1 );
      std::vector<Vec3d> out_points;
      auto remap = [&]( uint32_t v ) -> uint32_t
      {
        if( vertex_remap[v] < 0 ) { vertex_remap[v] = static_cast<int64_t>( out_points.size() ); out_points.push_back( mesh.vertices[v] ); }
        return static_cast<uint32_t>( vertex_remap[v] );
      };

      std::vector<std::vector<uint32_t>> connectivity( nb_cells );
      std::vector<int> cell_kind( nb_cells ); // 0=open line, 1=closed loop, 2=junction, 3=cross-rank junction
      std::vector<Vec3d> cell_burgers( nb_cells, Vec3d{0.,0.,0.} );
      // which physical dislocation this segment belongs to after merging (compute_dxa_circuit_sweep
      // only); -1 for junction cells and for other extractors that don't populate dislocation_id.
      const bool has_dislocation_id = dl.dislocation_id.size() == n_lines;
      std::vector<int32_t> cell_dislocation_id( nb_cells, -1 );

      const bool has_node_field = dxa_line_node_field_values.has_value() && dxa_line_node_field_values->value.size() == n_lines;
      std::vector<double> point_field_value; // parallel to out_points, only meaningful if has_node_field

      // compute_dxa_circuit_sweep's own line vertices are synthetic swept-circuit centroids with
      // no corresponding tessellation vertex -- it populates line_positions directly instead of
      // lines (see DXADislocationLines' own field comments). Prefer that when present.
      const bool use_positions = !dl.line_positions.empty();
      for(size_t li=0; li<n_lines; li++)
      {
        if( use_positions )
        {
          connectivity[li].reserve( dl.line_positions[li].size() );
          const bool line_has_field = has_node_field && dxa_line_node_field_values->value[li].size() == dl.line_positions[li].size();
          for(size_t k=0; k<dl.line_positions[li].size(); k++)
          {
            connectivity[li].push_back( static_cast<uint32_t>( out_points.size() ) );
            out_points.push_back( dl.line_positions[li][k] );
            if( has_node_field ) { point_field_value.push_back( line_has_field ? dxa_line_node_field_values->value[li][k] : 0.0 ); }
          }
        }
        else
        {
          connectivity[li].reserve( dl.lines[li].size() );
          for(uint32_t v : dl.lines[li])
          {
            const uint32_t out_idx = remap(v);
            connectivity[li].push_back( out_idx );
            if( has_node_field ) { point_field_value.resize( out_points.size(), 0.0 ); }
          }
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
      // Cross-rank N-way junctions reconstructed by compute_dxa_mpi_stitch_lines -- a real 3D
      // position, not a DelaunayTessellation vertex (see DXADislocationLines::junction_positions'
      // own doc comment for why), so appended directly to out_points rather than through remap().
      for(size_t ji=0; ji<n_xrank_junctions; ji++)
      {
        const size_t c = n_lines + n_junctions + ji;
        connectivity[c] = { static_cast<uint32_t>( out_points.size() ) };
        out_points.push_back( dl.junction_positions[ji] );
        cell_kind[c] = 3;
      }
      if( has_node_field ) { point_field_value.resize( out_points.size(), 0.0 ); } // cover any junction points appended above

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
        if( has_node_field )
        {
          vtu << "      <PointData Scalars=\""<<dxa_line_node_field_values->field_name<<"\">\n";
          vtu << "        <DataArray type=\"Float64\" Name=\""<<dxa_line_node_field_values->field_name<<"\" format=\"ascii\">\n          ";
          for(double v : point_field_value) { vtu << v << " "; }
          vtu << "\n        </DataArray>\n";
          vtu << "      </PointData>\n";
        }
        vtu << "      <Cells>\n";
        vtu << "        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
        for(const auto& c : connectivity) { vtu << "          "; for(uint32_t idx : c) { vtu << idx << " "; } vtu << "\n"; }
        vtu << "        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n          ";
        { uint64_t off=0; for(const auto& c : connectivity) { off += c.size(); vtu << off << " "; } }
        vtu << "\n        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
        for(int k : cell_kind) { vtu << ( (k==2||k==3) ? 1 : 4 ) << " "; } // 1=VTK_VERTEX, 4=VTK_POLY_LINE
        vtu << "\n        </DataArray>\n";
        vtu << "      </Cells>\n";
        vtu << "      <CellData Scalars=\"kind\">\n";
        vtu << "        <DataArray type=\"Int32\" Name=\"kind\" format=\"ascii\">\n          "; // 0=open line, 1=closed loop, 2=junction, 3=cross-rank junction
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
        if( has_node_field )
        {
          pvtu << "    <PPointData Scalars=\""<<dxa_line_node_field_values->field_name<<"\">\n";
          pvtu << "      <PDataArray type=\"Float64\" Name=\""<<dxa_line_node_field_values->field_name<<"\"/>\n";
          pvtu << "    </PPointData>\n";
        }
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
