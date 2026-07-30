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

#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exaStamp/delaunay/delaunay_tessellation.h>
#include <exaStamp/delaunay/dxa_edge_vectors.h>

#include <mpi.h>
#include <fstream>
#include <sstream>
#include <filesystem>

// Writes this rank's DelaunayTessellation (see compute_delaunay) as a VTK unstructured grid
// piece, plus (rank 0 only) a .pvtu master file referencing every rank's piece -- same
// subdirectory/master-file convention as exaNBody's write_grid_vtk.cpp, so opening the single
// .pvtu file in ParaView shows the whole system's tessellation across every MPI rank at once.
// Each rank only ever wrote tetrahedra it owns (see compute_delaunay's centroid-in-owned-cell
// rule), so the combined result has neither gaps nor duplicated cells at rank boundaries.
//
// color_fields (optional list) generalizes this to any per-particle scalar field (PTM's ptm_type
// via ptm_fields, or any other field::mk_generic_real): each tetrahedron gets one CellData value
// per listed field, the average of its 4 vertices' field value, looked up via
// DelaunayTessellation::vertex_particle_index (see compute_delaunay.cpp) -- writing several at
// once lets ParaView switch the active coloring between them without re-running the pipeline.
// A vertex that happens to be a ghost particle (a tet near this rank's boundary can have some of
// its 4 vertices be ghost copies) only has a meaningful value if something already synchronized
// that field's ghost copies; otherwise it reads as whatever default that ghost cell's field has
// (typically 0), which biases boundary tets' average slightly -- acceptable for visualization,
// not corrected here.
//
// dxa_tet_classification (optional, auto-wired) additionally writes DXA step (iv)'s good/bad
// classification (compute_dxa_tet_classification) as a direct "dxa_good" CellData array (0/1) --
// already per-tetrahedron, unlike color_fields, so no vertex averaging involved.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class WriteDelaunayVTK : public OperatorNode
  {
    ADD_SLOT( MPI_Comm             , mpi                  , INPUT , REQUIRED );
    ADD_SLOT( GridT                , grid                 , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation, INPUT , REQUIRED );
    ADD_SLOT( std::string          , filename             , INPUT , std::string("delaunay") );
    ADD_SLOT( std::vector<std::string> , color_fields     , INPUT , std::vector<std::string>{} , DocString{"Names of per-particle scalar fields (field::mk_generic_real, e.g. PTM's ptm_type written by ptm_fields) to average over each tetrahedron's 4 vertices and write as extra CellData arrays (one per field), so ParaView can switch the active coloring between them. Leave empty to skip."} );
    ADD_SLOT( DXATetClassification , dxa_tet_classification , INPUT , OPTIONAL , DocString{"DXA step (iv) good/bad tetrahedron classification (compute_dxa_tet_classification) -- auto-wires in if that operator ran earlier in the pipeline. Written as a direct 0/1 CellData array (\"dxa_good\"), no per-vertex averaging needed since it's already per-tetrahedron."} );

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
      const size_t nb_points = mesh.vertices.size();
      const size_t nb_cells = mesh.tetrahedra.size();

      // Optional per-tet coloring: average each listed field over a tet's 4 vertices, looked up
      // via vertex_particle_index. Needs the deep cells_accessor() (not the raw cells() pointer)
      // since mk_generic_real is a runtime/dynamic field, not a compile-time FieldId.
      std::vector<std::vector<double>> tet_color( color_fields->size() );
      if( ! color_fields->empty() )
      {
        auto cells = grid->cells_accessor();
        const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
        const size_t n_cells_grid = grid->number_of_cells();
        const size_t n_particles = grid->number_of_particles();

        for(size_t f=0; f<color_fields->size(); f++)
        {
          auto field_acc = grid->field_const_accessor( field::mk_generic_real( (*color_fields)[f] ) );

          std::vector<double> flat_values( n_particles, 0.0 );
          for(size_t c=0;c<n_cells_grid;c++)
          {
            const size_t np = cells[c].size();
            for(size_t p=0;p<np;p++) { flat_values[ cell_particle_offset[c] + p ] = cells[c][field_acc][p]; }
          }

          tet_color[f].resize( nb_cells );
          for(size_t c=0;c<nb_cells;c++)
          {
            double sum = 0.0;
            for(int lv=0;lv<4;lv++) { sum += flat_values[ mesh.vertex_particle_index[ mesh.tetrahedra[c][lv] ] ]; }
            tet_color[f][c] = sum * 0.25;
          }
        }
      }

      // per-rank piece file
      {
        std::ostringstream oss; oss << outputDir << "/piece" << rank << ".vtu";
        std::ofstream vtu( oss.str() );
        vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        vtu << "  <UnstructuredGrid>\n";
        vtu << "    <Piece NumberOfPoints=\""<<nb_points<<"\" NumberOfCells=\""<<nb_cells<<"\">\n";
        vtu << "      <Points>\n";
        vtu << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for(const auto& v : mesh.vertices) { vtu << "          " << v.x << " " << v.y << " " << v.z << "\n"; }
        vtu << "        </DataArray>\n";
        vtu << "      </Points>\n";
        vtu << "      <Cells>\n";
        vtu << "        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
        for(const auto& tet : mesh.tetrahedra) { vtu << "          " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << "\n"; }
        vtu << "        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n          ";
        for(size_t c=0;c<nb_cells;c++) { vtu << (4*(c+1)) << " "; }
        vtu << "\n        </DataArray>\n";
        vtu << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
        for(size_t c=0;c<nb_cells;c++) { vtu << "10 "; } // VTK_TETRA
        vtu << "\n        </DataArray>\n";
        vtu << "      </Cells>\n";
        vtu << "      <CellData Scalars=\"rank\">\n";
        vtu << "        <DataArray type=\"Int32\" Name=\"rank\" format=\"ascii\">\n          ";
        for(size_t c=0;c<nb_cells;c++) { vtu << rank << " "; }
        vtu << "\n        </DataArray>\n";
        for(size_t f=0; f<color_fields->size(); f++)
        {
          vtu << "        <DataArray type=\"Float64\" Name=\""<<(*color_fields)[f]<<"\" format=\"ascii\">\n          ";
          for(size_t c=0;c<nb_cells;c++) { vtu << tet_color[f][c] << " "; }
          vtu << "\n        </DataArray>\n";
        }
        if( dxa_tet_classification.has_value() )
        {
          vtu << "        <DataArray type=\"Float64\" Name=\"dxa_good\" format=\"ascii\">\n          ";
          for(size_t c=0;c<nb_cells;c++) { vtu << dxa_tet_classification->good[c] << " "; }
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
        pvtu << "    <PCellData Scalars=\"rank\">\n";
        pvtu << "      <PDataArray type=\"Int32\" Name=\"rank\"/>\n";
        for(const auto& f : *color_fields) { pvtu << "      <PDataArray type=\"Float64\" Name=\""<<f<<"\"/>\n"; }
        if( dxa_tet_classification.has_value() ) { pvtu << "      <PDataArray type=\"Float64\" Name=\"dxa_good\"/>\n"; }
        pvtu << "    </PCellData>\n";
        for(int i=0;i<np;i++) { pvtu << "    <Piece Source=\""<<basename<<"/piece"<<i<<".vtu\"/>\n"; }
        pvtu << "  </PUnstructuredGrid>\n";
        pvtu << "</VTKFile>\n";
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes a DelaunayTessellation (see compute_delaunay) to VTK unstructured grid format for ParaView.
Each MPI rank writes its own piece (<filename>/piece<rank>.vtu); rank 0 additionally writes a
master <filename>.pvtu referencing every piece. Open the .pvtu file in ParaView (not the individual
.vtu pieces) to see the whole system's tessellation across every rank at once, with no gaps or
duplicated cells at rank boundaries.

color_fields (optional list) additionally colors each tetrahedron by any number of per-particle
scalar grid fields (field::mk_generic_real -- e.g. PTM's ptm_type/ptm_rmsd, written by ptm_fields):
each tet gets one CellData array per listed field, the average of its 4 vertices' field value, so
ParaView can switch the active coloring between them without re-running the pipeline.

dxa_tet_classification (optional, auto-wired if compute_dxa_tet_classification ran earlier in the
pipeline) additionally writes DXA step (iv)'s good/bad classification as a direct "dxa_good"
CellData array (1=undistorted patch of the target lattice, 0=defect/grain-boundary/other phase) --
color by this in ParaView to see the "sane" vs "defective" mesh.

Usage example:

compute_delaunay: {}
compute_ptm: { rcut: 3.6 ang }
ptm_fields: {}
compute_dxa_edge_vectors: { target_structure: BCC }
compute_dxa_tet_classification: {}
write_delaunay_vtk: { filename: "paraview/delaunay", color_fields: [ ptm_type, ptm_rmsd ] }
# dxa_tet_classification auto-wires in; color by "dxa_good" in ParaView

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(write_delaunay_vtk)
  {
    OperatorNodeFactory::instance()->register_factory( "write_delaunay_vtk", make_grid_variant_operator< WriteDelaunayVTK > );
  }

}
