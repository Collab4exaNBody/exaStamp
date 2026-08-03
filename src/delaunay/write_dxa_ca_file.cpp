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

#include <fstream>
#include <string>

// Writes this operator's DXADislocationLines to OVITO's own "Crystal Analysis" (.ca) file format,
// so the result can be loaded directly into OVITO for visual comparison against its own DXA output
// -- reverse-engineered from OVITO's real (non-public) exporter/importer source, obtained directly
// by the user (CAExporter.cpp/CAImporter.cpp), not from public documentation.
//
// Deliberately minimal/simplified relative to a real OVITO-written file:
//   - Only one STRUCTURE_TYPE is declared (a bare "bcc" stub with one placeholder Burgers vector
//     family) -- enough for the importer's header parsing to succeed, not a faithful reproduction
//     of OVITO's own per-structure Burgers vector family tables.
//   - DISLOCATION_JUNCTIONS: this operator doesn't reconstruct real multi-way junction connectivity
//     (see compute_dxa_circuit_sweep.cpp's own "ponytail" note), so every dislocation is written as
//     a trivial self-referential 2-cycle (both of its own ends point back to each other) rather than
//     real cross-dislocation junction rings. OVITO will render every line as an independent,
//     unconnected segment -- exactly what this operator actually knows, no more.
//   - Single-rank only: writes whatever DXADislocationLines this rank holds directly, no cross-rank
//     gather/merge (unlike write_dxa_dislocation_lines.cpp's per-rank-piece + .pvtu convention --
//     the .ca format has no native multi-piece convention to hang that off of). Fine for the
//     single-MPI-rank validation runs this operator was built for; would need real gathering for a
//     genuine multi-rank use case.
//   - No DEFECT_MESH section (OVITO's own interface-mesh triangulation) -- write_interface_mesh.cpp
//     already covers that separately, in VTK format.
namespace exaStamp
{
  using namespace exanb;

  class WriteDXACAFile : public OperatorNode
  {
    ADD_SLOT( Domain               , domain                 , INPUT , REQUIRED );
    ADD_SLOT( DXADislocationLines  , dxa_dislocation_lines  , INPUT , REQUIRED );
    ADD_SLOT( std::string          , filename               , INPUT , std::string("dislocation_lines.ca") );

  public:
    inline void execute () override final
    {
      const DXADislocationLines& dl = *dxa_dislocation_lines;
      const size_t n = dl.line_positions.empty() ? dl.lines.size() : dl.line_positions.size();

      std::ofstream f( *filename );
      f << "CA_FILE_VERSION 6\n";
      f << "CA_LIB_VERSION 0.0.0\n";
      f << "STRUCTURE_TYPES 1\n";
      f << "STRUCTURE_TYPE 1\n";
      f << "NAME bcc\n";
      f << "FULL_NAME BCC\n";
      f << "COLOR 0.4 0.4 1.0\n";
      f << "TYPE LATTICE\n";
      f << "BURGERS_VECTOR_FAMILIES 1\n";
      f << "BURGERS_VECTOR_FAMILY ID 0\n";
      f << "Other\n";
      f << "0.0 0.0 0.0\n";
      f << "0.9 0.2 0.2\n";
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
      f << "CLUSTER_STRUCTURE 1\n";
      f << "CLUSTER_ORIENTATION\n";
      f << "1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n";
      f << "CLUSTER_COLOR 1.0 1.0 1.0\n";
      f << "CLUSTER_SIZE 0\n";
      f << "END_CLUSTER\n";
      f << "CLUSTER_TRANSITIONS 0\n";

      f << "DISLOCATIONS " << n << "\n";
      for(size_t i=0;i<n;i++)
      {
        const Vec3d& b = dl.burgers_vector[i];
        f << i << "\n";
        f << b.x << " " << b.y << " " << b.z << "\n";
        f << "1\n";
        if( !dl.line_positions.empty() )
        {
          const auto& pts = dl.line_positions[i];
          f << pts.size() << "\n";
          for(const Vec3d& p : pts) { f << p.x << " " << p.y << " " << p.z << " 0\n"; }
        }
        else
        {
          f << "0\n";
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
against OVITO's own DXA output. See this file's own header comment for the specific simplifications
(minimal structure-type stub, trivial self-referential junction records, single-rank only).

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
