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

#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/memory/allocator.h>
#include <onika/file_utils.h>

#include <fstream>
#include <iomanip>
#include <string>

// Plain-text export of compute_descriptor_pod's per-particle descriptor buffer, keyed by
// atom id and position, for cross-validation against an external reference (e.g. LAMMPS
// compute pod/atom) -- see data/regression_new/analysis_particle/test_pod_descriptors.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class DumpDescriptorPod : public OperatorNode
  {
    ADD_SLOT( GridT   , grid            , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , pod_descriptors , INPUT , REQUIRED , DocString{"Flat per-particle descriptor buffer (see compute_descriptor_pod)"} );
    ADD_SLOT( long    , ncoeff          , INPUT , REQUIRED , DocString{"Number of descriptor components per particle (see compute_descriptor_pod)"} );
    ADD_SLOT( std::string , filename    , INPUT , std::string("pod_descriptors.txt") , DocString{"Output file: one line per local particle, 'id x y z d0 d1 ... d_{ncoeff-1}'"} );

  public:
    inline void execute() override final
    {
      const auto * cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells = grid->number_of_cells();
      const long nc = *ncoeff;

      std::ofstream fout( onika::data_file_path(*filename) );
      fout << std::setprecision(17);
      for( size_t ci=0; ci<n_cells; ci++ )
      {
        if( grid->is_ghost_cell(ci) ) continue;
        const auto & cell = grid->cell(ci);
        const size_t n_particles = cell.size();
        for( size_t pi=0; pi<n_particles; pi++ )
        {
          const uint64_t id = cell[field::id][pi];
          const double x = cell[field::rx][pi];
          const double y = cell[field::ry][pi];
          const double z = cell[field::rz][pi];
          const size_t off = static_cast<size_t>(nc) * ( cell_particle_offset[ci] + pi );
          fout << id << " " << x << " " << y << " " << z;
          for( long k=0; k<nc; k++ ) fout << " " << (*pod_descriptors)[off+k];
          fout << "\n";
        }
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes every local (non-ghost) particle's id, position, and per-atom POD descriptor
(computed by compute_descriptor_pod) to a plain-text file, one atom per line:

  id x y z d0 d1 ... d_{ncoeff-1}

Not MPI-gathered: on a multi-rank run each rank writes only its own local particles --
use a distinct filename per rank, or run single-rank for a full-system dump.

Usage example:

pod_init: { parameters: { pod_file: "Ta_param.pod", coeff_file: "Ta_coefficients.pod" } }
compute_descriptor_pod
dump_descriptor_pod: { filename: "descriptors.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(dump_descriptor_pod)
  {
    OperatorNodeFactory::instance()->register_factory( "dump_descriptor_pod", make_grid_variant_operator< DumpDescriptorPod > );
  }

}
