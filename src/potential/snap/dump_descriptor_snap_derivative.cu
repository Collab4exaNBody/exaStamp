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

// Plain-text export of compute_descriptor_snap's compute_derivative CSR output, for
// cross-validation against an external reference (LAMMPS compute snad/atom, after
// reproducing its self+neighbor aggregation -- see
// data/regression_new/analysis_particle/test_snap_descriptors/compare_derivatives.py).
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class DumpDescriptorSnapDerivative : public OperatorNode
  {
    ADD_SLOT( GridT   , grid                    , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<long>    , bispectrum_deriv_offset , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double>  , bispectrum_deriv        , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<uint64_t>, bispectrum_deriv_nbh_id , INPUT , REQUIRED );
    ADD_SLOT( long    , ncoeff     , INPUT , REQUIRED );
    ADD_SLOT( std::string , filename , INPUT , std::string("descriptor_derivatives.txt") );

  public:
    inline void execute() override final
    {
      const auto * cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells = grid->number_of_cells();
      const long nc = *ncoeff;
      const long nc3 = nc * 3;

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
          const size_t p = cell_particle_offset[ci] + pi;
          const long row0 = (*bispectrum_deriv_offset)[p];
          const long row1 = (*bispectrum_deriv_offset)[p+1];
          fout << "ATOM " << id << " " << (row1-row0) << "\n";
          for( long row=row0; row<row1; row++ )
          {
            fout << (*bispectrum_deriv_nbh_id)[row];
            const double * const vec = bispectrum_deriv->data() + row*nc3;
            for( long i=0; i<nc3; i++ ) fout << " " << vec[i];
            fout << "\n";
          }
        }
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes compute_descriptor_snap's compute_derivative CSR output to a plain-text file, one
block per local particle:

  ATOM <id> <nrows>
  <neighbor_id> <dB0/dxJ> <dB0/dyJ> <dB0/dzJ> <dB1/dxJ> ... <dB_{ncoeff-1}/dzJ>
  ... (nrows lines, row 0 = the atom's own self/negative-sum term, neighbor_id == id)

Not MPI-gathered: on a multi-rank run each rank writes only its own local particles.

Usage example:

compute_descriptor_snap: { compute_derivative: true, parameters: { param: "W.snapparam", coef: "W.snapcoeff" } }
dump_descriptor_snap_derivative: { filename: "descriptor_derivatives.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(dump_descriptor_snap_derivative)
  {
    OperatorNodeFactory::instance()->register_factory( "dump_descriptor_snap_derivative", make_grid_variant_operator< DumpDescriptorSnapDerivative > );
  }

}
