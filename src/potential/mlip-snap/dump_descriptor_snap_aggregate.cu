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
#include <vector>

// Plain-text export of compute_descriptor_snap's LAMMPS compute-snad/atom-equivalent
// aggregate (compute_derivative: true), written in LAMMPS's own column order
// ([xyz-block][coeff], not our internal [coeff][xyz]) so the file is directly diffable
// against a `dump ... c_snad[*]` file, no reshaping needed -- see
// data/regression_new/analysis_particle/test_snap_descriptors/compare_aggregate.py.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class DumpDescriptorSnapAggregate : public OperatorNode
  {
    ADD_SLOT( GridT   , grid       , INPUT , REQUIRED );
    ADD_SLOT( long    , ncoeff     , INPUT , REQUIRED );
    ADD_SLOT( std::string , deriv_agg_field_prefix , INPUT , std::string("sda_")
            , DocString{"Must match compute_descriptor_snap's own deriv_agg_field_prefix. Run update_opt_from_ghost on these fields first on a multi-rank run -- see compute_descriptor_snap's documentation."} );
    ADD_SLOT( std::string , filename , INPUT , std::string("descriptor_derivative_aggregate.txt") );

  public:
    inline void execute() override final
    {
      const auto * cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells = grid->number_of_cells();
      const long nc = *ncoeff;
      const long nc3 = nc * 3;

      std::vector<const double*> comp_ptr( nc3 );
      for( long k=0; k<nc3; k++ )
      {
        comp_ptr[k] = grid->flat_array_data_nocreate( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
        if( comp_ptr[k] == nullptr )
        {
          fatal_error() << "dump_descriptor_snap_aggregate: field '"<<*deriv_agg_field_prefix<<k<<"' not found -- run compute_descriptor_snap with compute_derivative: true first" << std::endl;
        }
      }

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
          const size_t p = cell_particle_offset[ci] + pi;
          fout << id << " " << x << " " << y << " " << z;
          // reorder our [k*3+xyz] layout into LAMMPS compute snad/atom's [xyz*ncoeff+k]
          for( int c=0; c<3; c++ )
            for( long k=0; k<nc; k++ )
              fout << " " << comp_ptr[k*3+c][p];
          fout << "\n";
        }
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes compute_descriptor_snap's LAMMPS compute-snad/atom-equivalent aggregate
(compute_derivative: true) to a plain-text file, one line per local particle, in the SAME
column order as LAMMPS `compute snad/atom` (dBx for k=0..ncoeff-1, then dBy, then dBz):

  id x y z dB0/dx dB1/dx ... dB_{ncoeff-1}/dx  dB0/dy ... dB_{ncoeff-1}/dy  dB0/dz ... dB_{ncoeff-1}/dz

The aggregate lives in ncoeff*3 dynamically-named generic-real grid fields (see
compute_descriptor_snap's deriv_agg_field_prefix) rather than a private buffer, so it can be
reduced across MPI ranks via the generic update_opt_from_ghost operator -- run that (with
opt_fields matching the same prefix) right after compute_descriptor_snap, before this
operator, whenever more than one rank is used. Not itself MPI-gathered: each rank writes only
its own local particles.

Usage example:

compute_descriptor_snap: { compute_derivative: true, parameters: { param: "W.snapparam", coef: "W.snapcoeff" } }
update_opt_from_ghost: { opt_fields: [ "sda_.*" ] }   # only needed on more than one rank
dump_descriptor_snap_aggregate: { filename: "exastamp_snad.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(dump_descriptor_snap_aggregate)
  {
    OperatorNodeFactory::instance()->register_factory( "dump_descriptor_snap_aggregate", make_grid_variant_operator< DumpDescriptorSnapAggregate > );
  }

}
