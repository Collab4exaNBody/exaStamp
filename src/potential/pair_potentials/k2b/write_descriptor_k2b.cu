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
#include <onika/log.h>

// header-only, no link dependency (not gated by EXASTAMP_MLIP_POD_BUILD) -- k2b's real build path is
// src/potential/pair_potentials/k2b/ (no symlinking, see AddPairPotential in
// pair_potential_template/CMakeLists.txt), so this relative include resolves with no CMake change
#include "../../mlip-pod/include/npy_writer.h"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <string>
#include <vector>
#include <mpi.h>

// Single combined, MPI-gathered per-atom export of compute_descriptor_k2b's outputs -- direct
// structural copy of mlip-pod/write_descriptor_pod.cu (same hand-rolled MPI_Gather+MPI_Gatherv
// pattern). Regardless of how many MPI ranks own the simulation, this writes ONE file (from
// rank 0) with one row per globally-owned (non-ghost) atom, user-selectable columns, in either
// 'text' (default) or 'npy' (single combined float64 array, same convention as write_descriptor_pod
// and write_descriptor_snap) format.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class WriteDescriptorK2b : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi  , INPUT , REQUIRED );
    ADD_SLOT( GridT    , grid , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , k2b_descriptors , INPUT , OPTIONAL , DocString{"see compute_descriptor_k2b; required if 'descriptor' is in fields"} );
    ADD_SLOT( long     , ncoeff , INPUT , OPTIONAL , DocString{"see compute_descriptor_k2b; required if 'descriptor' or 'derivative' is in fields"} );
    ADD_SLOT( std::string , deriv_agg_field_prefix , INPUT , std::string("k2bda_")
            , DocString{"Must match compute_descriptor_k2b's own deriv_agg_field_prefix; required if 'derivative' is in fields. Run update_opt_from_ghost on these fields first on a multi-rank run -- see compute_descriptor_k2b's documentation."} );
    ADD_SLOT( std::vector<std::string> , fields , INPUT
            , std::vector<std::string>{"id","x","y","z","descriptor","derivative"}
            , DocString{"Columns to write, in this order. Choices: id, x, y, z, descriptor (ncoeff values), derivative (ncoeff*3 values -- the compact per-atom descriptor-derivative aggregate, NOT the full per-neighbor-pair Jacobian)."} );
    ADD_SLOT( std::string , filename , INPUT , std::string("k2b_descriptors.txt") , DocString{"Single combined output file, written once from rank 0 after gathering every rank's owned atoms. In 'npy' format this is used as a prefix (a trailing '.txt' is stripped) for the single combined '<prefix>.npy' file."} );
    ADD_SLOT( std::string , format , INPUT , std::string("text") , DocString{"Output format: 'text' (default, single combined plain-text file) or 'npy' (single combined '<prefix>.npy' file, one row per atom, columns per `fields` -- same convention as write_descriptor_pod/write_descriptor_snap), directly loadable via numpy.load())."} );

  public:
    inline void execute() override final
    {
      if( *format != "text" && *format != "npy" )
      {
        fatal_error() << "write_descriptor_k2b: unknown format '"<<*format<<"' (choices: text, npy)" << std::endl;
      }

      bool want_id=false, want_x=false, want_y=false, want_z=false, want_desc=false, want_deriv=false;
      for( const auto & f : *fields )
      {
             if( f == "id" )         want_id    = true;
        else if( f == "x" )          want_x     = true;
        else if( f == "y" )          want_y     = true;
        else if( f == "z" )          want_z     = true;
        else if( f == "descriptor" ) want_desc  = true;
        else if( f == "derivative" ) want_deriv = true;
        else fatal_error() << "write_descriptor_k2b: unknown field '"<<f<<"' (choices: id, x, y, z, descriptor, derivative)" << std::endl;
      }
      if( want_desc && ( ! k2b_descriptors.has_value() || ! ncoeff.has_value() ) )
      {
        fatal_error() << "write_descriptor_k2b: 'descriptor' requested but k2b_descriptors/ncoeff unavailable -- run compute_descriptor_k2b first" << std::endl;
      }
      if( want_deriv && ! ncoeff.has_value() )
      {
        fatal_error() << "write_descriptor_k2b: 'derivative' requested but ncoeff unavailable -- run compute_descriptor_k2b with compute_derivative: true first" << std::endl;
      }

      const long nc = ncoeff.has_value() ? *ncoeff : 0;

      std::vector<const double*> deriv_comp_ptr;
      if( want_deriv )
      {
        deriv_comp_ptr.resize( nc*3 );
        for( long k=0; k<nc*3; k++ )
        {
          deriv_comp_ptr[k] = grid->flat_array_data_nocreate( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
          if( deriv_comp_ptr[k] == nullptr )
          {
            fatal_error() << "write_descriptor_k2b: field '"<<*deriv_agg_field_prefix<<k<<"' not found -- run compute_descriptor_k2b with compute_derivative: true first" << std::endl;
          }
        }
      }

      // fixed internal per-atom payload layout (independent of the user-facing column order)
      int off_x=-1, off_y=-1, off_z=-1, off_desc=-1, off_deriv=-1, payload_width=0;
      if( want_x    ) { off_x    = payload_width; payload_width += 1; }
      if( want_y    ) { off_y    = payload_width; payload_width += 1; }
      if( want_z    ) { off_z    = payload_width; payload_width += 1; }
      if( want_desc ) { off_desc = payload_width; payload_width += nc; }
      if( want_deriv) { off_deriv= payload_width; payload_width += nc*3; }

      std::vector<uint64_t> local_ids;
      std::vector<double> local_payload;
      const auto * cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells = grid->number_of_cells();
      for( size_t ci=0; ci<n_cells; ci++ )
      {
        if( grid->is_ghost_cell(ci) ) continue;
        const auto & cell = grid->cell(ci);
        const size_t np = cell.size();
        for( size_t pi=0; pi<np; pi++ )
        {
          if( want_id ) local_ids.push_back( cell[field::id][pi] );
          if( payload_width == 0 ) continue;
          const size_t p = cell_particle_offset[ci] + pi;
          if( want_x ) local_payload.push_back( cell[field::rx][pi] );
          if( want_y ) local_payload.push_back( cell[field::ry][pi] );
          if( want_z ) local_payload.push_back( cell[field::rz][pi] );
          if( want_desc )
          {
            const double * const src = k2b_descriptors->data() + static_cast<size_t>(nc) * p;
            for( long k=0; k<nc; k++ ) local_payload.push_back( src[k] );
          }
          if( want_deriv )
          {
            for( long k=0; k<nc*3; k++ ) local_payload.push_back( deriv_comp_ptr[k][p] );
          }
        }
      }

      int rank=0, nprocs=1;
      MPI_Comm_rank( *mpi, &rank );
      MPI_Comm_size( *mpi, &nprocs );

      const int local_count = want_id ? static_cast<int>(local_ids.size())
                             : ( payload_width>0 ? static_cast<int>(local_payload.size()/payload_width) : 0 );
      std::vector<int> counts( nprocs, 0 );
      MPI_Gather( &local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, *mpi );

      std::vector<uint64_t> all_ids;
      if( want_id )
      {
        std::vector<int> id_counts, id_displs;
        if( rank == 0 )
        {
          id_counts = counts;
          id_displs.resize( nprocs );
          int total = 0;
          for( int r=0; r<nprocs; r++ ) { id_displs[r] = total; total += id_counts[r]; }
          all_ids.resize( total );
        }
        MPI_Gatherv( local_ids.data(), local_count, MPI_UNSIGNED_LONG_LONG
                   , rank==0 ? all_ids.data()   : nullptr
                   , rank==0 ? id_counts.data() : nullptr
                   , rank==0 ? id_displs.data() : nullptr
                   , MPI_UNSIGNED_LONG_LONG, 0, *mpi );
      }

      std::vector<double> all_payload;
      if( payload_width > 0 )
      {
        std::vector<int> pl_counts, pl_displs;
        if( rank == 0 )
        {
          pl_counts.resize( nprocs );
          pl_displs.resize( nprocs );
          int total = 0;
          for( int r=0; r<nprocs; r++ ) { pl_counts[r] = counts[r]*payload_width; pl_displs[r] = total; total += pl_counts[r]; }
          all_payload.resize( total );
        }
        MPI_Gatherv( local_payload.data(), local_count*payload_width, MPI_DOUBLE
                   , rank==0 ? all_payload.data() : nullptr
                   , rank==0 ? pl_counts.data()   : nullptr
                   , rank==0 ? pl_displs.data()   : nullptr
                   , MPI_DOUBLE, 0, *mpi );
      }

      if( rank != 0 ) return;

      const long total_atoms = std::accumulate( counts.begin(), counts.end(), 0L );

      if( *format == "npy" )
      {
        std::string prefix = *filename;
        static constexpr const char * TXT_SUFFIX = ".txt";
        if( prefix.size() >= 4 && prefix.compare(prefix.size()-4, 4, TXT_SUFFIX) == 0 ) prefix.resize(prefix.size()-4);

        // one combined 2-D float64 array, columns in the same order as `fields` (id included as
        // float64 -- exact for any realistic atom count) -- same convention as write_descriptor_pod
        // and write_descriptor_snap, not one file per field.
        long total_width = 0;
        for( const auto & f : *fields )
        {
               if( f == "id" || f == "x" || f == "y" || f == "z" ) total_width += 1;
          else if( f == "descriptor" ) total_width += nc;
          else if( f == "derivative" ) total_width += nc*3;
        }

        std::vector<double> buf( static_cast<size_t>(total_atoms) * total_width );
        for( long a=0; a<total_atoms; a++ )
        {
          double * const out = buf.data() + static_cast<size_t>(a)*total_width;
          long c = 0;
          for( const auto & f : *fields )
          {
                 if( f == "id" ) out[c++] = static_cast<double>( all_ids[a] );
            else if( f == "x"  ) out[c++] = all_payload[static_cast<size_t>(a)*payload_width + off_x];
            else if( f == "y"  ) out[c++] = all_payload[static_cast<size_t>(a)*payload_width + off_y];
            else if( f == "z"  ) out[c++] = all_payload[static_cast<size_t>(a)*payload_width + off_z];
            else if( f == "descriptor" ) { for( long k=0; k<nc;   k++ ) out[c++] = all_payload[static_cast<size_t>(a)*payload_width + off_desc + k]; }
            else if( f == "derivative" ) { for( long k=0; k<nc*3; k++ ) out[c++] = all_payload[static_cast<size_t>(a)*payload_width + off_deriv + k]; }
          }
        }
        write_npy( onika::data_file_path(prefix+".npy"), {static_cast<size_t>(total_atoms), static_cast<size_t>(total_width)}, "<f8", buf.data(), sizeof(double) );
        return;
      }

      std::ofstream fout( onika::data_file_path(*filename) );
      fout << std::setprecision(17);
      for( long a=0; a<total_atoms; a++ )
      {
        const double * const prow = ( payload_width>0 ) ? ( all_payload.data() + static_cast<size_t>(a)*payload_width ) : nullptr;
        bool first = true;
        for( const auto & f : *fields )
        {
          if( !first ) fout << " ";
          first = false;
               if( f == "id" ) fout << all_ids[a];
          else if( f == "x"  ) fout << prow[off_x];
          else if( f == "y"  ) fout << prow[off_y];
          else if( f == "z"  ) fout << prow[off_z];
          else if( f == "descriptor" ) { for( long k=0; k<nc;   k++ ) { if(k>0) fout << " "; fout << prow[off_desc+k]; } }
          else if( f == "derivative" ) { for( long k=0; k<nc*3; k++ ) { if(k>0) fout << " "; fout << prow[off_deriv+k]; } }
        }
        fout << "\n";
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes ONE combined plain-text file (from rank 0, after MPI-gathering every rank's owned,
non-ghost atoms -- independent of the number of MPI processes) with one row per atom and
user-selectable columns, in the order given by `fields`:

  - id           : field::id
  - x, y, z      : position
  - descriptor   : compute_descriptor_k2b's per-atom 2-body kernel descriptor vector, ncoeff
                   (= n_rbf) values
  - derivative   : compute_descriptor_k2b's compact per-atom descriptor-derivative aggregate,
                   ncoeff*3 values (self term + sum over every atom that has this one as a
                   neighbor), NOT the full per-neighbor-pair Jacobian. Read from the
                   deriv_agg_field_prefix-named dynamic fields -- run update_opt_from_ghost
                   on them first on a multi-rank run (see compute_descriptor_k2b's doc).

Setting `format: npy` writes a single `.npy` v1.0 file, `<prefix>.npy` (`<prefix>` is `filename`
with a trailing '.txt' stripped if present) -- a float64 array of shape (natoms, total_width), one
row per atom, columns in the same order as `fields` (id: 1 column, x/y/z: 1 column each,
descriptor: ncoeff columns, derivative: ncoeff*3 columns) -- same convention as
write_descriptor_pod and write_descriptor_snap. `id` is stored as float64 in this combined array
(exact for any realistic atom count). Directly loadable with `numpy.load()`, no custom parser needed.

Usage example:

compute_descriptor_k2b: { compute_derivative: true }
update_opt_from_ghost: { opt_fields: [ "k2bda_.*" ] }   # only needed on more than one rank
write_descriptor_k2b:
  filename: "descriptors.txt"
  fields: [ id, x, y, z, descriptor, derivative ]

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(write_descriptor_k2b)
  {
    OperatorNodeFactory::instance()->register_factory( "write_descriptor_k2b", make_grid_variant_operator< WriteDescriptorK2b > );
  }

}
