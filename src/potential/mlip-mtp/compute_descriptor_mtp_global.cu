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
#include <onika/log.h>

#include <cstdint>
#include <string>
#include <vector>
#include <mpi.h>

#include "include/mtp_config.h"

// Global linear-fitting design matrix for MTP -- same two-stage local-then-Allreduce design as
// snap/compute_descriptor_snap_global.cu and k2b/compute_descriptor_k2b_global.cu (NOT POD's
// single-rank, re-run-from-scratch pod_global_op.h pattern). Purely additive: does NOT re-run the
// descriptor+derivative pass -- reads compute_descriptor_mtp's existing output (the per-atom
// `mtp_descriptors` buffer and its `compute_derivative: true` aggregate fields).
//
// Column layout differs from POD's own per-element-block convention, and this IS the correct
// shape for MTP (not a simplification): MTP's fitted `linear_coeffs` (moment_coeffs in the file)
// is a SINGLE array shared across every species -- only `species_coeffs` (one scalar per species)
// varies by type. So:
//   size_array_cols = species_count + alpha_scalar_moments
//   columns [0, species_count)                -- one-hot species-count indicators (row 0 only,
//                                                zero gradient/virial contribution)
//   columns [species_count, species_count+alpha_scalar_moments) -- the shared B_k columns
// Dotting row 0 with [species_coeffs..., linear_coeffs...] gives the total configuration energy
// directly.
//
// Row layout (rows), matching the established convention (no reference-label column):
//   size_array_rows = 1 + 3*natoms + 6
//   row 0            -- summed one-hot species counts + summed MTP descriptor over every atom
//   rows 1..3*natoms -- per-atom true dB_k/dR_m (self term + every neighbor interaction),
//                       row = 1 + 3*field::id + xyz (exaStamp field::id is 0-indexed, no "-1")
//   rows 3N+1..3N+6  -- summed r_atom . dB_atom/dR_atom, Voigt order [xx,yy,zz,yz,xz,xy]
//
// Real ordering requirement: must run AFTER compute_descriptor_mtp: { compute_derivative: true }
// and BEFORE any update_opt_from_ghost call on its aggregate fields -- needs the raw, per-rank-
// local, UN-FOLDED aggregate. See compute_descriptor_snap_global.cu's header comment for the full
// reasoning (identical here).
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDescriptorMtpGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi     , INPUT , REQUIRED );
    ADD_SLOT( GridT    , grid    , INPUT , REQUIRED );
    ADD_SLOT( MtpContext , mtp_ctx , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , mtp_descriptors , INPUT , OPTIONAL , DocString{"see compute_descriptor_mtp; required"} );
    ADD_SLOT( long     , ncoeff , INPUT , OPTIONAL , DocString{"see compute_descriptor_mtp; required"} );
    ADD_SLOT( std::string , deriv_agg_field_prefix , INPUT , std::string("mda_")
            , DocString{"Must match compute_descriptor_mtp's own deriv_agg_field_prefix. Must be read before update_opt_from_ghost runs on these fields -- see this file's header comment."} );

    ADD_SLOT( onika::memory::CudaMMVector<double> , mtp_global , OUTPUT
            , DocString{"Row-major (1+3*natoms+6) x (species_count+ncoeff) global design matrix. Row 0 = summed one-hot species counts + summed descriptor; rows 1..3*natoms = per-atom gradient of the ncoeff B_k columns only (row=1+3*id+xyz); rows 3*natoms+1..+6 = virial, Voigt order [xx,yy,zz,yz,xz,xy]. See this file's header comment for why the column layout differs from POD's."} );
    ADD_SLOT( long , ncoeff_all , OUTPUT , DocString{"Number of columns (= species_count + ncoeff)"} );

  public:
    inline void execute() override final
    {
      if( ! mtp_descriptors.has_value() || ! ncoeff.has_value() )
      {
        fatal_error() << "compute_descriptor_mtp_global: mtp_descriptors/ncoeff unavailable -- run compute_descriptor_mtp first" << std::endl;
      }

      const long nc = *ncoeff;
      const int species_count = mtp_ctx->nspecies;
      const long nc_all = species_count + nc;
      *ncoeff_all = nc_all;
      const long nc3 = nc * 3;

      std::vector<const double*> agg_ptr( static_cast<size_t>(nc3), nullptr );
      for( long k=0; k<nc3; k++ )
      {
        agg_ptr[k] = grid->flat_array_data_nocreate( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
        if( agg_ptr[k] == nullptr )
        {
          fatal_error() << "compute_descriptor_mtp_global: field '"<<*deriv_agg_field_prefix<<k<<"' not found -- run compute_descriptor_mtp with compute_derivative: true first" << std::endl;
        }
      }

      const auto & type_map = mtp_ctx->type_map;
      const auto * cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells = grid->number_of_cells();

      // local owned (non-ghost) atom count -> global total via one small Allreduce, so every
      // rank allocates the same full-size array before the main per-atom Allreduce below
      long local_owned = 0;
      for( size_t ci=0; ci<n_cells; ci++ )
      {
        if( grid->is_ghost_cell(ci) ) continue;
        local_owned += static_cast<long>( grid->cell(ci).size() );
      }
      long natoms_global = 0;
      MPI_Allreduce( &local_owned, &natoms_global, 1, MPI_LONG, MPI_SUM, *mpi );

      const long nrows = 1 + 3*natoms_global + 6;
      const long virial_row0 = 1 + 3*natoms_global;

      mtp_global->clear();
      mtp_global->resize( static_cast<size_t>(nrows) * static_cast<size_t>(nc_all), 0.0 );
      double * const arr = mtp_global->data();

      for( size_t ci=0; ci<n_cells; ci++ )
      {
        const bool is_ghost = grid->is_ghost_cell(ci);
        const auto & cell = grid->cell(ci);
        const size_t np = cell.size();
        for( size_t pi=0; pi<np; pi++ )
        {
          const size_t p = cell_particle_offset[ci] + pi;
          const uint64_t id = cell[field::id][pi];
          const double rx = cell[field::rx][pi];
          const double ry = cell[field::ry][pi];
          const double rz = cell[field::rz][pi];

          if( ! is_ghost )
          {
            const int itype = cell[field::type][pi];
            const int mtp_sp = type_map[itype];
            arr[ mtp_sp ] += 1.0;   // one-hot species-count column, row 0
            const double * const src = mtp_descriptors->data() + static_cast<size_t>(nc) * p;
            for( long k=0; k<nc; k++ ) arr[species_count + k] += src[k];
          }

          const long grad_row0 = 1 + 3*static_cast<long>(id);
          for( long k=0; k<nc; k++ )
          {
            const double dx = agg_ptr[k*3+0][p];
            const double dy = agg_ptr[k*3+1][p];
            const double dz = agg_ptr[k*3+2][p];
            const long col = species_count + k;

            arr[ (grad_row0+0)*nc_all + col ] += dx;
            arr[ (grad_row0+1)*nc_all + col ] += dy;
            arr[ (grad_row0+2)*nc_all + col ] += dz;

            arr[ (virial_row0+0)*nc_all + col ] += dx*rx; // xx
            arr[ (virial_row0+1)*nc_all + col ] += dy*ry; // yy
            arr[ (virial_row0+2)*nc_all + col ] += dz*rz; // zz
            arr[ (virial_row0+3)*nc_all + col ] += dz*ry; // yz
            arr[ (virial_row0+4)*nc_all + col ] += dz*rx; // xz
            arr[ (virial_row0+5)*nc_all + col ] += dy*rx; // xy
          }
        }
      }

      MPI_Allreduce( MPI_IN_PLACE, arr, static_cast<int>(nrows*nc_all), MPI_DOUBLE, MPI_SUM, *mpi );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Global linear-fitting design matrix for MTP -- same row/rank convention as
compute_descriptor_snap_global, but a different COLUMN layout because MTP's fitted coefficients
are structured differently from POD/SNAP: `linear_coeffs` is a single array shared across every
species (only a per-species scalar offset varies by type), so:

  columns [0, species_count)                         -- one-hot species-count indicators (row 0 only)
  columns [species_count, species_count+ncoeff)       -- the shared per-basis-function B_k columns

Row 0 dotted with [species_coeffs..., linear_coeffs...] gives the total configuration energy
directly. Rows 1..3*natoms give the true dB_k/dR_m gradient (self term + every neighbor
interaction) at row 1+3*id+xyz (id = field::id, 0-indexed) -- only the trailing ncoeff columns are
populated (the species columns have zero spatial gradient). Rows 3N+1..3N+6 give the virial,
Voigt order [xx,yy,zz,yz,xz,xy], same trailing-columns-only rule.

Purely additive: reads compute_descriptor_mtp's existing output rather than re-running the
descriptor+derivative pass. Must run right after compute_descriptor_mtp (compute_derivative: true)
and BEFORE any update_opt_from_ghost call on its aggregate fields -- this operator needs the raw,
un-folded per-rank-local aggregate.

Usage example:

compute_descriptor_mtp: { compute_derivative: true }
compute_descriptor_mtp_global
update_opt_from_ghost: { opt_fields: [ "mda_.*" ] }   # only needed if the per-atom export below also runs
write_descriptor_mtp_global: { filename: "mtp_global.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_descriptor_mtp_global)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_descriptor_mtp_global", make_grid_variant_operator< ComputeDescriptorMtpGlobal > );
  }

}
