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
#include <exanb/core/domain.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/memory/allocator.h>
#include <onika/log.h>

#include <md/snap/snap_config.h>
#include <md/snap/snap_context.h>

#include <cstdint>
#include <string>
#include <vector>
#include <mpi.h>

// Global linear-fitting design matrix for SNAP, built from the output of compute_descriptor_snap
// (per-atom bispectrum + compute_derivative aggregate fields), without recomputing them.
//   rows             : 1 + 3*natoms + 6 (natoms = total atom count, all MPI ranks)
//   columns          : ncoeff * ntypes, one ncoeff-wide block per type (column = ncoeff*itype + k)
//   row 0            : summed bispectrum descriptor
//   rows 1..3*natoms : force-signed gradient w.r.t. atom id's x/y/z (row = 1 + 3*id + xyz),
//                      F = +coeff . row
//   last 6 rows      : sum over slots of r . gradient, Voigt order [xx,yy,zz,yz,xz,xy]
// The sda_* aggregate holds one ncoeff*3 slot per CENTRAL atom type, so the gradient rows loop over
// every type. Each rank accumulates its owned and ghost slots (ghosts weighted by their periodic-image
// position) before the final Allreduce: this operator must run BEFORE update_opt_from_ghost folds
// the ghost aggregates into their owners.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDescriptorSnapGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi  , INPUT , REQUIRED );
    ADD_SLOT( GridT    , grid , INPUT , REQUIRED );
    ADD_SLOT( Domain   , domain , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , bispectrum , INPUT , OPTIONAL , DocString{"see compute_descriptor_snap; required"} );
    ADD_SLOT( long     , ncoeff , INPUT , OPTIONAL , DocString{"see compute_descriptor_snap; required"} );
    ADD_SLOT( md::SnapXSContextRealT<double> , snap_ctx , INPUT , REQUIRED , DocString{"SNAP context built by snap_init (number of types)"} );
    ADD_SLOT( std::string , deriv_agg_field_prefix , INPUT , std::string("sda_")
            , DocString{"Must match compute_descriptor_snap's deriv_agg_field_prefix. Must run before update_opt_from_ghost on these fields."} );

    ADD_SLOT( onika::memory::CudaMMVector<double> , snap_global , OUTPUT
            , DocString{"Row-major (1+3*natoms+6) x (ncoeff*ntypes) global design matrix (natoms = total atom count across the whole simulation, all MPI ranks). Row 0 = summed descriptor; rows 1..3*natoms = per-atom gradient (row=1+3*id+xyz); rows 3*natoms+1..+6 = virial, Voigt order [xx,yy,zz,yz,xz,xy]."} );
    ADD_SLOT( long , ncoeff_all , OUTPUT , DocString{"Number of columns (= ncoeff * ntypes)"} );

  public:
    inline void execute() override final
    {
      if( ! bispectrum.has_value() || ! ncoeff.has_value() )
      {
        fatal_error() << "compute_descriptor_snap_global: bispectrum/ncoeff unavailable -- run compute_descriptor_snap first" << std::endl;
      }

      const long nc = *ncoeff;
      const long ntypes = static_cast<long>( snap_ctx->m_config.materials().size() );
      *ncoeff_all = nc * ntypes;
      const long ncoeff3 = nc * 3;
      const long nc3 = ncoeff3 * ntypes; // widened by ntypes, matches compute_descriptor_snap.cu

      std::vector<const double*> agg_ptr( static_cast<size_t>(nc3), nullptr );
      for( long k=0; k<nc3; k++ )
      {
        agg_ptr[k] = grid->flat_array_data_nocreate( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
        if( agg_ptr[k] == nullptr )
        {
          fatal_error() << "compute_descriptor_snap_global: field '"<<*deriv_agg_field_prefix<<k<<"' not found -- run compute_descriptor_snap with compute_derivative: true first" << std::endl;
        }
      }

      const Mat3d xform = domain->xform();
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

      const long ncols = nc * ntypes;
      const long nrows = 1 + 3*natoms_global + 6;
      const long virial_row0 = 1 + 3*natoms_global;

      snap_global->clear();
      snap_global->resize( static_cast<size_t>(nrows) * static_cast<size_t>(ncols), 0.0 );
      double * const arr = snap_global->data();

      for( size_t ci=0; ci<n_cells; ci++ )
      {
        const bool is_ghost = grid->is_ghost_cell(ci);
        const auto & cell = grid->cell(ci);
        const size_t np = cell.size();
        for( size_t pi=0; pi<np; pi++ )
        {
          const size_t p = cell_particle_offset[ci] + pi;
          const uint64_t id = cell[field::id][pi];
          // real-frame position of this slot (owned atom or ghost image), same frame as the
          // derivative fields, which were computed on xform-applied pair vectors
          const Vec3d r = xform * Vec3d{ cell[field::rx][pi], cell[field::ry][pi], cell[field::rz][pi] };
          const double rx = r.x;
          const double ry = r.y;
          const double rz = r.z;

          if( ! is_ghost )
          {
            // row 0 only needs this atom's own type
            const long itype0 = static_cast<long>( cell[field::type][pi] );
            const double * const src = bispectrum->data() + static_cast<size_t>(nc) * p;
            for( long k=0; k<nc; k++ ) arr[ nc*itype0 + k ] += src[k];
          }

          // gradient rows: the sda_* aggregate has one slot per CENTRAL-atom type, loop over all
          const long grad_row0 = 1 + 3*static_cast<long>(id);
          for( long itype=0; itype<ntypes; itype++ )
          {
            const long typeoffset = ncoeff3 * itype;
            for( long k=0; k<nc; k++ )
            {
              const long col = nc*itype + k;
              const double dx = agg_ptr[typeoffset+k*3+0][p];
              const double dy = agg_ptr[typeoffset+k*3+1][p];
              const double dz = agg_ptr[typeoffset+k*3+2][p];

              arr[ (grad_row0+0)*ncols + col ] += dx;
              arr[ (grad_row0+1)*ncols + col ] += dy;
              arr[ (grad_row0+2)*ncols + col ] += dz;

              arr[ (virial_row0+0)*ncols + col ] += dx*rx; // xx
              arr[ (virial_row0+1)*ncols + col ] += dy*ry; // yy
              arr[ (virial_row0+2)*ncols + col ] += dz*rz; // zz
              arr[ (virial_row0+3)*ncols + col ] += dz*ry; // yz
              arr[ (virial_row0+4)*ncols + col ] += dz*rx; // xz
              arr[ (virial_row0+5)*ncols + col ] += dy*rx; // xy
            }
          }
        }
      }

      MPI_Allreduce( MPI_IN_PLACE, arr, static_cast<int>(nrows*ncols), MPI_DOUBLE, MPI_SUM, *mpi );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Global linear-fitting design matrix for SNAP. Row-major (1+3*natoms+6) x (ncoeff*ntypes) array,
one ncoeff-wide column block per type (column = ncoeff*itype + k):

  row 0             -- summed bispectrum descriptor. coeff . row = total energy.
  rows 1..3*natoms  -- force-signed gradient w.r.t. atom m's x/y/z, at row 1+3*m+xyz
                       (m = particle id, 0-indexed). coeff . row = force component.
  rows 3N+1..3N+6   -- virial, Voigt order [xx,yy,zz,yz,xz,xy]. coeff . row = virial component.

The matrix has no reference-label column. Must run after compute_descriptor_snap
(compute_derivative: true) and before any update_opt_from_ghost on its aggregate fields.

Usage example:

init_parameters:
  - species
  - snap_init: { parameters: { param: "W.snapparam", coef: "W.snapcoeff" } }

compute_descriptor_snap: { compute_derivative: true }
compute_descriptor_snap_global
write_descriptor_snap_global: { filename: "snap_global.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_descriptor_snap_global)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_descriptor_snap_global", make_grid_variant_operator< ComputeDescriptorSnapGlobal > );
  }

}
