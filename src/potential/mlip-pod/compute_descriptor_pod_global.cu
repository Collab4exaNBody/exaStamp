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

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>

#include <cstdint>
#include <string>
#include <vector>
#include <mpi.h>

#include "pod_params.h"
#include "pod_config.h"

// Global linear-fitting design matrix for POD -- consumes compute_descriptor_pod's already-computed
// output (pod_descriptors + compute_derivative: true's pda_* aggregate) instead of re-running its
// own independent neighbor pass, matching SNAP/k2b/MTP's own compute_descriptor_<family>_global
// design exactly (see compute_descriptor_k2b_global.cu for the reference shape this mirrors).
// Purely additive: the existing per-atom compute_descriptor_pod stays fully intact and usable alone.
//
// Multi-species (nelements>1) support: PodGlobalOp's old per-pair scatter placed every gradient
// contribution's column using the CENTRAL atom's type (ti0) for BOTH the central (+=) and neighbor
// (-=) side of a pair. pod_descriptor_op.h's pda_* aggregate now preserves that same information by
// being WIDENED by nelements (one slot per possible central-atom species, see its own header
// comment) instead of collapsing every central-role type into one bin -- so a given atom's own
// aggregate row spans multiple ti0 slots (its own species for its central/self terms, plus one slot
// per OTHER species it was ever a neighbor of). Row 0 (one-body + descriptor sum) only ever needs
// the atom's OWN type (available directly at assembly time via field::type). The gradient rows,
// read from pda_*, loop over every ti0 in [0,nelements) to recover all of an atom's contributions --
// mono-species (nelements==1) is the trivial single-iteration special case of the same loop, not a
// separately-maintained path.
//
// Row/column layout: row 0 = system-wide per-element summed descriptor vector (incl. the nl1
// one-body atom-count term); rows 1..3*natoms = ALREADY FORCE-SIGNED gradient of row 0 w.r.t. atom
// id's x/y/z (F_atom = +coeff . row, not -coeff . row -- see the finite-difference-vs-energy note
// below, this is not what the central+=/neighbor-= scatter's naming naively suggests); rows
// 3*natoms+1..+6 = virial, Voigt order [xx,yy,zz,yz,xz,xy], same already-force-signed convention.
// Rows 0..3*natoms match LAMMPS compute pod/global exactly (single-rank); virial rows have no
// LAMMPS counterpart (POD fitting there doesn't use stress) but follow the same formula
// compute_descriptor_snap_global uses.
//
// Force/virial sign, verified independently (2026-09-22): compute_descriptor_pod's own
// documentation (and this operator's own gradient/virial rows, built from the exact same pda_*
// aggregate) previously claimed "F_atom = -coeff . row" -- a finite-difference-vs-energy check
// (perturb a small non-periodic cluster's positions/strain, central-difference row 0 against the
// gradient/virial rows directly, independent of any LAMMPS comparison) showed this is backwards:
// the central+=/neighbor-= scatter in pod_descriptor_op.h computes d(rij)/dr with rij=r_neighbor-
// r_central, so by the chain rule the stored aggregate is -dE/dr (already force-signed), not
// +dE/dr. LAMMPS's own compute pod/global apparently uses the identical convention (that's why the
// gradient-row VALUES still matched LAMMPS exactly in compare_global.py -- a plain values-vs-LAMMPS
// diff can't catch an overall sign both sides happen to share). Correct usage: F_atom = +coeff .
// row[1+3*id+xyz] directly, no extra negation. See
// data/regression_new/compute_descriptor/test_pod_descriptors/compare_global_strain_fd.py.
//
// Must run AFTER compute_descriptor_pod: { compute_derivative: true } and BEFORE any
// update_opt_from_ghost call on the pda_* fields -- this operator needs the raw, per-rank-local,
// UN-FOLDED aggregate (real and ghost slots each carry their own partial view); update_opt_from_ghost
// folding first would corrupt both the gradient-row and virial-row accumulation here. Same ordering
// rule as compute_descriptor_snap_global.cu/compute_descriptor_k2b_global.cu.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDescriptorPodGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi  , INPUT , REQUIRED );
    ADD_SLOT( GridT    , grid , INPUT , REQUIRED );
    ADD_SLOT( PodContext , pod_ctx , INPUT , REQUIRED , DocString{"still needed for Mdesc/nClusters/nCoeffPerElement/nl1/nelements"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , pod_descriptors , INPUT , OPTIONAL , DocString{"see compute_descriptor_pod; required"} );
    ADD_SLOT( long     , ncoeff , INPUT , OPTIONAL , DocString{"see compute_descriptor_pod; required (= Mdesc*nClusters)"} );
    ADD_SLOT( std::string , deriv_agg_field_prefix , INPUT , std::string("pda_")
            , DocString{"Must match compute_descriptor_pod's own deriv_agg_field_prefix. Must be read before update_opt_from_ghost runs on these fields -- see this file's header comment."} );

    ADD_SLOT( onika::memory::CudaMMVector<double> , pod_global , OUTPUT
            , DocString{"Row-major (1+3*natoms+6) x ncoeff_all global array, natoms = total atom count across every MPI rank. See this file's header comment for the exact row/column layout."} );
    ADD_SLOT( long , ncoeff_all , OUTPUT , DocString{"Number of columns = nCoeffPerElement*nelements"} );

  public:
    inline void execute() override final
    {
      if( ! pod_descriptors.has_value() || ! ncoeff.has_value() )
      {
        fatal_error() << "compute_descriptor_pod_global: pod_descriptors/ncoeff unavailable -- run compute_descriptor_pod first" << std::endl;
      }

      auto & eapod0 = *pod_ctx->m_eapod[0];

      const long Mdesc = eapod0.Mdesc;
      const long nClusters = eapod0.nClusters;
      const long nl1 = eapod0.nl1;
      const long nelements = eapod0.nelements;
      const long ncols = static_cast<long>(eapod0.nCoeffPerElement) * nelements;
      *ncoeff_all = ncols;
      const long nc = *ncoeff;   // == Mdesc*nClusters
      const long nc3 = nc * 3 * nelements;   // widened by nelements, matches compute_descriptor_pod.cu

      std::vector<const double*> agg_ptr( static_cast<size_t>(nc3), nullptr );
      for( long k=0; k<nc3; k++ )
      {
        agg_ptr[k] = grid->flat_array_data_nocreate( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
        if( agg_ptr[k] == nullptr )
        {
          fatal_error() << "compute_descriptor_pod_global: field '"<<*deriv_agg_field_prefix<<k<<"' not found -- run compute_descriptor_pod with compute_derivative: true first" << std::endl;
        }
      }

      const auto * cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells = grid->number_of_cells();

      // local owned (non-ghost) atom count -> global total via one small Allreduce, so every rank
      // allocates the same full-size array before the main per-atom Allreduce below
      long local_owned = 0;
      for( size_t ci=0; ci<n_cells; ci++ )
      {
        if( grid->is_ghost_cell(ci) ) continue;
        local_owned += static_cast<long>( grid->cell(ci).size() );
      }
      long natoms_global = 0;
      MPI_Allreduce( &local_owned, &natoms_global, 1, MPI_LONG, MPI_SUM, *mpi );

      const long grad_rows = 1 + 3*natoms_global;
      const long virial_row0 = grad_rows;
      const long rows = grad_rows + 6;

      pod_global->clear();
      pod_global->resize( static_cast<size_t>(rows) * static_cast<size_t>(ncols), 0.0 );
      double * const arr = pod_global->data();

      for( size_t ci=0; ci<n_cells; ci++ )
      {
        const bool is_ghost = grid->is_ghost_cell(ci);
        const auto & cell = grid->cell(ci);
        const size_t np = cell.size();
        for( size_t pi=0; pi<np; pi++ )
        {
          const size_t p = cell_particle_offset[ci] + pi;
          const uint64_t id = cell[field::id][pi];

          if( ! is_ghost )
          {
            // Row 0 (one-body + descriptor sum) only ever needs THIS atom's own species -- no
            // per-contribution type ambiguity here, unlike the gradient rows below (pod_descriptors
            // itself was never widened, see compute_descriptor_pod.cu's own comment on why not).
            const long ti0 = pod_ctx->type_map[ cell[field::type][pi] ] - 1;
            if( nl1 > 0 ) arr[ eapod0.nCoeffPerElement*ti0 ] += 1.0; // one-body atom-count term, row 0
            const double * const src = pod_descriptors->data() + static_cast<size_t>(nc) * p;
            for( long m=0; m<Mdesc; m++ )
            for( long j=0; j<nClusters; j++ )
            {
              const long col = eapod0.nCoeffPerElement*ti0 + nl1 + m + j*Mdesc;
              arr[ col ] += src[ m + Mdesc*j ];
            }
          }

          // Gradient rows: this atom's pda_* aggregate spans one slot per possible CENTRAL-atom
          // species it was ever involved with (its own, for central/self terms; every OTHER species
          // it was ever a neighbor of, for neighbor terms) -- loop every ti0 to recover all of it.
          // Mono-species (nelements==1) is the trivial single-iteration case of this same loop.
          const long grad_row0 = 1 + 3*static_cast<long>(id);
          for( long ti0=0; ti0<nelements; ti0++ )
          for( long m=0; m<Mdesc; m++ )
          for( long j=0; j<nClusters; j++ )
          {
            const long col = eapod0.nCoeffPerElement*ti0 + nl1 + m + j*Mdesc;
            const long comp = (m + Mdesc*j + Mdesc*nClusters*ti0)*3;
            const double dx = agg_ptr[comp+0][p];
            const double dy = agg_ptr[comp+1][p];
            const double dz = agg_ptr[comp+2][p];
            arr[ (grad_row0+0)*ncols + col ] += dx;
            arr[ (grad_row0+1)*ncols + col ] += dy;
            arr[ (grad_row0+2)*ncols + col ] += dz;
          }
        }
      }

      // Combine every rank's local partial contribution to rows 0..3*natoms_global.
      MPI_Allreduce( MPI_IN_PLACE, arr, static_cast<int>(grad_rows*ncols), MPI_DOUBLE, MPI_SUM, *mpi );

      // Virial rows: Σ r_atom . dDescriptor_atom/dr_atom, Voigt order [xx,yy,zz,yz,xz,xy] --
      // computed from the now-COMPLETE (post-Allreduce) gradient rows above, each rank summing only
      // over its own owned atoms (every atom is owned by exactly one rank, so no double counting and
      // nothing missed); one more (small) Allreduce combines every rank's partial virial sum.
      for( size_t ci=0; ci<n_cells; ci++ )
      {
        if( grid->is_ghost_cell(ci) ) continue;
        const auto & cell = grid->cell(ci);
        const size_t np = cell.size();
        for( size_t pi=0; pi<np; pi++ )
        {
          const uint64_t id = cell[field::id][pi];
          const double rx = cell[field::rx][pi];
          const double ry = cell[field::ry][pi];
          const double rz = cell[field::rz][pi];
          const long grad_row0 = 1 + 3*static_cast<long>(id);
          for( long k=0; k<ncols; k++ )
          {
            const double dx = arr[ (grad_row0+0)*ncols + k ];
            const double dy = arr[ (grad_row0+1)*ncols + k ];
            const double dz = arr[ (grad_row0+2)*ncols + k ];
            arr[ (virial_row0+0)*ncols + k ] += dx*rx; // xx
            arr[ (virial_row0+1)*ncols + k ] += dy*ry; // yy
            arr[ (virial_row0+2)*ncols + k ] += dz*rz; // zz
            arr[ (virial_row0+3)*ncols + k ] += dz*ry; // yz
            arr[ (virial_row0+4)*ncols + k ] += dz*rx; // xz
            arr[ (virial_row0+5)*ncols + k ] += dy*rx; // xy
          }
        }
      }
      MPI_Allreduce( MPI_IN_PLACE, arr + static_cast<size_t>(virial_row0)*ncols, static_cast<int>(6*ncols), MPI_DOUBLE, MPI_SUM, *mpi );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Global-array analogue of LAMMPS's compute pod/global, extended with 6 virial rows LAMMPS's own
pod/global doesn't have (POD fitting there just doesn't use stress, not a structural limitation):
a single, row-major (1+3*natoms+6) x ncoeff_all array (ncoeff_all = nCoeffPerElement*nelements).
Multi-species (nelements>1) supported -- see this file's header comment for the pda_* widening
this relies on.

  row 0             -- system-wide per-element summed descriptor vector, including the nl1
                        one-body atom-count term. Dot this with a coefficient vector to get the
                        total configuration energy.
  rows 1..3*natoms  -- ALREADY FORCE-SIGNED gradient of row 0 w.r.t. atom id's x/y/z (verified by
                        finite difference against row 0, see this file's header comment). Dot row
                        (1+3*id+xyz) with the same coefficient vector directly to get that atom's
                        force component: F = +coeff . row (no extra negation).
  rows 3N+1..3N+6   -- sum over atoms of position . gradient row, Voigt order
                        [xx,yy,zz,yz,xz,xy], same already-force-signed convention. Dot with the
                        same coefficient vector to get the virial/stress tensor component.

Purely additive: reads compute_descriptor_pod's existing output (pod_descriptors + compute_derivative:
true's pda_* aggregate) rather than re-running the descriptor+derivative pass, so the existing
per-atom descriptor capability stays fully intact and usable on its own. Must run right after
compute_descriptor_pod (compute_derivative: true) and BEFORE any update_opt_from_ghost call on its
aggregate fields -- this operator needs the raw, un-folded per-rank-local aggregate.

Usage example:

init_parameters:
  - species
  - pod_init: { parameters: { pod_file: "Ta_param.pod", coeff_file: "Ta_coefficients.pod" } }

compute_descriptor_pod: { compute_derivative: true }
compute_descriptor_pod_global
write_descriptor_pod_global: { filename: "pod_global.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_descriptor_pod_global)
  {
    OperatorNodeFactory::instance()->register_factory("compute_descriptor_pod_global", make_grid_variant_operator<ComputeDescriptorPodGlobal>);
  }

}
