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
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <cstdint>
#include <memory>
#include <vector>
#include <mpi.h>

#include "pod_params.h"
#include "pod_config.h"
#include "pod_force_op.h"    // PodComputeBuffer, CopyParticleType
#include "pod_global_op.h"   // PodGlobalOp

// Global-array descriptor+gradient+virial pass, analogous to LAMMPS's compute pod/global (plus 6
// virial rows LAMMPS's own pod/global doesn't have -- POD fitting there just doesn't use stress,
// not a structural limitation; added here to match compute_descriptor_snap_global's shape) -- see
// pod_global_op.h for the descriptor+gradient algorithm and why rows are indexed by atom id-1.
// Single MPI rank only, matching compute_pod_global.cpp's own restriction (no MPI reduction exists
// there either) -- this also means the virial pass below needs no Allreduce/pre-fold timing care,
// unlike compute_descriptor_snap_global.cu's own virial pass.
namespace exaStamp
{

  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_type >
    >
  class ComputeDescriptorPodGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                  , mpi             , INPUT , REQUIRED );
    ADD_SLOT( double                    , rcut_max        , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{}, DocString{"neighbor list"} );
    ADD_SLOT( bool                      , ghost           , INPUT , false );
    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain          , INPUT , REQUIRED );
    ADD_SLOT( PodContext                , pod_ctx         , INPUT , REQUIRED );

    ADD_SLOT( onika::memory::CudaMMVector<double>, pod_global, OUTPUT,
               DocString{"Row-major (1+3*natoms+6) x ncoeff_all global array, natoms = number of owned (non-ghost) particles: row 0 = system-wide per-element summed descriptor vector (incl. one-body atom-count term); rows 1..3*natoms = gradient of row 0 w.r.t. atom (id-1)'s x/y/z; rows 3*natoms+1..+6 = virial, Voigt order [xx,yy,zz,yz,xz,xy]. Rows 0..3*natoms match LAMMPS compute pod/global exactly; the virial rows have no LAMMPS counterpart (POD fitting there doesn't use stress) but follow the same formula compute_descriptor_snap_global uses. Single MPI rank only."} );
    ADD_SLOT( long, ncoeff_all, OUTPUT,
               DocString{"Number of columns = nCoeffPerElement*nelements"} );

    static constexpr bool UseWeights   = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights, UseNeighbors, PodComputeBuffer, CopyParticleType>;
    static constexpr FieldSet<field::_type> compute_global_field_set{};

  public:

    inline void execute() override final
    {
      int nprocs = 1;
      MPI_Comm_size(*mpi, &nprocs);
      if (nprocs > 1)
      {
        fatal_error() << "compute_descriptor_pod_global: only supported on a single MPI rank"
                      << " (matches LAMMPS compute pod/global's own restriction -- it has no MPI reduction either)" << std::endl;
      }

      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      const size_t nt = omp_get_max_threads();
      if (nt > pod_ctx->m_eapod.size())
      {
        lerr << "POD: omp_get_max_threads() grew from " << pod_ctx->m_eapod.size()
             << " to " << nt << " after init -- some threads lack an EAPOD context."
             << " Re-run with the correct OMP_NUM_THREADS set before launch." << std::endl;
        fatal_error() << "POD thread context size mismatch" << std::endl;
      }

      if (grid->number_of_cells() == 0) { *ncoeff_all = 0; return; }

      auto& eapod0 = *pod_ctx->m_eapod[0];
      const long ncols = static_cast<long>(eapod0.nCoeffPerElement) * eapod0.nelements;
      *ncoeff_all = ncols;

      const size_t natoms = grid->number_of_particles() - grid->number_of_ghost_particles();
      const size_t rows = 1 + 3*natoms + 6;
      pod_global->clear();
      pod_global->resize( rows * ncols, 0.0 );

      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto global_buf = make_compute_pair_buffer<ComputeBuffer>();
      LinearXForm cp_xform{ domain->xform() };
      ComputePairOptionalLocks<false> cp_locks{};

      PodGlobalOp global_op{ pod_ctx->m_eapod, pod_ctx->type_map, ncols, pod_global->data() };
      compute_cell_particle_pairs(
          *grid, *rcut_max, *ghost,
          make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks),
          global_buf, global_op, compute_global_field_set,
          parallel_execution_context());

      // Virial rows: Σ r_atom . dDescriptor_atom/dr_atom over owned atoms, Voigt order
      // [xx,yy,zz,yz,xz,xy] -- same formula as compute_descriptor_snap_global.cu's virial pass.
      // Unlike SNAP, no local/pre-fold timing subtlety here: this operator is single-MPI-rank
      // only (see above), and the (1+3N)-row block above is already the complete, final answer
      // (ghost contributions already collapsed onto their owner's id-indexed row by the fused
      // pair-loop's own atomic scatter) -- so this is just a plain finalization pass reading rows
      // already written above, no separate reduction step needed.
      {
        double * const arr = pod_global->data();
        const long virial_row0 = 1 + 3*static_cast<long>(natoms);
        const size_t n_cells = grid->number_of_cells();
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
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Global-array analogue of LAMMPS's compute pod/global, extended with 6 virial rows LAMMPS's own
pod/global doesn't have (POD fitting there just doesn't use stress, not a structural limitation):
a single, row-major (1+3*natoms+6) x ncoeff_all array (ncoeff_all = nCoeffPerElement*nelements):

  row 0             -- system-wide per-element summed descriptor vector, including the nl1
                        one-body atom-count term. Dot this with a coefficient vector to get the
                        total configuration energy.
  rows 1..3*natoms  -- the gradient of row 0 w.r.t. atom id's x/y/z. Dot row (1+3*id+xyz)
                        with the same coefficient vector and negate to get that atom's force
                        component: F = -coeff . row.
  rows 3N+1..3N+6   -- sum over atoms of position . gradient row, Voigt order
                        [xx,yy,zz,yz,xz,xy]. Dot with the same coefficient vector to get the
                        virial/stress tensor component.

This is the design-matrix structure needed to fit a linear POD potential against total energy plus
per-atom forces plus virial (stack these rows across many training configurations, stack the
corresponding energy/force/virial targets, solve by least squares).

Single MPI rank only -- matches compute_pod_global.cpp's own restriction (it has no MPI reduction
either). Rows are indexed by atom id rather than internal particle order: this makes row order
directly comparable to LAMMPS's own output (for rows 0..3*natoms), and makes ghost/real folding
automatic (a ghost carries the same field::id as its real counterpart), so no update_opt_from_ghost
step is needed here.

Usage example:

pod_init: { parameters: { pod_file: "Ta_param.pod", coeff_file: "Ta_coefficients.pod" } }
compute_descriptor_pod_global
write_descriptor_pod_global: { filename: "pod_global.txt" }

)EOF";
    }
  };

  template<class GridT> using ComputeDescriptorPodGlobalTmpl = ComputeDescriptorPodGlobal<GridT>;

  ONIKA_AUTORUN_INIT(compute_descriptor_pod_global)
  {
    OperatorNodeFactory::instance()->register_factory("compute_descriptor_pod_global", make_grid_variant_operator<ComputeDescriptorPodGlobalTmpl>);
  }

}
