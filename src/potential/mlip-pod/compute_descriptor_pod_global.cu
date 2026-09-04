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

#include <memory>
#include <vector>
#include <mpi.h>

#include "pod_params.h"
#include "pod_config.h"
#include "pod_force_op.h"    // PodComputeBuffer, CopyParticleType
#include "pod_global_op.h"   // PodGlobalOp

// Global-array descriptor+gradient pass, analogous to LAMMPS's compute pod/global -- see
// pod_global_op.h for the exact algorithm and why rows are indexed by atom id-1. Single MPI rank
// only, matching compute_pod_global.cpp's own restriction (no MPI reduction exists there either).
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
               DocString{"Row-major (1+3*natoms) x ncoeff_all global array, natoms = number of owned (non-ghost) particles: row 0 = system-wide per-element summed descriptor vector (incl. one-body atom-count term); rows 1..3*natoms = gradient of row 0 w.r.t. atom (id-1)'s x/y/z. Matches LAMMPS compute pod/global exactly. Single MPI rank only."} );
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
      const size_t rows = 1 + 3*natoms;
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
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Global-array analogue of LAMMPS's compute pod/global: a single, row-major
(1+3*natoms) x ncoeff_all array (ncoeff_all = nCoeffPerElement*nelements):

  row 0            -- system-wide per-element summed descriptor vector, including the nl1
                       one-body atom-count term. Dot this with a coefficient vector to get the
                       total configuration energy.
  rows 1..3*natoms -- the gradient of row 0 w.r.t. atom (id-1)'s x/y/z. Dot row (1+3*(id-1)+xyz)
                       with the same coefficient vector and negate to get that atom's force
                       component: F = -coeff . row.

This is the design-matrix structure needed to fit a linear POD potential against total energy plus
per-atom forces (stack these rows across many training configurations, stack the corresponding
energy/force targets, solve by least squares).

Single MPI rank only -- matches compute_pod_global.cpp's own restriction (it has no MPI reduction
either). Rows are indexed by atom id-1 rather than internal particle order: this makes row order
directly comparable to LAMMPS's own output, and makes ghost/real folding automatic (a ghost carries
the same field::id as its real counterpart), so no update_opt_from_ghost step is needed here.

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
