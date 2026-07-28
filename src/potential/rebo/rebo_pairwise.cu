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

#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>

#include <exanb/core/make_grid_variant_operator.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <exanb/core/particle_type_id.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <memory>
#include <mpi.h>
#include <vector>

#include "rebo_force_op_cached.h"
#include "rebo_params.h"

namespace exaStamp
{
  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  // Pure two-body pass: Brenner repulsive term VR only (REBOPairwiseForceOp,
  // see rebo_force_op_cached.h). No coordination-number fields needed. Split
  // out from the many-body pass (rebo_manybody.cu) so the two can be run and
  // profiled as separate kernels.
  template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz>> class ReboPairwiseForce : public OperatorNode
  {
    ADD_SLOT(MPI_Comm, mpi, INPUT, REQUIRED);
    ADD_SLOT(double, rcut_max, INPUT_OUTPUT, 0.0);
    // REBOPairwiseForceOp (VR only) is exactly zero beyond rcmax — no need for
    // the wider rebo_cutoff (3*rcmax) that ManyBody's lm-subloop requires.
    // bondorder_cutoff (= rcmax_max, see rebo_init.cu) is the correctly-sized
    // filter here, same as rebo_bond_order/rebo_conjugation already use. The
    // underlying neighbor list is still built out to rcut_max regardless —
    // this only trims how many of those candidates get decoded/transformed/
    // buffered per particle for this specific pass.
    ADD_SLOT(double, bondorder_cutoff, INPUT, REQUIRED);
    ADD_SLOT(exanb::GridChunkNeighbors, chunk_neighbors, INPUT, exanb::GridChunkNeighbors{}, DocString{"neighbor list"});
    ADD_SLOT(bool, ghost, INPUT, false);
    ADD_SLOT(bool, conv_coef_units, INPUT, true);
    ADD_SLOT(bool, trigger_thermo_state, INPUT, OPTIONAL);
    ADD_SLOT(GridT, grid, INPUT_OUTPUT);
    ADD_SLOT(Domain, domain, INPUT, REQUIRED);
    ADD_SLOT(GridParticleLocks, particle_locks, INPUT, OPTIONAL, DocString{"particle spin locks"});
    ADD_SLOT(long, timestep, INPUT, REQUIRED, DocString{"Iteration number"});
    ADD_SLOT(ParticleSpecies, species, INPUT, REQUIRED);
    ADD_SLOT(ReboParams, parameters, INPUT, REQUIRED);

    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBufferRebo = ComputePairBuffer2<UseWeights, UseNeighbors, ReboComputeBuffer, CopyParticleType, REBO_MAX_NEIGHBORS>;

    using CellParticles = typename GridT::CellParticles;
    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;

    using ComputeFieldsWithoutVirial = FieldSet<field::_ep, field::_fx, field::_fy, field::_fz, field::_type>;
    using ComputeFieldsWithVirial = FieldSet<field::_ep, field::_fx, field::_fy, field::_fz, field::_type, field::_virial>;
    using ComputeFields = std::conditional_t<has_virial_field, ComputeFieldsWithVirial, ComputeFieldsWithoutVirial>;
    static constexpr ComputeFields compute_force_field_set{};

  public:
    inline void execute() override final
    {
      assert(chunk_neighbors->number_of_cells() == grid->number_of_cells());
      if (grid->number_of_cells() == 0)
        return;

      // bool eflag = false;
      if (trigger_thermo_state.has_value())
      {
        ldbg << "trigger_thermo_state = " << *trigger_thermo_state << std::endl;
        // eflag = *trigger_thermo_state;
      }

      // Read Only Rebo parameters. Allocated in CUDA managed memory (not on the
      // host stack) so the pointer stashed in the force op below stays valid
      // when dereferenced from GPU kernel code.
      onika::memory::CudaMMVector<ReboParamsRO> params_ro_mm(1, ReboParamsRO{*parameters});
      const ReboParamsRO *params_ro = params_ro_mm.data();

      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{*chunk_neighbors};
      LinearXForm cp_xform{domain->xform()};

      ComputePairNullWeightIterator cp_weight{};
      auto compute_buf_rebo = make_compute_pair_buffer<ComputeBufferRebo>();

      const bool use_locks = omp_get_max_threads() > 1 && particle_locks.has_value();
      auto pairwise_op_fields = make_field_tuple_from_field_set(ComputeFields{});

      auto with_locks = [&](auto cp_locks) { compute_cell_particle_pairs(*grid, *bondorder_cutoff, *ghost, make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks), compute_buf_rebo, REBOPairwiseForceOp{params_ro}, pairwise_op_fields, parallel_execution_context()); };
      if (use_locks)
        with_locks(ComputePairOptionalLocks<true>{particle_locks->data()});
      else
        with_locks(ComputePairOptionalLocks<false>{});
    }
  };

  template <class GridT> using ReboPairwiseForceTmpl = ReboPairwiseForce<GridT>;

  ONIKA_AUTORUN_INIT(rebo_pairwise) { OperatorNodeFactory::instance()->register_factory("rebo_pairwise", make_grid_variant_operator<ReboPairwiseForceTmpl>); }

} // namespace exaStamp
