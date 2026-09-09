/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/
#include <mpi.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/reduce_cell_particles.h>
#include <memory>

namespace exaStamp {

  ONIKA_HOST_DEVICE_FUNC inline void ATOMIC_ADD(Mat3d& a, const Mat3d& b) {
    ONIKA_CU_ATOMIC_ADD(a.m11, b.m11); ONIKA_CU_ATOMIC_ADD(a.m12, b.m12); ONIKA_CU_ATOMIC_ADD(a.m13, b.m13);
    ONIKA_CU_ATOMIC_ADD(a.m21, b.m21); ONIKA_CU_ATOMIC_ADD(a.m22, b.m22); ONIKA_CU_ATOMIC_ADD(a.m23, b.m23);
    ONIKA_CU_ATOMIC_ADD(a.m31, b.m31); ONIKA_CU_ATOMIC_ADD(a.m32, b.m32); ONIKA_CU_ATOMIC_ADD(a.m33, b.m33);
  }

  struct FdotrValue {
    Mat3d vir_tot;
  };

  struct ReduceFdotrFunctor {

    ONIKA_HOST_DEVICE_FUNC inline void operator()(FdotrValue& local, const double fx, const double fy, const double fz,
                                                   const double rx, const double ry, const double rz, reduce_thread_local_t = {}) const {
      local.vir_tot += tensor(Vec3d{fx, fy, fz}, Vec3d{rx, ry, rz});
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator()(FdotrValue& global, const FdotrValue& local, reduce_thread_block_t) const {
      ATOMIC_ADD(global.vir_tot, local.vir_tot);
    }

    ONIKA_HOST_DEVICE_FUNC inline void operator()(FdotrValue& global, const FdotrValue& local, reduce_global_t) const {
      ATOMIC_ADD(global.vir_tot, local.vir_tot);
    }
  };

}

namespace exanb
{

  template <>
  struct ReduceCellParticlesTraits< exaStamp::ReduceFdotrFunctor >
  {
    static inline constexpr bool CudaCompatible = true;
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool RequiresCellParticleIndex = false;
  };

}

// Global virial/stress tensor via LAMMPS's cheap "fdotr" trick: Sum_i (F_i (x) r_i), the outer
// product of each particle's TOTAL force with its ABSOLUTE position, summed over every particle --
// no per-pair bookkeeping, no per-atom field::_virial needed. Mirrors LAMMPS's
// Pair::virial_fdotr_compute() (src/pair.cpp) exactly: a single O(N) pass over owned+ghost atoms
// using only fields (_fx/_fy/_fz, _rx/_ry/_rz) that already exist for any MD run, regardless of
// potential -- same idea as how this codebase's own momentum reduction (mass*velocity, in
// simulation_thermodynamic_state.cpp) needs no potential-specific per-atom field either.
//
// CORRECTNESS REQUIREMENT (not a tunable, not a performance knob): this sum MUST include ghost
// atoms with their true periodic-image coordinates (ghost slot defaults to true, unlike every
// other operator of this shape in the codebase), and MUST run BEFORE any ghost-force fold-back
// (update_force_energy_from_ghost) -- a ghost's force contribution has to still be "in place" on
// the ghost's own (correctly wrapped) position for the telescoping Newton's-third-law cancellation
// across periodic boundaries to work out, exactly as documented for LAMMPS's own virial_fdotr_compute
// (see LAMMPS doc/src/Developer_flow.rst and Developer_write_pair.rst). Typical placement:
//
//   compute_force_epilog:
//     - compute_sum_fdotr
//     - update_force_energy_from_ghost
//     - force_to_accel
//
// Uses exanb::reduce_cell_particles (same GPU/CUDA-capable reduction machinery as
// src/compute/sum_forces.cu), with enable_ghosts=true. NOTE: this exercised a real bug in
// exaNBody's reduce_cell_particles.h -- when enable_ghosts=true (m_ghost_layers==0), the flat
// cell index was only ever computed inside an `if (m_ghost_layers != 0)` guard, so it was NEVER
// assigned in the ghost-inclusive case, leaving it at its uninitialized sentinel value (silently
// reading zero-valued/garbage cells instead of asserting or crashing in a release build). Fixed
// upstream in exaNBody/src/compute/include/exanb/compute/reduce_cell_particles.h to always compute
// the flat index, matching compute_cell_particles.h's own (correct, unconditional) sibling code --
// no other caller of reduce_cell_particles in this codebase passed enable_ghosts=true before, so
// this path was previously untested.
//
// Standalone operator only for now -- NOT wired into simulation_thermodynamic_state.cpp.
namespace exaStamp {

  using namespace exanb;

  template <
    typename GridT,
    class = AssertGridHasFields<GridT, field::_fx, field::_fy, field::_fz>
    >
  class ComputeSumFdotr : public OperatorNode
  {
    using ReduceFields = FieldSet<field::_fx, field::_fy, field::_fz, field::_rx, field::_ry, field::_rz>;
    static constexpr ReduceFields reduce_field_set{};

    ADD_SLOT(MPI_Comm, mpi, INPUT, MPI_COMM_WORLD);
    ADD_SLOT(GridT, grid, INPUT, REQUIRED);
    ADD_SLOT(bool, ghost, INPUT, true,
             DocString{"Include ghost atoms in the sum. Must stay true for correctness under periodic "
                       "boundary conditions -- this is the same requirement LAMMPS's virial_fdotr_compute() "
                       "has (owned+ghost, before ghost-force fold-back), not a performance knob."});
    ADD_SLOT(Mat3d, out, OUTPUT, DocString{"Sum_i (F_i (x) r_i) -- the global virial/stress tensor, "
                                            "computed via the O(N) LAMMPS-style 'fdotr' trick instead of "
                                            "a per-pair virial tally or a per-atom field::_virial field."});

    inline std::string documentation() const final {
      return R"EOF(
        Computes the global virial/stress tensor as Sum_i (F_i (x) r_i) over every particle
        (owned + ghost), mirroring LAMMPS's Pair::virial_fdotr_compute() -- an O(N) reduction over
        fields (force, position) that already exist for any MD run, needing no per-pair virial
        bookkeeping in any force kernel and no per-atom field::_virial field.

        Must run in compute_force_epilog BEFORE update_force_energy_from_ghost (needs the raw,
        un-folded ghost force contributions still on their own wrapped positions -- see this file's
        header comment for the full correctness argument). Standalone operator for now, not wired
        into simulation_thermodynamic_state.

        YAML example:

          compute_force_epilog:
            - compute_sum_fdotr
            - update_force_energy_from_ghost
            - force_to_accel
      )EOF";
    }

    // A pure-OUTPUT slot with no downstream consumer gets silently pruned from the execution graph
    // (this operator would simply never run) unless explicitly marked as a sink -- same reason
    // print_thermodynamic_state.cpp/grid_clear.cpp do this. This operator's whole purpose (for now)
    // is its printed value, not feeding another operator, so it must always be a sink.
    inline bool is_sink() const override final { return true; }

  public:
    inline void execute() final {

      FdotrValue value = {};
      ReduceFdotrFunctor func;
      reduce_cell_particles(*grid, *ghost, func, value, reduce_field_set, parallel_execution_context());

      double local[9] = { value.vir_tot.m11, value.vir_tot.m12, value.vir_tot.m13,
                           value.vir_tot.m21, value.vir_tot.m22, value.vir_tot.m23,
                           value.vir_tot.m31, value.vir_tot.m32, value.vir_tot.m33 };
      double global[9] = {0.,0.,0., 0.,0.,0., 0.,0.,0.};
      MPI_Allreduce(&local, &global, 9, MPI_DOUBLE, MPI_SUM, *mpi);
      *out = Mat3d{ global[0], global[1], global[2], global[3], global[4], global[5], global[6], global[7], global[8] };

      lout << "compute_sum_fdotr: Sum_i(F_i (x) r_i) = " << *out << std::endl;
    }
  };

  template<class GridT> using ComputeSumFdotrTmpl = ComputeSumFdotr<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_sum_fdotr) {
    OperatorNodeFactory::instance()->register_factory("compute_sum_fdotr", make_grid_variant_operator<ComputeSumFdotrTmpl>);
  }

}
