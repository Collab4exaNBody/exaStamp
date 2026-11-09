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
#include <exanb/core/domain.h>

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

// Global virial tensor Sum_i (F_i (x) r_i), summed over owned AND ghost particles in a single O(N)
// pass over the force and position fields: no per-pair virial tally, no per-atom virial field.
// It must include ghosts (with their periodic-image positions) and run BEFORE the ghost force
// fold-back (update_force_energy_from_ghost): the ghost force contributions must still sit on the
// ghost positions for the Newton's-third-law cancellation across periodic boundaries to hold.
// Positions are stored in grid space; the real-frame result is (Sum_i F_i (x) r_i) . xform^T.
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
    ADD_SLOT(Domain, domain, INPUT, REQUIRED);
    ADD_SLOT(bool, ghost, INPUT, true,
             DocString{"Include ghost atoms in the sum. Must stay true for correctness under periodic "
                       "boundary conditions (owned + ghost, before the ghost force fold-back)."});
    ADD_SLOT(Mat3d, out, OUTPUT, DocString{"Sum_i (F_i (x) r_i): the global virial tensor, computed in a single "
                                            "O(N) pass over forces and positions."});

    inline std::string documentation() const final {
      return R"EOF(
        Computes the global virial tensor Sum_i (F_i (x) r_i) over every particle (owned + ghost),
        in a single O(N) reduction over the force and position fields. It needs no per-pair virial
        bookkeeping in the force kernels and no per-atom virial field, so it works with any potential.

        Must run in compute_force_epilog BEFORE update_force_energy_from_ghost: the ghost force
        contributions must still sit on the ghost positions. The result is printed at each call.

        YAML example:

          compute_force_epilog:
            - compute_sum_fdotr
            - update_force_energy_from_ghost
            - force_to_accel
      )EOF";
    }

    // pure-output operator whose result is only printed: mark it as a sink so it is not pruned
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
      // positions are in grid space: F (x) (X r) = (F (x) r) X^T
      *out = Mat3d{ global[0], global[1], global[2], global[3], global[4], global[5], global[6], global[7], global[8] } * transpose( domain->xform() );

      lout << "compute_sum_fdotr: Sum_i(F_i (x) r_i) = " << *out << std::endl;
    }
  };

  template<class GridT> using ComputeSumFdotrTmpl = ComputeSumFdotr<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_sum_fdotr) {
    OperatorNodeFactory::instance()->register_factory("compute_sum_fdotr", make_grid_variant_operator<ComputeSumFdotrTmpl>);
  }

}
