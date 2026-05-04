/*
exaDEMLicensed to the Apache Software Foundation (ASF) under one
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
#include <exanb/core/particle_type_properties.h>

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>

namespace exaStamp {
  
  ONIKA_HOST_DEVICE_FUNC inline void ATOMIC_ADD(Vec3d& a, const Vec3d& b) {
    ONIKA_CU_ATOMIC_ADD(a.x, b.x);
    ONIKA_CU_ATOMIC_ADD(a.y, b.y);
    ONIKA_CU_ATOMIC_ADD(a.z, b.z);
  }

  struct ForceValue {
    Vec3d f_tot;
  };
  
  struct ReduceForceFunctor {

    ONIKA_HOST_DEVICE_FUNC inline void operator()(ForceValue& local, const double fx, const double fy, const double fz, reduce_thread_local_t = {}) const {
      Vec3d f = {fx, fy, fz};
      local.f_tot += f;
    }
    
    ONIKA_HOST_DEVICE_FUNC inline void operator()(ForceValue& global, const ForceValue& local, reduce_thread_block_t) const {
      ATOMIC_ADD(global.f_tot, local.f_tot);
    }
    
    ONIKA_HOST_DEVICE_FUNC inline void operator()(ForceValue& global, const ForceValue& local, reduce_global_t) const {
      ATOMIC_ADD(global.f_tot, local.f_tot);
    }
  };
  
}

namespace exanb
{
  
  template <>
  struct ReduceCellParticlesTraits< exaStamp::ReduceForceFunctor >
  {
    static inline constexpr bool CudaCompatible = true;
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool RequiresCellParticleIndex = false;
  };
  
}

namespace exaStamp {
  
  using namespace exanb;

  template <
    typename GridT,
    class = AssertGridHasFields<GridT, field::_fx, field::_fy, field::_fz>
    >
  class SumForces : public OperatorNode
  {
    using ReduceFields = FieldSet<field::_fx, field::_fy, field::_fz>;
    static constexpr ReduceFields reduce_field_set{};
    
    ADD_SLOT(MPI_Comm, mpi, INPUT, MPI_COMM_WORLD);
    ADD_SLOT(GridT, grid, INPUT, REQUIRED);
    ADD_SLOT(Vec3d, out, OUTPUT, DocString("Sum[f_i]"));

    inline std::string documentation() const final {
      return R"EOF(
        This operator returns out = Sum_{particles i}(f_i).
        Remark: do not hesite to use rebind to rename the output variable

        YAML example:

          opex:
            rebind:
              out: total_force
            body:
              - sum_forces
        )EOF";
    }
    
  public:
    inline void execute() final {
      
      // Reduce over the subdomain
      ForceValue value = { Vec3d{0.0, 0.0, 0.0 }};
      ReduceForceFunctor func;
      reduce_cell_particles(*grid, false, func, value, reduce_field_set, parallel_execution_context() );
      
      // Reduce over MPI processes
      double local[3] = {value.f_tot.x, value.f_tot.y, value.f_tot.z};
      double global[3] = {0.0, 0.0, 0.0};
      MPI_Allreduce(&local, &global, 3, MPI_DOUBLE, MPI_SUM, *mpi);
      Vec3d sum_forces = {global[0], global[1], global[2]};
      *out = sum_forces;
      static constexpr double conv  = EXASTAMP_CONST_QUANTITY( 1./eV );
      ldbg << " Total force in system (no ghosts) = (" << sum_forces * conv << ") eV/ang" << std::endl;
    }
  };
  
  // === register factories ===
  ONIKA_AUTORUN_INIT(sum_forces) {
    OperatorNodeFactory::instance()->register_factory("sum_forces", make_grid_variant_operator<SumForces>);
  }
  
}  // namespace exaDEM
