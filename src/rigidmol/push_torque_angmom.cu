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





#include <string>
#include <numeric>

#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/log.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/parallel/random.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/physics/units.h>

#include <onika/cuda/cuda.h>
#include <exanb/compute/compute_cell_particles.h>

//#include "quaternion_rotation.h"
#include <onika/math/quaternion_operators.h>

namespace exaStamp
{
  using namespace exanb;

  struct PushTorqueAngmomComputeFunc
  {
    const double dt;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Vec3d& couple, Vec3d& angmom ) const
    {
      angmom += couple * dt / 2.0;
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits< exaStamp::PushTorqueAngmomComputeFunc >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_angmom, field::_couple >
    >
  class PushTorqueAngmomRigidMol : public OperatorNode
  {
    //ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( double          , dt           , INPUT, REQUIRED );

    static constexpr FieldSet< field::_couple , field::_angmom > compute_field_set{};

  public:
    inline void execute () override final
    {
      compute_cell_particles( *grid , false , PushTorqueAngmomComputeFunc{*dt} , compute_field_set , parallel_execution_context() );
    }

  };

  // === register factories ===
  template<class GridT> using PushTorqueAngmomRigidMolTmpl = PushTorqueAngmomRigidMol<GridT>;

  ONIKA_AUTORUN_INIT(push_torque_angmom)
  {
    OperatorNodeFactory::instance()->register_factory("push_torque_angmom", make_grid_variant_operator< PushTorqueAngmomRigidMolTmpl >);
  }

}
