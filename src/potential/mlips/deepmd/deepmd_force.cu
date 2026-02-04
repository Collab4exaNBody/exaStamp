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

#include <memory>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/source_term.h>

#include <mpi.h>
#include <iomanip>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz , field::_fx, field::_fy, field::_fz >
    >
  class DeepMDForce : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi          , INPUT , MPI_COMM_WORLD );
    
    ADD_SLOT( ParticleSpecies, species      , INPUT , REQUIRED );
    ADD_SLOT( double         , rcut_max     , INPUT_OUTPUT , 0.0 );  // neighborhood distance, in grid space

    ADD_SLOT( long           , grid_subdiv  , INPUT , 3 );
    ADD_SLOT( GridCellValues , grid_cell_values      , INPUT_OUTPUT );

    ADD_SLOT( GridT  , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain , domain       , INPUT , REQUIRED );
    
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {

      std::cout << "DeepMD force computation" << std::endl;
      
    }

  };

  template<class GridT> using DeepMDForceTmpl = DeepMDForce<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(deepmd_force)
  {
   OperatorNodeFactory::instance()->register_factory("deepmd_force", make_grid_variant_operator< DeepMDForceTmpl > );
  }

}
