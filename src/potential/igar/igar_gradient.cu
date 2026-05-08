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
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>

#include <exaStamp/compute/physics_functors.h>

#include <exanb/core/grid_particle_field_accessor.h>

#include <exaStamp/compute/physics_functors.h>
#include <exanb/compute/field_combiners.h>
#include <exaStamp/compute/field_combiners.h>
#include <exaStamp/particle_species/particle_specie.h>

//#include <exanb/analytics/particle_cell_projection.h>

#include "energy_gradient.h"

#include <memory>
#include <vector>
#include <mpi.h>

namespace exaStamp
{

  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class IgarGradient : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                  , mpi                 , INPUT        , REQUIRED );
    ADD_SLOT( ParticleSpecies           , species             , INPUT        , REQUIRED );
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( GridCellValues            , grid_cell_values    , INPUT_OUTPUT );

  public:

    inline void execute() override final
    {
      using namespace ParticleCellProjectionTools;

      if( grid->number_of_cells() == 0 ) return;
        
      int rank=0;
      MPI_Comm_rank(*mpi, &rank);

      ldbg << "Igar gradient from grid_cell_values" << std::endl;
      compute_energy_gradient( ldbg, *grid, *grid_cell_values );      
      ldbg << "DONE." << std::endl;
    }

  };

  template<class GridT> using IgarGradientTmpl = IgarGradient<GridT>;

  ONIKA_AUTORUN_INIT(igar_gradient)
  {
    OperatorNodeFactory::instance()->register_factory("igar_gradient", make_grid_variant_operator<IgarGradientTmpl>);
  }

}
