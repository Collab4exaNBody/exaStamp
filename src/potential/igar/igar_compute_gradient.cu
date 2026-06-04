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

#include "igar_compute_gradient.h"

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
  class IgarComputeGradient : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                  , mpi                 , INPUT        , REQUIRED );
    ADD_SLOT( ParticleSpecies           , species             , INPUT        , REQUIRED );
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( GridCellValues            , grid_cell_values    , INPUT_OUTPUT );

  public:
    
    inline std::string documentation() const override final
    {
      return R"EOF(
Precomputes the spatial gradient of the IGAR energy field igar_ep using
second-order centered finite differences at every subcell of the grid.

Produces three new fields in grid_cell_values:
  igar_dEdx, igar_dEdy, igar_dEdz

each at the same subdiv resolution as igar_ep. These fields are consumed
by igar_force_from_gradient to apply forces to particles.

This operator should typically be called once in setup_system when igar_ep
is a static reference field (e.g. loaded from a VTK file). Calling it every
timestep is correct but expensive: it sweeps all subdiv^3 subcells per cell
regardless of particle count.

Example (LJ_igar.msp — gradient precomputed once at setup):

  setup_system:
    - read_cell_values:
        field_name: "igar_ep"
        field_dim: 1
        grid_subdiv: 64
        file_name: "igar_test.vtk"
    - igar_compute_gradient

  compute_force:
    - lj_compute_force
    - igar_force_from_gradient: { energy_factor: 20000.0 }
)EOF";
    }

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

  template<class GridT> using IgarComputeGradientTmpl = IgarComputeGradient<GridT>;

  ONIKA_AUTORUN_INIT(igar_compute_gradient)
  {
    OperatorNodeFactory::instance()->register_factory("igar_compute_gradient", make_grid_variant_operator<IgarComputeGradientTmpl>);
  }

}
