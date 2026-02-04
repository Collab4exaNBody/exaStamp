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
#include <exanb/analytics/particle_cell_projection.h>
#include <exanb/core/grid_particle_field_accessor.h>

#include <exaStamp/compute/physics_functors.h>
#include <exanb/compute/field_combiners.h>
#include <exaStamp/compute/field_combiners.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/mechanical/cell_particles_local_mechanical_metrics.h>

#include <mpi.h>
#include <regex>

namespace exaStamp
{
  using namespace exanb;

  template< class GridT >
  class MechanicalCellProjection : public OperatorNode
  {    
    using StringList = std::vector<std::string>;
    using has_field_type_t = typename GridT:: template HasField < field::_type >;
    static constexpr bool has_field_type = has_field_type_t::value;

    using KineticEnergyCombiner = std::conditional_t< has_field_type , MultimatKineticEnergyCombiner , MonomatKineticEnergyCombiner >;
    using MassCombiner = std::conditional_t< has_field_type , MultimatMassCombiner , MonomatMassCombiner >;
    using MomentumCombiner = std::conditional_t< has_field_type , MultimatMomentumCombiner , MonomatMomentumCombiner >;
    using KineticEnergyTensorCombiner = std::conditional_t< has_field_type , MultimatKineticEnergyTensorCombiner , MonomatKineticEnergyTensorCombiner >;

    ADD_SLOT( MPI_Comm    , mpi             , INPUT );
    ADD_SLOT( ParticleSpecies, species                    , INPUT , REQUIRED );    
    ADD_SLOT( GridT          , grid              , INPUT , REQUIRED );
    ADD_SLOT( double         , splat_size        , INPUT , REQUIRED );
    ADD_SLOT( StringList  , fields            , INPUT , StringList({".*"}) , DocString{"List of regular expressions to select fields to project"} );

    ADD_SLOT( GridParticleLocalMechanicalMetrics, local_mechanical_data , INPUT, OPTIONAL );

    ADD_SLOT( long           , grid_subdiv       , INPUT_OUTPUT , 1 );
    ADD_SLOT( GridCellValues , grid_cell_values  , INPUT_OUTPUT );
    
  public:

    // -----------------------------------------------
    inline void execute ()  override final
    {
      using namespace ParticleCellProjectionTools;

      if( grid->number_of_cells() == 0 ) return;
        
      int rank=0;
      MPI_Comm_rank(*mpi, &rank);

      VelocityNormCombiner vnorm = {};
      ParticleCountCombiner count = {};
      KineticEnergyCombiner mv2 = { { species->data() , 0 } };
      MassCombiner mass = { { species->data() , 0 } };
      MomentumCombiner momentum = { { species->data() , 0 } };
      KineticEnergyTensorCombiner mv2tensor = { { species->data() , 0 } };

      // mechanical fields
      const CellParticleLocalMechanicalMetrics * __restrict__ mech_data = nullptr;
      if( local_mechanical_data.has_value() ) mech_data = local_mechanical_data->data();
      
      auto defgrad = mechanical_field(mech_data,field::defgrad);
      auto greenlag = mechanical_field(mech_data,field::greenlag);
      auto rot     = mechanical_field(mech_data,field::rot);
      auto stretch = mechanical_field(mech_data,field::stretch);
      
      auto proj_fields = make_field_tuple_from_field_set( grid->field_set, count, vnorm, mv2, mass, momentum, mv2tensor, defgrad, greenlag, rot, stretch );
      auto field_selector = [flist = *fields] ( const std::string& name ) -> bool { for(const auto& f:flist) if( std::regex_match(name,std::regex(f)) ) return true; return false; } ;
      project_particle_fields_to_grid( ldbg, *grid, *grid_cell_values, *grid_subdiv, *splat_size, field_selector, proj_fields );
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(project atom quantities onto a regular grid)EOF";
    }    

  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(mechanical_cell_projection)
  {
    OperatorNodeFactory::instance()->register_factory("mechanical_cell_projection", make_grid_variant_operator< MechanicalCellProjection > );
  }

}
