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

#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/io/write_xyz.h>
#include <exanb/compute/field_combiners.h>
#include <exaStamp/compute/field_combiners.h>


namespace exaStamp
{

  template<class GridT >
  class WriteXYZ : public OperatorNode
  {    
    using StringList = std::vector<std::string>;
    using has_field_type_t = typename GridT:: template HasField < field::_type >;
    static constexpr bool has_field_type = has_field_type_t::value;
    using KineticEnergyCombiner = std::conditional_t< has_field_type , MultimatKineticEnergyCombiner , MonomatKineticEnergyCombiner >;
    using MassCombiner = std::conditional_t< has_field_type , MultimatMassCombiner , MonomatMassCombiner >;
    using MomentumCombiner = std::conditional_t< has_field_type , MultimatMomentumCombiner , MonomatMomentumCombiner >;
    using KineticEnergyTensorCombiner = std::conditional_t< has_field_type , MultimatKineticEnergyTensorCombiner , MonomatKineticEnergyTensorCombiner >;
        
    ADD_SLOT( MPI_Comm        , mpi      , INPUT );
    ADD_SLOT( GridT           , grid     , INPUT );
    ADD_SLOT( Domain          , domain   , INPUT );
    ADD_SLOT( bool            , ghost    , INPUT , false );
    ADD_SLOT( std::string     , filename , INPUT , "output"); // default value for backward compatibility
    ADD_SLOT( StringList      , fields   , INPUT , StringList({".*"}) , DocString{"List of regular expressions to select fields to project"} );
    ADD_SLOT( ParticleSpecies , species  , INPUT , REQUIRED );
      
    template<class... fid>
    inline void execute_on_field_set( FieldSet<fid...> ) 
    {
      int rank=0;
      MPI_Comm_rank(*mpi, &rank);
      
      VelocityNorm2Combiner vnorm2 = {};
      ProcessorRankCombiner processor_id = { {rank} };
      KineticEnergyCombiner mv2 = { { species->data() , 0 } };
      MassCombiner mass = { { species->data() , 0 } };
      MomentumCombiner momentum = { { species->data() , 0 } };
      
      const auto& sp = *species;
      auto particle_type_func = [&sp](auto cells, size_t c, size_t pos) -> const char*
      {
        using has_field_type_t = typename GridT:: template HasField < field::_type >;
        static constexpr bool has_field_type = has_field_type_t::value;
        if constexpr ( has_field_type )
        {
          return sp[ cells[c][field::type][pos] ].m_name;
        }
        else
        {
          return "XX";
        }
      };
      
      write_xyz_details::write_xyz_grid_fields( ldbg, *mpi, *grid, *domain, *fields, *filename, particle_type_func, *ghost, vnorm2, processor_id, mv2, mass, momentum, onika::soatl::FieldId<fid>{} ... );
    }

    public:
    inline void execute() override
    {
      execute_on_field_set(grid->field_set);
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "write_xyz", make_grid_variant_operator< WriteXYZ > );
  }

}
