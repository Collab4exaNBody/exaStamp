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

#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/physics/units.h>
#include <exanb/io/write_xyz.h>
#include <exanb/compute/field_combiners.h>
#include <exaStamp/compute/field_combiners.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/mechanical/cell_particles_local_mechanical_metrics.h>
#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>
#include <exaStamp/mechanical/average_local_field.h>

namespace exaStamp
{

  template<class GridT >
  class WriteXYZ : public OperatorNode
  {    
    using StringList = std::vector<std::string>;
    using StringMap = std::map<std::string,std::string>;
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
    ADD_SLOT( StringMap       , units    , INPUT , StringMap( { {"position","ang"} , {"velocity","m/s"} , {"force","m/s/g"} } ) , DocString{"Units to be used for specific fields."} );
    ADD_SLOT( ParticleSpecies , species  , INPUT , REQUIRED );
    ADD_SLOT( double          , physical_time  , INPUT , 0.0 );
    ADD_SLOT( GridParticleLocalStructuralMetrics, local_structural_data , INPUT );
      
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
        //else // commented out to avoid intel compiler warnings
        //{
          return "XX";
        //}
      };

      std::unordered_map<std::string,double> conv_scale;
      for(const auto& umap : *units)
      {
        const auto s = "1.0 " + umap.second;
        bool conv_ok = false;
        auto q = onika::physics::quantity_from_string( s , conv_ok );
        if( ! conv_ok ) { fatal_error() << "Failed to parse unit string '"<<s<<"'"<<std::endl; }
        conv_scale[umap.first] = q.convert();
      }
      write_xyz_details::DefaultFieldFormatter field_formatter = { *units , conv_scale };
      
      PositionVec3Combiner position = {};
      VelocityVec3Combiner velocity = {};
      ForceVec3Combiner    force    = {};
      
      // position is required, and must be the first field
      StringList flist = { "position" };
      for(const auto& f : *fields) { if( f != "position" ) flist.push_back(f); }
      
      // property name for position must be 'Position'
      field_formatter.m_field_name_map["position"] = "pos";
      field_formatter.m_field_name_map["velocity"] = "velo";
      field_formatter.m_field_name_map["virial"] = "stress";
      field_formatter.m_field_name_map["orient"] = "orientation";
      field_formatter.m_field_name_map["angmom"] = "angular_momentum";
      field_formatter.m_field_name_map["couple"] = "torque";
      field_formatter.m_field_name_map["idmol"] = "molecule";

      // structural fields
      const CellParticleLocalStructuralMetrics * __restrict__ struct_data = nullptr;      
      if( local_structural_data.has_value() ) struct_data = local_structural_data->data();
      auto nneigh = structural_field(struct_data,field::numneighbors);
      auto csp    = structural_field(struct_data,field::csp);
      auto local_field = grid->field_accessor( field::local_field );
        
      write_xyz_details::write_xyz_grid_fields( ldbg, *mpi, *grid, *domain, flist, *filename, particle_type_func, field_formatter, *ghost, *physical_time
                                                , position, velocity, force, processor_id, vnorm2, mv2, mass, momentum, nneigh, csp, local_field, onika::soatl::FieldId<fid>{} ... );
    }

    public:
    inline void execute() override
    {
      execute_on_field_set(grid->field_set);
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(write_xyz)
  {
    OperatorNodeFactory::instance()->register_factory( "write_xyz", make_grid_variant_operator< WriteXYZ > );
  }

}
