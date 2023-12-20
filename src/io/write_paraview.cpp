#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/log.h>
#include <exanb/core/domain.h>

#include <exanb/compute/field_combiners.h>
#include <exaStamp/compute/field_combiners.h>
#include <exaStamp/mechanical/cell_particles_local_mechanical_metrics.h>

#include <exanb/io/vtk_writer.h>
#include <exanb/io/vtk_writer_binary.h>
#include <exanb/io/vtk_writer_ascii.h>
#include <exanb/io/write_paraview.h>

#include <mpi.h>
#include <string>
#include <regex>

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT>
  class ParaviewWriter : public OperatorNode
  {
    using StringList = std::vector<std::string>;

    using has_field_type_t = typename GridT:: template HasField < field::_type >;
    static constexpr bool has_field_type = has_field_type_t::value;
    using KineticEnergyCombiner = std::conditional_t< has_field_type , MultimatKineticEnergyCombiner , MonomatKineticEnergyCombiner >;
    using MassCombiner = std::conditional_t< has_field_type , MultimatMassCombiner , MonomatMassCombiner >;
    using MomentumCombiner = std::conditional_t< has_field_type , MultimatMomentumCombiner , MonomatMomentumCombiner >;
    using KineticEnergyTensorCombiner = std::conditional_t< has_field_type , MultimatKineticEnergyTensorCombiner , MonomatKineticEnergyTensorCombiner >;

    ADD_SLOT( MPI_Comm    , mpi                  , INPUT );
    ADD_SLOT( ParticleSpecies, species           , INPUT , REQUIRED );    
    ADD_SLOT( GridT       , grid                 , INPUT );
    ADD_SLOT( Domain      , domain               , INPUT );
    ADD_SLOT( bool        , binary_mode          , INPUT , true);
    ADD_SLOT( bool        , write_box            , INPUT , true);
    ADD_SLOT( bool        , write_external_box   , INPUT , false);
    ADD_SLOT( bool        , write_ghost          , INPUT , false);
    ADD_SLOT( std::string , compression          , INPUT , "default");
    ADD_SLOT( std::string , filename             , INPUT , "output"); // default value for backward compatibility
    ADD_SLOT( StringList  , fields               , INPUT , StringList({".*"}) , DocString{"List of regular expressions to select fields to project"} );

    ADD_SLOT( GridParticleLocalMechanicalMetrics, local_mechanical_data , INPUT, OPTIONAL );

    template<class... GridFields>
    inline void execute_on_fields( const GridFields& ... grid_fields) 
    {
      {
        int s=0;
        ldbg << "Paraview writer available fields:";
        ( ... , ( ldbg<< (((s++)==0)?' ':',') <<grid_fields.short_name() ) ) ;
        ldbg << std::endl;
      }
    
      const auto& flist = *fields;
      auto field_selector = [&flist] ( const std::string& name ) -> bool { for(const auto& f:flist) if( std::regex_match(name,std::regex(f)) ) return true; return false; } ;

      GridParticleMechanicalAccessor<typename GridT::CellParticles *> gridacc = { grid->cells() , nullptr };
      if( local_mechanical_data.has_value() ) gridacc.mech_data = local_mechanical_data->data();
      
      ParaviewWriteTools::write_particles(ldbg,*mpi,*grid,gridacc,*domain,*filename,field_selector,*compression,*binary_mode,*write_box,*write_external_box,*write_ghost, grid_fields ... );
    }

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
      //KineticEnergyTensorCombiner mv2tensor = { { species->data() , 0 } }; // too big for full outputs

      execute_on_fields( processor_id, vnorm2, mv2, mass, momentum,
                         mechanical::defgrad , mechanical::greenlag , mechanical::rot , mechanical::stretch ,
                         onika::soatl::FieldId<fid>{} ... );
    }

  public:
    inline void execute() override final
    {
      execute_on_field_set(grid->field_set);
    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "write_paraview",make_grid_variant_operator<ParaviewWriter>);
  }

}
