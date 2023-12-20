#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>

#include <exaStamp/compute/physics_functors.h>
#include <exanb/grid_cell_particles/particle_cell_projection.h>
#include <exanb/grid_cell_particles/grid_particle_field_accessor.h>

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
  class AtomCellProjection : public OperatorNode
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
    
    template<class... GridFields>
    inline void execute_on_fields( const GridFields& ... grid_fields) 
    {
      using namespace ParticleCellProjectionTools;

      {
        int s=0;
        ldbg << "Atom cell projection available fields:";
        ( ... , ( ldbg<< (((s++)==0)?' ':',') <<grid_fields.short_name() ) ) ;
        ldbg << std::endl;
      }

      const auto& flist = *fields;
      auto field_selector = [&flist] ( const std::string& name ) -> bool { for(const auto& f:flist) if( std::regex_match(name,std::regex(f)) ) return true; return false; } ;

      // create cell value fields
      std::vector<AddCellFieldInfo> fields_to_add;
      CollectCellValueFieldToAdd collect_fields = {*grid_cell_values,fields_to_add,field_selector,*grid_subdiv};
      apply_grid_fields   ( *grid, collect_fields , grid_fields ... );
      ldbg << "add "<< fields_to_add.size() << " cell fields :"<<std::endl;
      for(const auto& f:fields_to_add) ldbg << "\t" << f.m_name<<std::endl;
      grid_cell_values->add_fields( fields_to_add );

      // project particle quantities to cells
      //using ParticleAcessor = GridParticleFieldAccessor<typename GridT::CellParticles *>;
      using ParticleAcessor = GridParticleMechanicalAccessor<typename GridT::CellParticles *>;
      ParticleAcessor gridacc = { grid->cells() , nullptr };
      if( local_mechanical_data.has_value() ) gridacc.mech_data = local_mechanical_data->data();

      ProjectCellValueField<ParticleAcessor> project_fields = { gridacc , *grid_cell_values,field_selector,*splat_size,*grid_subdiv};
      apply_grid_fields   ( *grid, project_fields , grid_fields ... );
    }

    template<class... fid>
    inline void execute_on_field_set( FieldSet<fid...> ) 
    {
      int rank=0;
      MPI_Comm_rank(*mpi, &rank);

      //VelocityNorm2Combiner vnorm2 = {};
      VelocityNormCombiner vnorm = {};
      ParticleCountCombiner count = {};
      ProcessorRankCombiner processor_id = { {rank} };
      KineticEnergyCombiner mv2 = { { species->data() , 0 } };
      MassCombiner mass = { { species->data() , 0 } };
      MomentumCombiner momentum = { { species->data() , 0 } };
      KineticEnergyTensorCombiner mv2tensor = { { species->data() , 0 } };

      execute_on_fields( count, processor_id, vnorm, mv2, mass, momentum, mv2tensor,
                         mechanical::defgrad , mechanical::greenlag , mechanical::rot , mechanical::stretch ,
                         onika::soatl::FieldId<fid>{} ... );
    }

  public:

    // -----------------------------------------------
    inline void execute ()  override final
    {
      if( grid->number_of_cells() > 0 )
      {
        execute_on_field_set(grid->field_set);
      }
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(project atom quantities onto a regular grid)EOF";
    }    

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory("atom_cell_projection", make_grid_variant_operator< AtomCellProjection > );
  }

}
