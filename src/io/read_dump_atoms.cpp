#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/log.h>
#include <exanb/core/domain.h>
#include <exanb/core/file_utils.h>

#include <iostream>
#include <fstream>
#include <string>

#include <exanb/io/sim_dump_reader.h>
#include <exaStamp/io/atom_dump_filter.h>

namespace exaStamp
{
  using namespace exanb;

  using BoolVector = std::vector<bool>;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_id, field::_type >
    >
  class ReadDumpAtoms : public OperatorNode
  {
    ADD_SLOT( MPI_Comm    , mpi             , INPUT );
    ADD_SLOT( std::string , filename , INPUT );
    ADD_SLOT( long        , timestep      , INPUT , DocString{"Iteration number"} );
    ADD_SLOT( double      , physical_time , INPUT , DocString{"Physical time"} );

    ADD_SLOT( double      , scale_cell_size , INPUT ,OPTIONAL , DocString{"if set, change cell size stored in file by scaling it with given factor"} );
    ADD_SLOT( BoolVector  , periodicity     , INPUT ,OPTIONAL , DocString{"if set, overrides domain's periodicity stored in file with this value"}  );
    ADD_SLOT( bool        , expandable      , INPUT ,OPTIONAL , DocString{"if set, override domain expandability stored in file"} );
    ADD_SLOT( AABB        , bounds          , INPUT ,OPTIONAL , DocString{"if set, override domain's bounds, filtering out particle outside of overriden bounds"} );
    ADD_SLOT( bool        , shrink_to_fit   , INPUT ,OPTIONAL , DocString{"if set to true and bounds was wpecified, try to reduce domain's grid size to the minimum size enclosing fixed bounds"} );

    ADD_SLOT( GridT       , grid     , INPUT_OUTPUT );
    ADD_SLOT( Domain      , domain   , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , ParticleSpecies{} );

  public:
    inline void execute () override final
    {
      static constexpr FieldSet< field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_id, field::_type > dump_field_set={};
      std::string file_name = data_file_path( *filename );
      *physical_time = 0.0;
      *timestep = 0;
      auto dump_filter = make_atom_dump_filter(*grid,*species,ldbg,dump_field_set);
      if( scale_cell_size.has_value() )
      {
        dump_filter.scale_cell_size = *scale_cell_size;
        ldbg << "force cell size scaling to "<<dump_filter.scale_cell_size<<std::endl;
      }
      if( periodicity.has_value() )
      {
        dump_filter.override_periodicity = true;
        if( periodicity->size() >= 1 ) dump_filter.periodic_x = periodicity->at(0);
        if( periodicity->size() >= 2 ) dump_filter.periodic_y = periodicity->at(1);
        if( periodicity->size() >= 3 ) dump_filter.periodic_z = periodicity->at(2);
        ldbg << "force periodicity to ("<<std::boolalpha<<dump_filter.periodic_x<<","<<dump_filter.periodic_y<<","<<dump_filter.periodic_z<<")" <<std::endl;
      }
      if( expandable.has_value() )
      {
        ldbg << "force expandability to "<<std::boolalpha<< *expandable << std::endl;
        dump_filter.override_expandable = true;
        dump_filter.expandable = *expandable;
      }
      if( bounds.has_value() )
      {
        ldbg << "force domain bounds to "<< *bounds << std::endl;
        dump_filter.override_domain_bounds = true;
        dump_filter.domain_bounds = *bounds;
        if( shrink_to_fit.has_value() )
        {
          dump_filter.shrink_to_fit = *shrink_to_fit;
        }
      }
      
      exanb::read_dump( *mpi, ldbg, *grid, *domain, *physical_time, *timestep, file_name, dump_field_set , dump_filter );
    }
  };

  template<class GridT> using ReadDumpAtomsTmpl = ReadDumpAtoms<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "read_dump_atoms" , make_grid_variant_operator<ReadDumpAtomsTmpl> );
  }

}

