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

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_charge, field::_virial, field::_id, field::_idmol, field::_cmol, field::_type >
    >
  class ReadDumpMolecule : public OperatorNode
  {
    ADD_SLOT( MPI_Comm    , mpi             , INPUT );
    ADD_SLOT( GridT       , grid     , INPUT );
    ADD_SLOT( Domain      , domain   , INPUT );
    ADD_SLOT( std::string , filename , INPUT );
    ADD_SLOT( long        , timestep      , INPUT , DocString{"Iteration number"} );
    ADD_SLOT( double      , physical_time , INPUT , DocString{"Physical time"} );

    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );

  public:
    inline void execute () override final
    {
      using DumpFieldSet = FieldSet< field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_charge, field::_virial, field::_id, field::_idmol, field::_cmol, field::_type >;
      std::string file_name = data_file_path( *filename );
      exanb::read_dump( *mpi, ldbg, *grid, *domain, *physical_time, *timestep, file_name, DumpFieldSet{} , make_atom_dump_filter(*grid,*species,ldbg,DumpFieldSet{}) );
    }
  };

  template<class GridT> using ReadDumpMoleculeTmpl = ReadDumpMolecule<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "read_dump_molecule" , make_grid_variant_operator<ReadDumpMoleculeTmpl> );
  }

}

