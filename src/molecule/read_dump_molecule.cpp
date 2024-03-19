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
#include <exaStamp/molecule/molecule_species.h>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ReadDumpMolecule : public OperatorNode
  {
    ADD_SLOT( MPI_Comm    , mpi             , INPUT );
    ADD_SLOT( std::string , filename , INPUT );
    ADD_SLOT( long        , timestep      , INPUT , DocString{"Iteration number"} );
    ADD_SLOT( double      , physical_time , INPUT , DocString{"Physical time"} );

    ADD_SLOT( GridT       , grid     , INPUT_OUTPUT );
    ADD_SLOT( Domain      , domain   , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT_OUTPUT, MoleculeSpeciesVector{} , DocString{"Molecule descriptions"} );
    ADD_SLOT(double , bond_max_dist     , INPUT_OUTPUT , 0.0 , DocString{"molecule bond max distance, in physical space"} );
    ADD_SLOT(double , bond_max_stretch  , INPUT_OUTPUT , 0.0 , DocString{"fraction of bond_max_dist."} );

  public:
    inline void execute () override final
    {
      using DumpFieldSet = FieldSet< field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_charge, field::_virial, field::_id, field::_idmol, field::_cmol, field::_type >;
      std::string file_name = data_file_path( *filename );

      AtomDumpFilter<GridT,DumpFieldSet,decltype(ldbg),MoleculeOptionalHeaderIO> dump_filter = { *species, ldbg , { *bond_max_dist , *bond_max_stretch , nullptr } };
      if( molecules.has_value() ) dump_filter.optional_header_io.m_molecules = molecules.get_pointer();

      exanb::read_dump( *mpi, ldbg, *grid, *domain, *physical_time, *timestep, file_name, DumpFieldSet{} , dump_filter );
    }
  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "read_dump_molecule" , make_grid_variant_operator<ReadDumpMolecule> );
  }

}

