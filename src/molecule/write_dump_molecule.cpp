#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/log.h>
#include <exanb/core/domain.h>

#include <iostream>
#include <fstream>
#include <string>

#include <exanb/io/sim_dump_writer.h>
#include <exaStamp/io/atom_dump_filter.h>
#include <exaStamp/molecule/molecule_species.h>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_charge, field::_virial, field::_id, field::_idmol, field::_cmol, field::_type >
    >
  class WriteDumpMolecule : public OperatorNode
  {
    ADD_SLOT( MPI_Comm    , mpi             , INPUT );
    ADD_SLOT( GridT       , grid     , INPUT , REQUIRED );
    ADD_SLOT( Domain      , domain   , INPUT , REQUIRED );
    ADD_SLOT( std::string , filename , INPUT , "molecule.dump" );
    ADD_SLOT( long        , timestep      , INPUT , DocString{"Iteration number"} );
    ADD_SLOT( double      , physical_time , INPUT , DocString{"Physical time"} );
    ADD_SLOT( long        , compression_level , INPUT , 6 , DocString{"Zlib compression level"} );
    ADD_SLOT( long        , max_part_size , INPUT , -1 , DocString{"Maximum file partition size. set -1 for system default value"} );
    ADD_SLOT( ParticleSpecies , species , INPUT , REQUIRED );
    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT, OPTIONAL , DocString{"Molecule descriptions"} );

    ADD_SLOT(double , bond_max_dist     , INPUT , 0.0 , DocString{"molecule bond max distance, in physical space"} );
    ADD_SLOT(double , bond_max_stretch  , INPUT , 0.0 , DocString{"fraction of bond_max_dist."} );

  public:
    inline void execute () override final
    {
      using DumpFieldSet = FieldSet< field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_charge, field::_virial, field::_id, field::_idmol, field::_cmol, field::_type >; 
      size_t mps = MpiIO::DEFAULT_MAX_FILE_SIZE;
      if( *max_part_size > 0 ) mps = *max_part_size;

      AtomDumpFilter<GridT,DumpFieldSet,decltype(ldbg),MoleculeOptionalHeaderIO> dump_filter = { *species, ldbg , { *bond_max_dist , *bond_max_stretch , nullptr } };
      if( molecules.has_value() ) dump_filter.optional_header_io.m_molecules = molecules.get_pointer();
      
      exanb::write_dump( *mpi, ldbg, *grid, *domain, *physical_time, *timestep, *filename, *compression_level, DumpFieldSet{} , dump_filter , mps );
    }
  };

  template<class GridT> using WriteDumpMoleculeTmpl = WriteDumpMolecule<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "write_dump_molecule" , make_grid_variant_operator<WriteDumpMoleculeTmpl> );
  }

}

