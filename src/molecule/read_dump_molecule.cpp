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
#include <exaStamp/molecule/molecule_optional_header_io.h>

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

    ADD_SLOT( GridT           , grid     , INPUT_OUTPUT );
    ADD_SLOT( Domain          , domain   , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );
    
    ADD_SLOT( MoleculeSpeciesVector            , molecules                , INPUT_OUTPUT, MoleculeSpeciesVector{} , DocString{"Molecule descriptions"} );
    ADD_SLOT( BondsPotentialParameters         , potentials_for_bonds     , INPUT_OUTPUT, BondsPotentialParameters{} );
    ADD_SLOT( BendsPotentialParameters         , potentials_for_angles    , INPUT_OUTPUT, BendsPotentialParameters{} );
    ADD_SLOT( TorsionsPotentialParameters      , potentials_for_torsions  , INPUT_OUTPUT, TorsionsPotentialParameters{} );
    ADD_SLOT( ImpropersPotentialParameters     , potentials_for_impropers , INPUT_OUTPUT, ImpropersPotentialParameters{} );    
    ADD_SLOT( LJExp6RFMultiParms               , potentials_for_pairs     , INPUT_OUTPUT, LJExp6RFMultiParms{} );
    ADD_SLOT( IntramolecularPairWeighting      , mol_pair_weights         , INPUT_OUTPUT, IntramolecularPairWeighting{} );

    ADD_SLOT(double , bond_max_dist     , INPUT_OUTPUT , 0.0 , DocString{"molecule bond max distance, in physical space"} );
    ADD_SLOT(double , bond_max_stretch  , INPUT_OUTPUT , 0.0 , DocString{"fraction of bond_max_dist."} );

  public:
    inline void execute () override final
    {
      using DumpFieldSet = FieldSet< field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_charge, field::_virial, field::_id, field::_idmol, field::_cmol, field::_type >;
      using MolIOExt = MoleculeOptionalHeaderIO<decltype(ldbg)>;
      std::string file_name = data_file_path( *filename );

      MolIOExt molecule_io = {
        *bond_max_dist ,
        *bond_max_stretch ,
        molecules.get_pointer() ,
        ldbg,
        potentials_for_bonds.get_pointer() ,
        potentials_for_angles.get_pointer() ,
        potentials_for_torsions.get_pointer() ,
        potentials_for_impropers.get_pointer() ,
        potentials_for_pairs.get_pointer() ,
        mol_pair_weights.get_pointer() };

      AtomDumpFilter<GridT,DumpFieldSet,decltype(ldbg),MolIOExt> dump_filter = { *species, ldbg , molecule_io };
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

