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
    using BoolVector = std::vector<bool>;
    using StringVector = std::vector<std::string>;

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
    ADD_SLOT(double , bond_max_stretch  , INPUT_OUTPUT , 1.0 , DocString{"fraction of bond_max_dist."} );

    ADD_SLOT( double      , scale_cell_size , INPUT ,OPTIONAL , DocString{"if set, change cell size stored in file by scaling it with given factor"} );
    ADD_SLOT( BoolVector  , periodic        , INPUT ,OPTIONAL , DocString{"if set, overrides domain's periodicity stored in file with this value"}  );
    ADD_SLOT( StringVector, mirror          , INPUT ,OPTIONAL , DocString{"if set, overrides domain's boundary mirror flags in file with provided values"}  );
    ADD_SLOT( bool        , expandable      , INPUT ,OPTIONAL , DocString{"if set, override domain expandability stored in file"} );
    ADD_SLOT( AABB        , bounds          , INPUT ,OPTIONAL , DocString{"if set, override domain's bounds, filtering out particles outside of overriden bounds"} );
    ADD_SLOT( bool        , shrink_to_fit   , INPUT ,OPTIONAL , DocString{"if set to true and bounds was wpecified, try to reduce domain's grid size to the minimum size enclosing fixed bounds"} );

  public:
    inline void execute () override final
      {
      using DumpFieldSet = FieldSet< field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz, field::_charge, field::_virial, field::_id, field::_idmol, field::_cmol, field::_type >;
      using MolIOExt = MoleculeOptionalHeaderIO<decltype(ldbg)>;
      std::string file_name = data_file_path( *filename );

      double reader_bond_max_dist = 0.0;
      double reader_bond_max_stretch = 0.0;
      MolIOExt molecule_io = {
        reader_bond_max_dist ,
        reader_bond_max_stretch ,
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
      
      if( scale_cell_size.has_value() )
      {
        dump_filter.scale_cell_size = *scale_cell_size;
        ldbg << "force cell size scaling to "<<dump_filter.scale_cell_size<<std::endl;
      }
      if( periodic.has_value() )
      {
        dump_filter.override_periodicity = true;
        if( periodic->size() >= 1 ) dump_filter.periodic_x = periodic->at(0);
        if( periodic->size() >= 2 ) dump_filter.periodic_y = periodic->at(1);
        if( periodic->size() >= 3 ) dump_filter.periodic_z = periodic->at(2);
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
      if( mirror.has_value() )
      {
        dump_filter.override_mirroring = true;
        for(auto m : *mirror)
        {
          if( exanb::str_tolower(m) == "x-" ) { dump_filter.mirror_x_min=true; }
          if( exanb::str_tolower(m) == "x+" ) { dump_filter.mirror_x_max=true; }
          if( exanb::str_tolower(m) == "x" )  { dump_filter.mirror_x_min=true; dump_filter.mirror_x_max=true; }
          if( exanb::str_tolower(m) == "y-" ) { dump_filter.mirror_y_min=true; }
          if( exanb::str_tolower(m) == "y+" ) { dump_filter.mirror_y_max=true; }
          if( exanb::str_tolower(m) == "y" )  { dump_filter.mirror_y_min=true; dump_filter.mirror_y_max=true; }
          if( exanb::str_tolower(m) == "z-" ) { dump_filter.mirror_z_min=true; }
          if( exanb::str_tolower(m) == "z+" ) { dump_filter.mirror_z_max=true; }
          if( exanb::str_tolower(m) == "z" )  { dump_filter.mirror_z_min=true; dump_filter.mirror_z_max=true; }
        }
      }

      exanb::read_dump( *mpi, ldbg, *grid, *domain, *physical_time, *timestep, file_name, DumpFieldSet{} , dump_filter );
      
      ldbg << "--- Molecule compute parameters ---" << std::endl;
      molecule_io.print( ldbg , *species );

      MPI_Allreduce( MPI_IN_PLACE, &reader_bond_max_dist , 1 , MPI_DOUBLE , MPI_MAX , *mpi );
      MPI_Allreduce( MPI_IN_PLACE, &reader_bond_max_stretch , 1 , MPI_DOUBLE , MPI_MAX , *mpi );
      *bond_max_dist = std::max( *bond_max_dist , reader_bond_max_dist );
      *bond_max_stretch = std::max( *bond_max_stretch , reader_bond_max_stretch );
      ldbg << "bond max dist = "<< *bond_max_dist << " , bond_max_stretch = " << *bond_max_stretch << std::endl;
    }
  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "read_dump_molecule" , make_grid_variant_operator<ReadDumpMolecule> );
  }

}

