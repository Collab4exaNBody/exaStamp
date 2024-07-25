#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exanb/core/basic_types.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/log.h>
#include <exanb/core/particle_type_id.h>

#include <exaStamp/particle_species/particle_specie.h>

#include <exaStamp/molecule/mol_connectivity.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/molecule_compute_param.h>

#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/impropers_potentials_parameters.h>
#include <exaStamp/molecule/id_map.h>
#include <exanb/core/particle_id_codec.h>

#include <exaStamp/molecule/intramolecular_pair_weight.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>


#include <onika/oarray_stream.h>

#include <unordered_set>
#include <vector>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct AtomNode
  {
    int64_t id[5] = { -1 , -1 , -1 , -1 , -1 };
    uint64_t cell_particle = 0;
    int type = -1;
    inline bool operator < (const AtomNode & other) const
    {
      for(int i=0;i<5;i++)
      {
        if( id[i]<other.id[i] ) return true;
        else if( id[i]>other.id[i] ) return false;
      }
      return false;
    }
  };

  struct MoleculeAtoms
  {
    std::vector<AtomNode> atoms;
    int moltype = -1;
    inline bool operator < (const MoleculeAtoms & other) const
    {
      return atoms < other.atoms;
    }
  };

  template<class CellT, class FieldIdT, class FieldTypeT>
  struct MoleculeParser
  {
    CellT m_cells = {};
    FieldIdT m_field_cmol = {};
    FieldTypeT m_field_type = {};
    const IdMap& id_map;
    std::unordered_set< uint64_t > m_visited_atoms;
    
    inline void explore_molecule ( std::vector<AtomNode>& molecule_atoms , int64_t id )
    {
      if( m_visited_atoms.find(id) != m_visited_atoms.end() ) return;
      auto it = id_map.find(id);
      if( it != id_map.end() )
      {
        m_visited_atoms.insert( id );
        size_t cell_i=0, p_i=0;
        decode_cell_particle( it->second , cell_i , p_i );
        auto cmol = m_cells[cell_i][m_field_cmol][p_i];
        molecule_atoms.push_back( { { id , int64_t(cmol[0]) , int64_t(cmol[1]) , int64_t(cmol[2]) , int64_t(cmol[3]) } , encode_cell_particle(cell_i,p_i) , m_cells[cell_i][m_field_type][p_i] } );
        for(auto c : cmol)
        {
          if( c != std::numeric_limits<uint64_t>::max() ) explore_molecule( molecule_atoms , c );
        }
      }
    }
    
  };

  
  template< class GridT >
  class IntramolecularSetup : public OperatorNode
  {

    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , MPI_COMM_WORLD);

    ADD_SLOT( ParticleSpecies              , species           , INPUT , REQUIRED );
    ADD_SLOT( ParticleTypeMap              , particle_type_map , INPUT , REQUIRED );
    
    ADD_SLOT( BondsPotentialParameters     , potentials_for_bonds     , INPUT, OPTIONAL );
    ADD_SLOT( BendsPotentialParameters     , potentials_for_angles    , INPUT, OPTIONAL );
    ADD_SLOT( TorsionsPotentialParameters  , potentials_for_torsions  , INPUT, OPTIONAL );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT, OPTIONAL );
    ADD_SLOT( LJExp6RFMultiParms           , potentials_for_pairs , INPUT , LJExp6RFMultiParms{} );
    
    ADD_SLOT( IntramolecularPairWeighting  , mol_pair_weights  , INPUT , IntramolecularPairWeighting{} );

    ADD_SLOT( GridT                        , grid              , INPUT_OUTPUT);
    ADD_SLOT( MoleculeSpeciesVector        , molecules         , INPUT_OUTPUT , REQUIRED , DocString{"Molecule descriptions"} );
    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );

  public:
    inline bool is_sink() const override final { return true; }
    
    inline void execute ()  override final
    {
      //static constexpr MoleculeGenericFuncParam null_param = {0.0,0.0,0.0,0.0};

      const unsigned int nmol = molecules->m_molecules.size();
      ldbg << "Number of molecule species : "<<nmol<<std::endl;
      if( nmol == 1 && molecules->m_molecules.front().name().empty() ) molecules->m_molecules.front().set_name("MOL");
      for(unsigned int m=0;m<nmol;m++)
      {
        ldbg << "molecule #"<<m<<" : '"<<molecules->m_molecules.at(m).name()<<"'" << std::endl;
        if( ! molecules->m_molecules[m].has_connectivity() )
        {
          ldbg << "Update connectivity for molecule #"<<m<< std::endl;
          molecules->m_molecules[m].update_connectivity();
          molecules->m_molecules[m].print( ldbg , *species );
        }
      }

//      molecule_compute_parameters->m_molecules.assign( nmol , MoleculeComputeParams{} );   
      const auto & tmap = *particle_type_map;
      auto str2type = [&tmap]( const std::string& s ) -> int
        {
          auto it=tmap.find(s);
          if(it!=tmap.end()) return it->second;
          else return -1;
        };


      /********* global energy and virial correction ***************/
      const size_t n_cells = grid->number_of_cells();
      const unsigned int n_type_pairs = unique_pair_count( species->size() );      
      // count number of atoms of each type
      std::vector<long> type_count( species->size() , 0 );
      for(size_t ci=0;ci<n_cells;ci++)
      {
        if( ! grid->is_ghost_cell(ci) )
        {
          const auto& cell = grid->cell(ci);
          size_t n_particles = cell.size();
          for(size_t pi=0;pi<n_particles;pi++)
          {
            int t = cell[field::type][pi];
            ++ type_count[t];
          }
        }
      }
      MPI_Allreduce(MPI_IN_PLACE,type_count.data(), type_count.size(), MPI_LONG, MPI_SUM, *mpi );
      for(size_t t=0;t<species->size();t++) ldbg << "Type "<<species->at(t).name()<<" has "<<type_count[t]<<" atoms"<<std::endl;

      // compute global energy correction term
      bool all_lj = true;
      bool all_exp6 = true;
      std::vector<double> energy_correction( species->size() , 0.0 );
      std::vector<Mat3d> virial_correction( species->size() , Mat3d{} );

      for(const auto& potelem : potentials_for_pairs->m_potentials)
      {
        const auto & pot = potelem.m_params;
        all_lj = all_lj && pot.is_lj();
        all_exp6 = all_exp6 && pot.is_exp6();
      }
      
      for(const auto& potelem : potentials_for_pairs->m_potentials)
      {
        const auto & pot = potelem.m_params;
        if( str2type(potelem.m_type_a)==-1 ) { fatal_error()<<"unknown type "<<potelem.m_type_a<<" in potential description"<<std::endl; }
        if( str2type(potelem.m_type_b)==-1 ) { fatal_error()<<"unknown type "<<potelem.m_type_b<<" in potential description"<<std::endl; }
        const unsigned int ta = str2type(potelem.m_type_a);
        const unsigned int tb = str2type(potelem.m_type_b); 
        const size_t na = type_count[ta]; // number of atom of type A
        const size_t nb = type_count[tb]; // number of atom of type B
        if( all_lj )
        {
          const double rcut = pot.m_rcut;
          const double epsilon = pot.m_C_EPSILON;
          const double sigma = pot.m_D_SIGMA;
          energy_correction[ta] += 0.0; // so something here ...
          virial_correction[ta] += Mat3d{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // so something here ...
        }
        else if( all_exp6 )
        {
          const double rc = pot.m_rf.rc;
          const double A = pot.m_A;
          const double B = pot.m_B_ISEXP6;
          const double C = pot.m_C_EPSILON;
          const double D = pot.m_D_SIGMA;
          energy_correction[ta] += 0.0; // so something here ...
          virial_correction[ta] += Mat3d{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // so something here ...
        }
        else fatal_error() << "pair potentials must be all LJ or all Exp6" << std::endl;
      }
      /************************************************************/        

     
      auto& intramol_param_map = molecule_compute_parameters->m_intramol_param_map; // types to functional parameters index
      intramol_param_map.clear();      
      std::map< MoleculeGenericFuncParam , int > intramol_param_id_map ; // functional parameters to its index

      unsigned int parameter_id = 0;
      
      if( potentials_for_bonds.has_value() )
      for(const auto& bond : potentials_for_bonds->m_bond_desc)
      {
        int a = str2type( bond.species[0] );
        int b = str2type( bond.species[1] );
        if( a > b ) std::swap( a , b );        
        onika::oarray_t<int,4> types = {a,b,-1,-1};
        ldbg << "read bond "<<bond.species[0]<<","<<bond.species[1]<<" -> "<<types[0]<<","<<types[1]<<std::endl;
        auto param = bond.potential->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"bond: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key ("<<types[0]<<","<<types[1]<<","<<types[2]<<","<<types[3]<<")"<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
      }
      
      if( potentials_for_angles.has_value() )
      for(const auto& angle : potentials_for_angles->m_potentials)
      {
        int a = str2type( angle.species[0] );
        int b = str2type( angle.species[1] );
        int c = str2type( angle.species[2] );
        if( a > c ) std::swap( a , c );
        onika::oarray_t<int,4> types = {a,b,c,-1};
        ldbg << "read angle "<<angle.species[0]<<","<<angle.species[1]<<","<<angle.species[2]<<" -> "<<types[0]<<","<<types[1]<<","<<types[2]<<std::endl;
        auto param = angle.m_potential_function->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"angle: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key ("<<types[0]<<","<<types[1]<<","<<types[2]<<","<<types[3]<<")"<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
      }
      
      if( potentials_for_torsions.has_value() )
      for(const auto& torsion : potentials_for_torsions->m_potentials)
      {
        int a = str2type( torsion.species[0] );
        int b = str2type( torsion.species[1] );
        int c = str2type( torsion.species[2] );
        int d = str2type( torsion.species[3] );
        if( a > d ) { std::swap( a , d ); std::swap( b , c ); }
        onika::oarray_t<int,4> types = {a,b,c,d};
        ldbg << "read torsion "<<torsion.species[0]<<","<<torsion.species[1]<<","<<torsion.species[2]<<","<<torsion.species[3]<<" -> "<<types[0]<<","<<types[1]<<","<<types[2]<<","<<types[3]<<std::endl;
        auto param = torsion.m_potential_function->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"torsion: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key ("<<types[0]<<","<<types[1]<<","<<types[2]<<","<<types[3]<<")"<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
      }
      
      if( potentials_for_impropers.has_value() )
      for(const auto& improper : potentials_for_impropers->m_potentials)
      {
        int a = str2type( improper.species[0] );
        int b = str2type( improper.species[1] );
        int c = str2type( improper.species[2] );
        int d = str2type( improper.species[3] );
        if( b > c ) { std::swap( b , c ); }
        if( c > d ) { std::swap( c , d ); }
        if( b > c ) { std::swap( b , c ); }
        onika::oarray_t<int,4> types = {a,b,c,d};
        ldbg << "read improper "<<improper.species[0]<<","<<improper.species[1]<<","<<improper.species[2]<<","<<improper.species[3]<<" -> "<<types[0]<<","<<types[1]<<","<<types[2]<<","<<types[3]<<std::endl;
        auto param = improper.m_potential_function->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"improper: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key ("<<types[0]<<","<<types[1]<<","<<types[2]<<","<<types[3]<<")"<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
      }

      unsigned int nbparams = intramol_param_id_map.size();
      ldbg << nbparams << " different parameter sets"<<std::endl;
      molecule_compute_parameters->m_func_params.resize( nbparams );
      for(const auto &p : intramol_param_id_map)
      {
        unsigned int pidx = p.second;
        assert( pidx < nbparams );
        molecule_compute_parameters->m_func_params[ pidx ] = p.first;
      }
            
      // look for and initialize pair potentials
      if( mol_pair_weights.has_value() )
      {
        ldbg << "Intramolecular pair weighting map :"<<std::endl;      
        for( const auto & p : mol_pair_weights->m_molecule_weight )
        {
          ldbg << p.first << " : bond="<<p.second.m_bond_weight<<" , bond_rf="<<p.second.m_rf_bond_weight
                          <<" , bend="<<p.second.m_bend_weight<<" , bend_rf="<<p.second.m_rf_bend_weight
                          <<" , torsion="<<p.second.m_torsion_weight<<" , torsion_rf="<<p.second.m_rf_torsion_weight<<std::endl;
        }
      }

      for(unsigned int m=0;m<nmol;m++)
      {
        if( mol_pair_weights.has_value() )
        {
          if( mol_pair_weights->m_molecule_weight.find( molecules->m_molecules.at(m).name() ) == mol_pair_weights->m_molecule_weight.end() )
          {
            fatal_error() << "mol_pair_weights has no entry for molecule name '"<<molecules->m_molecules.at(m).name()<<"'" << std::endl;
          }
        }
      }
      
      std::map< std::pair<int,int> , LJExp6RFMultiParmsPair > pair_param_map;
      for( const auto & pp : potentials_for_pairs->m_potentials )
      {
        ldbg << "molecule pair : " << pp.m_type_a <<" / "<<pp.m_type_b<<std::endl;
        int ta = str2type(pp.m_type_a);
        int tb = str2type(pp.m_type_b);
        if( ta > tb ) std::swap(ta,tb);
        pair_param_map[ {ta,tb} ] = pp;
      }
      
      ldbg << molecule_compute_parameters->m_pair_params.size()<<" pair parameters and "<<molecule_compute_parameters->m_func_params.size()<<" intra parameters"<<std::endl;     
    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "intramolecular_setup", make_grid_variant_operator< IntramolecularSetup > );
  }

}

