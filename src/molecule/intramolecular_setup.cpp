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

#include <exaStamp/molecule/pair_potential_parameters.h>

#include <onika/oarray_stream.h>

#include <unordered_set>
#include <vector>

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
    using WeightMap = std::map<std::string, std::vector<double> >;

    ADD_SLOT( GridT    , grid   , INPUT_OUTPUT);

    ADD_SLOT( ParticleSpecies       , species           , INPUT , REQUIRED );
    ADD_SLOT( ParticleTypeMap       , particle_type_map , INPUT , REQUIRED );
    
    ADD_SLOT( IdMap                 , id_map            , INPUT , OPTIONAL );
    ADD_SLOT( MoleculeSpeciesVector , molecules         , INPUT_OUTPUT , REQUIRED , DocString{"Molecule descriptions"} );
    ADD_SLOT( bool                  , rebuild_molecules , INPUT , true , DocString{"Rebuild molecule description from molecular connectivity"} );    

    ADD_SLOT( BondsPotentialParameters     , potentials_for_bonds     , INPUT, OPTIONAL );
    ADD_SLOT( BendsPotentialParameters     , potentials_for_angles    , INPUT, OPTIONAL );
    ADD_SLOT( TorsionsPotentialParameters  , potentials_for_torsions  , INPUT, OPTIONAL );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT, OPTIONAL );

    ADD_SLOT( IntramolecularPairUserParamVector, potentials_for_pairs , INPUT , OPTIONAL );
    ADD_SLOT( WeightMap          , weight            , INPUT , OPTIONAL );

    ADD_SLOT( MoleculeSetComputeParams     , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );

  public:
    inline bool is_sink() const override final { return true; }
    
    inline void execute ()  override final
    {
      static constexpr MoleculeGenericFuncParam null_param = {0.0,0.0,0.0,0.0};


      /* rebuild molecule species from connectivity here */
      if( *rebuild_molecules )
      {
        if( ! grid->has_allocated_field(field::cmol) )
        {
          fatal_error() << "Impossible to rebuild molecular information without atom connectivity" << std::endl;
        }

        ldbg << "rebuild molecule species from cmol and idmol" << std::endl;
        if( ! id_map.has_value() )
        {
          fatal_error() << "id_map is needed to rebuild molecule species" << std::endl;
        }
        
        
        const size_t n_cells = grid->number_of_cells();
        auto cells = grid->cells_accessor();
        auto field_cmol = grid->field_const_accessor( field::cmol );
        auto field_id = grid->field_accessor( field::idmol );
        auto field_type = grid->field_const_accessor( field::type );  

        for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          size_t n_particles = cells[cell_i].size();
          for(size_t p_i=0;p_i<n_particles;p_i++)
          {
            cells[cell_i][field_id][p_i] = std::numeric_limits<uint64_t>::max();
          }
        }
        
        MoleculeParser<decltype(cells),decltype(field_cmol),decltype(field_type)> molecule_parser = { cells, field_cmol, field_type, *id_map };

        std::set< MoleculeAtoms > molecule_graphs;
        std::vector< std::pair<int,int> > atom_mol_place;
        atom_mol_place.assign( grid->number_of_particles() , {-1,-1} );
                
        std::vector<AtomNode> molecule_atoms;
        MoleculeAtoms moldesc;
        size_t mol_count = 0;
        for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          if( ! grid->is_ghost_cell(cell_i) )
          {
            size_t n_particles = cells[cell_i].size();
            for(size_t p_i=0;p_i<n_particles;p_i++)
            {
              molecule_atoms.clear();
              molecule_parser.explore_molecule( molecule_atoms , cells[cell_i][field_id][p_i] );
              if( ! molecule_atoms.empty() )
              {
                std::sort( molecule_atoms.begin() , molecule_atoms.end() );
                auto c0 = molecule_atoms[0].id[0];
                for(auto & an : molecule_atoms)
                {
                  for(int i=0;i<5;i++) if(an.id[i]!=-1) an.id[i]-=c0;
                }
                moldesc = MoleculeAtoms{ std::move(molecule_atoms) , -1 };
                auto it = molecule_graphs.find(moldesc);
                int moltype = molecule_graphs.size();
                if( it == molecule_graphs.end() )
                {
                  moldesc.moltype = moltype;
                  molecule_graphs.insert( moldesc );
                }
                else
                {
                  moltype = it->moltype;
                  moldesc.moltype = moltype;
                }
                for(size_t ma=0;ma<moldesc.atoms.size();ma++)
                {
                  size_t cell_j=0, p_j=0;
                  decode_cell_particle( moldesc.atoms[ma].cell_particle , cell_j, p_j );
                  assert( moltype < molecule_graphs.size() );
                  cells[cell_j][field_id][p_j] = make_molecule_id( mol_count , ma , moltype );
                }
              }
              ++ mol_count;
            }
          }
        }
        
        ldbg << "found "<< molecule_graphs.size() << " molecules"<<std::endl;
        for(const auto & moldesc : molecule_graphs)
        {
          ldbg<<"\t";
          for(auto an:moldesc.atoms)
          {
            ldbg <<" "<<species->at(an.type).name()<<":";
            for(int i=0;i<5;i++) if(an.id[i]!=-1) ldbg << ( (i==0)?"":"," ) <<an.id[i];
          }
          ldbg << std::endl;
        }
        
        molecules->m_bridge_molecules.clear();
        molecules->m_molecules.resize( molecule_graphs.size() );
        int m = 0;
        for(const auto & moldesc : molecule_graphs)
        {
          int n_atoms = moldesc.atoms.size();
          molecules->m_molecules[m].m_nb_atoms = n_atoms;
          for(int a=0;a<n_atoms;a++)
          {
            molecules->m_molecules[m].m_atom_type[a] = moldesc.atoms[a].type;
            for(int k=0;k<4;k++) molecules->m_molecules[m].m_atom_connectivity[a][k] = moldesc.atoms[a].id[k+1];
          }
          if( molecules->m_molecules[m].name().empty() )
          {
            std::ostringstream oss;
            oss << "mol_"<<m;
            molecules->m_molecules[m].set_name( oss.str() );
          }
          molecules->m_molecules[m].update_connectivity();
          molecules->m_molecules[m].print( ldbg , *species );
          ++ m;
        }
        
        //fatal_error() << "not fully implemented" << std::endl;
      }

      const unsigned int nmol = molecules->m_molecules.size();    
      molecule_compute_parameters->m_molecules.assign( nmol , MoleculeComputeParams{} );
      
      const auto & tmap = *particle_type_map;
      auto str2type = [&tmap]( const std::string& s ) -> int
        {
          auto it=tmap.find(s);
          if(it!=tmap.end()) return it->second;
          else return -1;
        };
      std::map< onika::oarray_t<int,4> , int > parameter_map; // types to functional parameters index
      std::map< MoleculeGenericFuncParam , int > parameter_id_map; // functional parameters to its index
      unsigned int parameter_id = 0;
      
      if( potentials_for_bonds.has_value() )
      for(const auto& bond : potentials_for_bonds->m_bond_desc)
      {
        ldbg << "read bond "<<bond.species[0]<<","<<bond.species[1]<<std::endl;
        int a = str2type( bond.species[0] );
        int b = str2type( bond.species[1] );
        if( a > b ) std::swap( a , b );
        
        onika::oarray_t<int,4> types = {a,b,-1,-1};
        auto param = bond.potential->generic_parameters();
        if( param != null_param )
        {
          if( parameter_id_map.find(param)==parameter_id_map.end() )
          {
            ldbg<<"add parameter pack "<<format_array(param)<<" -> "<<parameter_id<<" , for types "<<format_array(types)<<std::endl;
            parameter_id_map[param] = parameter_id++;
          }
          parameter_map[ types ] = parameter_id_map[param];
        }
      }
      
      if( potentials_for_angles.has_value() )
      for(const auto& angle : potentials_for_angles->m_potentials)
      {
        ldbg << "read angle "<<angle.species[0]<<","<<angle.species[1]<<","<<angle.species[2]<<std::endl;
        int a = str2type( angle.species[0] );
        int b = str2type( angle.species[1] );
        int c = str2type( angle.species[2] );
        if( a > c ) std::swap( a , c );
        
        onika::oarray_t<int,4> types = {a,b,c,-1};
        auto param = angle.m_potential_function->generic_parameters();
        if( param != null_param )
        {
          if( parameter_id_map.find(param)==parameter_id_map.end() )
          {
            ldbg<<"add parameter pack "<<format_array(param)<<" -> "<<parameter_id<<" , for types "<<format_array(types)<<std::endl;
            parameter_id_map[param] = parameter_id++;
          }
          parameter_map[ types ] = parameter_id_map[param];
        }
      }
      
      if( potentials_for_torsions.has_value() )
      for(const auto& torsion : potentials_for_torsions->m_potentials)
      {
        ldbg << "read torsion "<<torsion.species[0]<<","<<torsion.species[1]<<","<<torsion.species[2]<<","<<torsion.species[3]<<std::endl;
        int a = str2type( torsion.species[0] );
        int b = str2type( torsion.species[1] );
        int c = str2type( torsion.species[2] );
        int d = str2type( torsion.species[3] );
        if( a > d ) { std::swap( a , d ); std::swap( b , c ); }
                
        onika::oarray_t<int,4> types = {a,b,c,d};
        auto param = torsion.m_potential_function->generic_parameters();
        if( param != null_param )
        {
          if( parameter_id_map.find(param)==parameter_id_map.end() )
          {
            ldbg<<"add parameter pack "<<format_array(param)<<" -> "<<parameter_id<<" , for types "<<format_array(types)<<std::endl;
            parameter_id_map[param] = parameter_id++;
          }
          parameter_map[ types ] = parameter_id_map[param];
        }
      }
      
      if( potentials_for_impropers.has_value() )
      for(const auto& improper : potentials_for_impropers->m_potentials)
      {
        ldbg << "read improper "<<improper.species[0]<<","<<improper.species[1]<<","<<improper.species[2]<<","<<improper.species[3]<<std::endl;
        int a = str2type( improper.species[0] );
        int b = str2type( improper.species[1] );
        int c = str2type( improper.species[2] );
        int d = str2type( improper.species[3] );
        if( b > c ) { std::swap( b , c ); }
        if( c > d ) { std::swap( c , d ); }
        if( b > c ) { std::swap( b , c ); }

        onika::oarray_t<int,4> types = {a,b,c,d};
        auto param = improper.m_potential_function->generic_parameters();
        if( param != null_param )
        {
          if( parameter_id_map.find(param)==parameter_id_map.end() )
          {
            ldbg<<"add parameter pack "<<format_array(param)<<" -> "<<parameter_id<<" , for types "<<format_array(types)<<std::endl;
            parameter_id_map[param] = parameter_id++;
          }
          parameter_map[ types ] = parameter_id_map[param];
        }
      }

      unsigned int nbparams = parameter_id_map.size();
      ldbg << nbparams << " different parameter sets"<<std::endl;
      molecule_compute_parameters->m_func_params.resize( nbparams );
      for(const auto &p : parameter_id_map)
      {
        unsigned int pidx = p.second;
        assert( pidx<nbparams);
        molecule_compute_parameters->m_func_params[ pidx ] = p.first;
      }
      
      for(unsigned int i=0;i<nmol;i++)
      {
      
        for(unsigned int j=0;j<molecules->m_molecules.at(i).m_nb_bonds;j++)
        {
          uint32_t a = molecules->m_molecules.at(i).m_bonds[j][0];
          uint32_t b = molecules->m_molecules.at(i).m_bonds[j][1];
          int ta = molecules->m_molecules.at(i).m_atom_type[a];
          int tb = molecules->m_molecules.at(i).m_atom_type[b];
          if( ta > tb ) std::swap(ta,tb);
          onika::oarray_t<int,4> types = {ta,tb,-1,-1};
          auto types_it = parameter_map.find( types );
          if( types_it != parameter_map.end() )
          {
            uint32_t pidx = types_it->second;
            assert( pidx>=0 && pidx<nbparams );
            molecule_compute_parameters->m_molecules[i].m_bonds[ molecule_compute_parameters->m_molecules[i].m_nb_bonds ++ ] = (a<<16) | (b<<8) | pidx ;
            ldbg<<"molecule"<<i<<".bond"<<j<<" : A{"<<a<<","<<b<<"} / T{"<<ta<<","<<tb<<"} -> "<<pidx<<std::endl;
          }
          else
          {
            ldbg<<"molecule"<<i<<".bond"<<j<<" : A{"<<a<<","<<b<<"} ignored"<<std::endl;
          }
        }
        
        for(unsigned int j=0;j<molecules->m_molecules.at(i).m_nb_bends;j++)
        {
          uint32_t a = molecules->m_molecules.at(i).m_bends[j][0];
          uint32_t b = molecules->m_molecules.at(i).m_bends[j][1];
          uint32_t c = molecules->m_molecules.at(i).m_bends[j][2];
          int ta = molecules->m_molecules.at(i).m_atom_type[a];
          int tb = molecules->m_molecules.at(i).m_atom_type[b];
          int tc = molecules->m_molecules.at(i).m_atom_type[c];
          if( ta > tc ) std::swap(ta,tc);
          onika::oarray_t<int,4> types = {ta,tb,tc,-1};
          auto types_it = parameter_map.find( types );
          if( types_it != parameter_map.end() )
          {
            uint32_t pidx = types_it->second;
            assert( pidx>=0 && pidx<nbparams );
            molecule_compute_parameters->m_molecules[i].m_bends[ molecule_compute_parameters->m_molecules[i].m_nb_bends ++ ] = (a<<24) | (b<<16) | (c<<8) | pidx ;
            ldbg<<"molecule"<<i<<".angle"<<j<<" : A{"<<a<<","<<b<<","<<c<<"} / T{"<<ta<<","<<tb<<","<<tc<<"} -> "<<pidx<<std::endl;
          }
          else
          {
            ldbg<<"molecule"<<i<<".angle"<<j<<" : A{"<<a<<","<<b<<","<<c<<"} ignored"<<std::endl;
          }

        }

        for(unsigned int j=0;j<molecules->m_molecules.at(i).m_nb_torsions;j++)
        {
          uint64_t a = molecules->m_molecules.at(i).m_torsions[j][0];
          uint64_t b = molecules->m_molecules.at(i).m_torsions[j][1];
          uint64_t c = molecules->m_molecules.at(i).m_torsions[j][2];
          uint64_t d = molecules->m_molecules.at(i).m_torsions[j][3];
          int ta = molecules->m_molecules.at(i).m_atom_type[a];
          int tb = molecules->m_molecules.at(i).m_atom_type[b];
          int tc = molecules->m_molecules.at(i).m_atom_type[c];
          int td = molecules->m_molecules.at(i).m_atom_type[d];
          if( ta > td ) { std::swap(ta,td); std::swap(tb,tc); }
          onika::oarray_t<int,4> types = {ta,tb,tc,td};
          auto types_it = parameter_map.find( types );
          if( types_it != parameter_map.end() )
          {
            uint64_t pidx = types_it->second;
            assert( pidx>=0 && pidx<nbparams );
            molecule_compute_parameters->m_molecules[i].m_torsions[ molecule_compute_parameters->m_molecules[i].m_nb_torsions ++ ] = (a<<32) | (b<<24) | (c<<16) | (d<<8) | pidx ;
            ldbg<<"molecule"<<i<<".torsion"<<j<<" : A{"<<a<<","<<b<<","<<c<<","<<d<<"} / T{"<<ta<<","<<tb<<","<<tc<<","<<td<<"} -> "<<pidx<<std::endl;
          }
          else
          {
            ldbg<<"molecule"<<i<<".torsion"<<j<<" : A{"<<a<<","<<b<<","<<c<<","<<d<<"} ignored"<<std::endl;
          }
        }
        
        for(unsigned int j=0;j<molecules->m_molecules.at(i).m_nb_impropers;j++)
        {
          uint64_t a = molecules->m_molecules.at(i).m_impropers[j][0];
          uint64_t b = molecules->m_molecules.at(i).m_impropers[j][1];
          uint64_t c = molecules->m_molecules.at(i).m_impropers[j][2];
          uint64_t d = molecules->m_molecules.at(i).m_impropers[j][3];
          int ta = molecules->m_molecules.at(i).m_atom_type[a];
          int tb = molecules->m_molecules.at(i).m_atom_type[b];
          int tc = molecules->m_molecules.at(i).m_atom_type[c];
          int td = molecules->m_molecules.at(i).m_atom_type[d];
          if( tb > tc ) { std::swap(tb,tc); }
          if( tc > td ) { std::swap(tc,td); }
          if( tb > tc ) { std::swap(tb,tc); }
          onika::oarray_t<int,4> types = {ta,tb,tc,td};
          auto types_it = parameter_map.find( types );
          if( types_it != parameter_map.end() )
          {
            uint64_t pidx = types_it->second;
            assert( pidx>=0 && pidx<nbparams );
            molecule_compute_parameters->m_molecules[i].m_impropers[ molecule_compute_parameters->m_molecules[i].m_nb_impropers ++ ] = (a<<32) | (b<<24) | (c<<16) | (d<<8) | pidx ;
            ldbg<<"molecule"<<i<<".improper"<<j<<" : A{"<<a<<","<<b<<","<<c<<","<<d<<"} / T{"<<ta<<","<<tb<<","<<tc<<","<<td<<"} -> "<<pidx<<std::endl;
          }
          else
          {
            ldbg<<"molecule"<<i<<".improper"<<j<<" : A{"<<a<<","<<b<<","<<c<<","<<d<<"} ignored"<<std::endl;
          }
        }

      }
      
      
      // look for and initialize pair potentials
      if( weight.has_value() )
      {
        for(const auto& pw:*weight)
        {
          ldbg <<pw.first<<" pair weighting : bond="<<pw.second[0]<<" angle="<<pw.second[1]<<" torsion="<<pw.second[2]<<std::endl;
        }
      }


      std::map< std::pair<int,int> , IntramolecularPairUserParam > pair_param_map;
      for( const auto & pp : *potentials_for_pairs )
      {
        ldbg << "molecule pair : " << pp.m_type_a <<" / "<<pp.m_type_b<<std::endl;
        int ta = str2type(pp.m_type_a);
        int tb = str2type(pp.m_type_b);
        if( ta > tb ) std::swap(ta,tb);
        pair_param_map[ {ta,tb} ] = pp;
      }

      // build pairs with weighting
      if( potentials_for_pairs.has_value() )
      {
        int nb_rf_params = 0;
        int nb_ljexp6_params = 0;
        std::map< IntramolecularRFParam , int > pair_rf_id_map;
        std::map< IntramolecularLJExp6Param , int > pair_ljexp6_id_map;

        for(unsigned int m=0;m<nmol;m++)
        {          
          std::set< std::array<unsigned int,2> > bonds;
          std::set< std::array<unsigned int,2> > bends;
          std::set< std::array<unsigned int,2> > torsions;
          
          for(unsigned int i=0;i<molecules->m_molecules.at(m).m_nb_bonds;i++)
          {
            const unsigned int a = molecules->m_molecules.at(m).m_bonds[i][0];
            const unsigned int b = molecules->m_molecules.at(m).m_bonds[i][1];
            bonds.insert( {a,b} );
          }
          for(unsigned int i=0;i<molecules->m_molecules.at(m).m_nb_bends;i++)
          {
            const unsigned int a = molecules->m_molecules.at(m).m_bends[i][0];
            const unsigned int b = molecules->m_molecules.at(m).m_bends[i][2];
            bends.insert( {a,b} );
          }
          for(unsigned int i=0;i<molecules->m_molecules.at(m).m_nb_torsions;i++)
          {
            const unsigned int a = molecules->m_molecules.at(m).m_torsions[i][0];
            const unsigned int b = molecules->m_molecules.at(m).m_torsions[i][3];
            torsions.insert( {a,b} );
          }
          
          molecule_compute_parameters->m_molecules[m].m_nb_pairs = 0;
          for(unsigned int i=0;i<molecules->m_molecules.at(m).m_nb_atoms;i++)
          {
            int type_a = molecules->m_molecules.at(m).m_atom_type[i];
            for(unsigned int j=i+1;j<molecules->m_molecules.at(m).m_nb_atoms;j++)
            {
              int type_b = molecules->m_molecules.at(m).m_atom_type[j];
              int ta=type_a , tb=type_b;
              if( ta > tb ) std::swap(ta,tb);
              auto pair_param_it = pair_param_map.find( {ta,tb} );
              if( pair_param_it != pair_param_map.end() )
              {
                auto pair_param = pair_param_it->second;
                pair_param.m_rf.m_c1 = species->at(type_a).m_charge;
                pair_param.m_rf.m_c2 = species->at(type_b).m_charge;
                double w = 0.0;
                const char* type="pair";
                if( weight.has_value() )
                {
                  if( torsions.find({i,j}) != torsions.end() )
                  {
                    w = weight->at(molecules->m_molecules.at(m).name())[2];
                    type="torsion";
//                    if(w==0.0) ldbg<<"pair "<<i<<" / "<<j<<" skipped because torsion weight=0"<<std::endl;
                  }
                  else if( bends.find({i,j}) != bends.end() )
                  {
                    w = weight->at(molecules->m_molecules.at(m).name())[1];
                    type="bend";
//                    if(w==0.0) ldbg<<"pair "<<i<<" / "<<j<<" skipped because bend weight=0"<<std::endl;
                  }
                  else if( bonds.find({i,j}) != bonds.end() )
                  {
                    w = weight->at(molecules->m_molecules.at(m).name())[0];
                    type="bond";
//                    if(w==0.0) ldbg<<"pair "<<i<<" / "<<j<<" skipped because bond weight=0"<<std::endl;
                  }
                  else
                  {
                    w = 1.0;
                  }
                }
                if( w > 0.0 )
                {
                  if( molecule_compute_parameters->m_molecules[m].m_nb_pairs >= MAX_MOLECULE_PAIRS )
                  {
                    fatal_error() << "intramolecular pair overflow"<<std::endl;
                  }
                  if( w != 1.0 && w != 0.5 )
                  {
                    fatal_error() << "intramolecular pair potential weighting only supports 1 or 0.5 weight factor"<<std::endl;
                  }

                  ldbg << molecules->m_molecules.at(m).name() <<" : pair#"<<molecule_compute_parameters->m_molecules[m].m_nb_pairs<<" , "<<i<<" / "<<j<<" , type="<<type<<" , weight="<<w<<std::endl;
                  
                  if( pair_rf_id_map.find(pair_param.m_rf) == pair_rf_id_map.end() )
                  {
                    ldbg << "insert RF pair param #"<<nb_rf_params <<" for types "<<type_a<<"/"<<type_b
                         << " ,rcut="<<pair_param.m_rf.m_param.rc<<" ,c1="<<pair_param.m_rf.m_c1<<" ,c2="<<pair_param.m_rf.m_c2<<std::endl;
                    pair_rf_id_map[pair_param.m_rf] = nb_rf_params++;
                  }
                  
                  if( pair_ljexp6_id_map.find(pair_param.m_ljexp6) == pair_ljexp6_id_map.end() )
                  {
                    ldbg << "insert LJ/Exp6 pair param #"<<nb_ljexp6_params <<" for types "<<type_a<<"/"<<type_b<< " ,rcut="<<pair_param.m_ljexp6.m_rcut<<std::endl;
                    pair_ljexp6_id_map[pair_param.m_ljexp6] = nb_ljexp6_params++;
                  }

                  uint64_t pidx_rf = pair_rf_id_map[pair_param.m_rf];
                  uint64_t pidx_ljexp6 = pair_ljexp6_id_map[pair_param.m_ljexp6];

                  uint64_t wi = (w==1.0) ? 2 : 1;
                  //lout << i<<" / "<<j<<" pair_weight = "<<w<<" , types = "<< int(molecules->m_molecules.at(m).m_atom_type[i]) <<" / "<< int(molecules->m_molecules.at(m).m_atom_type[j]) <<std::endl;
                  molecule_compute_parameters->m_molecules[m].m_pairs[ molecule_compute_parameters->m_molecules[m].m_nb_pairs ] = (wi<<32) | (i<<24) | (j<<16) | (pidx_rf<<8) | pidx_ljexp6;
                  ++ molecule_compute_parameters->m_molecules[m].m_nb_pairs;
                }
              }
              else
              {
                ldbg<<"pair "<<i<<" / "<<j<<" skipped because no pair parameter for types "<<type_a<<" / "<<type_b<<std::endl;
              }
            }
          }
        }

        molecule_compute_parameters->m_rf_params.resize(nb_rf_params);
        for(const auto& pp:pair_rf_id_map)
        {
          molecule_compute_parameters->m_rf_params[pp.second] = pp.first;
        }
        molecule_compute_parameters->m_ljexp6_params.resize(nb_ljexp6_params);
        for(const auto& pp:pair_ljexp6_id_map)
        {
          molecule_compute_parameters->m_ljexp6_params[pp.second] = pp.first;
          molecule_compute_parameters->m_ljexp6_params[pp.second].update_ecut();
        }
        
      }
      
      ldbg << molecule_compute_parameters->m_rf_params.size()<<" RF parameters, "<<molecule_compute_parameters->m_ljexp6_params.size()<<" LJ/Exp6 parameters and "<<molecule_compute_parameters->m_func_params.size()<<" intra parameters"<<std::endl;     
    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "intramolecular_setup", make_grid_variant_operator< IntramolecularSetup > );
  }

}

