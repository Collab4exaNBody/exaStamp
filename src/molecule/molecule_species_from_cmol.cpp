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

#include <exaStamp/molecule/id_map.h>
#include <exanb/core/particle_id_codec.h>

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
  class MoleculeSpeciesFromCMol : public OperatorNode
  {    
    ADD_SLOT( GridT    , grid   , INPUT_OUTPUT);

    ADD_SLOT( ParticleSpecies       , species           , INPUT , REQUIRED );
    ADD_SLOT( ParticleTypeMap       , particle_type_map , INPUT , REQUIRED );
    
    ADD_SLOT( IdMap                 , id_map            , INPUT , OPTIONAL );
    ADD_SLOT( MoleculeSpeciesVector , molecules         , INPUT_OUTPUT , REQUIRED , DocString{"Molecule descriptions"} );

  public:
    inline bool is_sink() const override final { return true; }
    
    inline void execute ()  override final
    {
      static constexpr uint64_t null_idmol = std::numeric_limits<uint64_t>::max();
    
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
      auto field_id = grid->field_accessor( field::id );
      auto field_idmol = grid->field_accessor( field::idmol );
      auto field_type = grid->field_const_accessor( field::type );  

      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        size_t n_particles = cells[cell_i].size();
        for(size_t p_i=0;p_i<n_particles;p_i++)
        {
          cells[cell_i][field_idmol][p_i] = null_idmol;
        }
      }
      
      MoleculeParser<decltype(cells),decltype(field_cmol),decltype(field_type)> molecule_parser = { cells, field_cmol, field_type, *id_map };

      std::map< std::vector<AtomNode> , int > molecule_graphs;
      std::vector< std::pair<int,int> > atom_mol_place;
      atom_mol_place.assign( grid->number_of_particles() , {-1,-1} );
              
      std::vector<AtomNode> molecule_atoms;
      std::vector<AtomNode> moldesc;
      
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
              moldesc = molecule_atoms;
              auto c0 = molecule_atoms[0].id[0];
              for(auto & an : molecule_atoms)
              {
                for(int i=0;i<5;i++) if(an.id[i]!=-1) an.id[i]-=c0;
              }
              int moltype = molecule_graphs.size();
              auto it = molecule_graphs.find( molecule_atoms );
              if( it == molecule_graphs.end() )
              {
                molecule_graphs.insert( { molecule_atoms , moltype } );
              }
              else
              {
                moltype = it->second;
              }

              for(size_t ma=0;ma<moldesc.size();ma++)
              {
                size_t cell_j=0, p_j=0;
                decode_cell_particle( moldesc[ma].cell_particle , cell_j, p_j );
                assert( moltype < molecule_graphs.size() );
                assert( cells[cell_j][field_idmol][p_j] == null_idmol );
                cells[cell_j][field_idmol][p_j] = make_molecule_id( mol_count , ma , moltype );
              }
            }
            ++ mol_count;
          }
        }
      }

      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        if( ! grid->is_ghost_cell(cell_i) )
        {
          size_t n_particles = cells[cell_i].size();
          for(size_t p_i=0;p_i<n_particles;p_i++)
          {
            assert( cells[cell_i][field_idmol][p_i] != null_idmol );
          }
        }
      }
      
      ldbg << "found "<< molecule_graphs.size() << " molecules"<<std::endl;
      for(const auto & moldesc : molecule_graphs)
      {
        ldbg<<"\tmoltype="<<moldesc.second<<" : ";
        for(const auto& an : moldesc.first)
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
        int n_atoms = moldesc.first.size();
        molecules->m_molecules[m].m_nb_atoms = n_atoms;
        for(int a=0;a<n_atoms;a++)
        {
          molecules->m_molecules[m].m_atom_type[a] = moldesc.first[a].type;
          for(int k=0;k<4;k++) molecules->m_molecules[m].m_atom_connectivity[a][k] = moldesc.first[a].id[k+1];
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
      
    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "molecule_species_from_cmol", make_grid_variant_operator< MoleculeSpeciesFromCMol > );
  }

}

