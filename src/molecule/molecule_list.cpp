#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/molecule/molecule_list.h>
#include <exaStamp/molecule/molecule_species.h>

#include <unordered_map>

#include <omp.h>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_idmol>
    >
  class MoleculeList : public OperatorNode
  {
    ADD_SLOT( GridT    , grid   , INPUT_OUTPUT);
    ADD_SLOT( MoleculeSpeciesVector , molecules  , INPUT, REQUIRED, DocString{"Molecule descriptions"} );
    ADD_SLOT( MoleculeLists , molecule_list , INPUT_OUTPUT );

  public:
    inline void execute ()  override final
    {
      auto cells = grid->cells();

      std::vector< std::unordered_map<uint64_t,uint64_t> > mol_id_maps;
      const int max_nt = omp_get_max_threads();
      mol_id_maps.resize( max_nt );

      //set an id of molecules
      //we follow the connectivity list
      const size_t n_cells = grid->number_of_cells();
      
#     pragma omp parallel
      {      
        const int tid = omp_get_thread_num();
        assert( omp_get_num_threads() <= max_nt );
        assert( tid < omp_get_num_threads() );
        auto & id_map = mol_id_maps[tid];
        
#       pragma omp for schedule(guided)        
        for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          const bool is_ghost = grid->is_ghost_cell( cell_i );
          const size_t n_particles = cells[cell_i].size();
          for(size_t i=0;i<n_particles;i++)
          {
            const auto idmol = cells[cell_i][field::idmol][i];
            assert( is_ghost || idmol!=std::numeric_limits<uint64_t>::max() );
            if( idmol != std::numeric_limits<uint64_t>::max() )
            {
              uint64_t mid = molecule_instance_from_id( idmol );
              uint64_t mtype = molecule_type_from_id( idmol );
              assert( mtype < molecules->m_molecules.size() );
              if( id_map.find(mid) == id_map.end() )
              {
                id_map.insert( { mid , mtype } );
              }
            }
          }
        }
      }

      // merge maps from al threads
      auto& id_map = mol_id_maps[0];
      for(int i=1;i<max_nt;i++)
      {
        id_map.merge( mol_id_maps[i] );
        mol_id_maps[i].clear();
      }
      ldbg << "local proc has " << id_map.size() << " unique molecules" << std::endl;
      for(const auto& mol : id_map)
      {
        assert( mol.second >= 0 && mol.second < molecules->m_molecules.size() );
      }

      // now compute molecule offsets
      std::vector< std::pair<uint64_t,uint64_t> > sorted_molids( id_map.begin() , id_map.end() );
      id_map.clear();
      std::sort( sorted_molids.begin() , sorted_molids.end() ,
        [](const std::pair<uint64_t,uint64_t>& a, const std::pair<uint64_t,uint64_t>& b)->bool
        {
          return a.second < b.second;
        });
        
      const size_t nmol = sorted_molids.size();
      size_t mol_data_size = 0;
      molecule_list->m_offset.resize( nmol );
      for(size_t i=0;i<nmol;i++)
      {
        const uint64_t mid = sorted_molids[i].first;
        const unsigned int mtype = sorted_molids[i].second;
        if( mtype >= molecules->m_molecules.size() )
        {
          fatal_error() << "molecule #"<<i<<" has type "<<mtype<<" but molecules size is "<< molecules->m_molecules.size() << std::endl;
        }
        id_map[ mid ] = i;
        molecule_list->m_offset[i] = mol_data_size;
        mol_data_size += 1 + molecules->m_molecules.at(mtype).m_nb_atoms;
      }
      ldbg << "molecule data size = " << mol_data_size << std::endl;
      
      // resize molecule data
      molecule_list->m_data.assign( mol_data_size , -1 );

      // fill in molecule types
      for(size_t m=0;m<nmol;m++)
      {
        molecule_list->m_data[ molecule_list->m_offset[m] ] = sorted_molids[m].second;
      }
      sorted_molids.clear();
      sorted_molids.shrink_to_fit();

      // fill molecule data (atom localizations)
#     pragma omp parallel
      {      
        
#       pragma omp for schedule(dynamic)        
        for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( grid->is_ghost_cell(cell_i) )
        {
          const size_t n_particles = cells[cell_i].size();
          for(size_t p_i=0;p_i<n_particles;p_i++)
          {
            const auto idmol = cells[cell_i][field::idmol][p_i];
            const uint64_t mid = molecule_instance_from_id( idmol );
            assert( id_map.find(mid) != id_map.end() );
            const uint64_t mplace = molecule_place_from_id( idmol );
            assert( mplace < molecules->m_molecules.at( molecule_type_from_id(idmol) ).m_nb_atoms );
            const uint64_t molidx = id_map[mid];
            assert( molecule_list->m_data[ molecule_list->m_offset[molidx] ] == molecule_type_from_id(idmol) );
            molecule_list->m_data[ molecule_list->m_offset[molidx] + 1 + mplace ] = encode_cell_particle( cell_i , p_i );
          }
        }

#       pragma omp barrier
        
#       pragma omp for schedule(dynamic)        
        for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
        {
          const size_t n_particles = cells[cell_i].size();
          for(size_t p_i=0;p_i<n_particles;p_i++)
          {
            const auto idmol = cells[cell_i][field::idmol][p_i];
            const uint64_t mid = molecule_instance_from_id( idmol );
            assert( id_map.find(mid) != id_map.end() );
            const uint64_t mplace = molecule_place_from_id( idmol );
            assert( mplace < molecules->m_molecules.at( molecule_type_from_id(idmol) ).m_nb_atoms );
            const uint64_t molidx = id_map[mid];
            assert( molecule_list->m_data[ molecule_list->m_offset[molidx] ] == molecule_type_from_id(idmol) );
            molecule_list->m_data[ molecule_list->m_offset[molidx] + 1 + mplace ] = encode_cell_particle( cell_i , p_i );
          }
        }
        
      }

#     ifndef NDEBUG
      // check if everything is ok
      for(size_t m=0;m<nmol;m++)
      {
        int mtype = molecule_list->m_data[ molecule_list->m_offset[m] ];
        int natoms = molecules->m_molecules.at(mtype).m_nb_atoms;
        for(int a=0;a<natoms;a++)
        {
          if( molecule_list->m_data[ molecule_list->m_offset[m] + 1 + a] == -1 )
          {
            ldbg << "Molecule #"<<m<<" (type="<<mtype<<") misses atom #"<<a<<std::endl;
          }
        }
      }
#     endif

    }

  };

  template<class GridT> using MoleculeListTmpl = MoleculeList<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("molecule_list", make_grid_variant_operator< MoleculeListTmpl > );
  }

}
