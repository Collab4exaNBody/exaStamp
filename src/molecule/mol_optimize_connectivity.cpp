#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/log.h>
#include <exaStamp/molecule/id_map.h>
#include <exanb/core/thread.h>
#include <exanb/core/particle_id_codec.h>

#include <exaStamp/molecule/mol_connectivity.h>

#include <chrono>
#include <unordered_set>


#include <mpi.h>

namespace exaStamp
{

  template< class GridT >
  class MolOptimizeConnectivity : public OperatorNode
  {
    ADD_SLOT( GridT             , grid               , INPUT_OUTPUT );

    ADD_SLOT( ChemicalBonds     , chemical_bonds     , INPUT_OUTPUT );
    ADD_SLOT( GridAtomsToChemicalChain  , atoms_to_bonds     , INPUT_OUTPUT );

    ADD_SLOT( ChemicalAngles    , chemical_angles    , INPUT_OUTPUT );
#   ifdef XSTAMP_USE_BEND_FORCE_BUFFER
    ADD_SLOT( GridAtomsToChemicalChain  , atoms_to_angles     , INPUT_OUTPUT );
#   endif
    
    ADD_SLOT( ChemicalTorsions  , chemical_torsions  , INPUT_OUTPUT );
    ADD_SLOT( ChemicalImpropers , chemical_impropers , INPUT_OUTPUT );

    ADD_SLOT( IdMap             , id_map             , INPUT  );
    ADD_SLOT( IdMapGhosts       , id_map_ghosts      , INPUT  );

  public:
    inline void execute () override final
    {
      auto cells = grid->cells();
      size_t n_cells = grid->number_of_cells();
      const size_t n_bonds = chemical_bonds->size();
      const size_t n_bends = chemical_angles->size();
      const size_t n_torsions = chemical_torsions->size();
      const size_t n_impropers = chemical_impropers->size();
      
      std::vector< std::vector< std::pair<uint64_t,uint64_t> > > particle_to_bond( n_cells );
      std::vector< std::vector< std::pair<uint64_t,uint64_t> > > particle_to_angle( n_cells );

      spin_mutex_array cell_locks;
      cell_locks.resize( n_cells );

#     pragma omp parallel
      {
#       pragma omp for schedule(static)
        for(size_t i=0;i<n_bonds;i++)
        {
          auto& b = (*chemical_bonds)[i];
          
          // bond atom 0
          b[0] = atom_from_idmap( b[0] , *id_map , *id_map_ghosts );
          size_t cell1 = std::numeric_limits<size_t>::max();
          size_t particle1 = std::numeric_limits<size_t>::max();
          unsigned int type1 = std::numeric_limits<unsigned int>::max();
          decode_cell_particle(b[0], cell1, particle1, type1);
          cell_locks[cell1].lock();
          particle_to_bond[cell1].push_back( { particle1 , i*4+0 } ); // position 0 in bond i => i*4+0
          cell_locks[cell1].unlock();

          // bond atom 1
          b[1] = atom_from_idmap( b[1] , *id_map , *id_map_ghosts );
          size_t cell2 = std::numeric_limits<size_t>::max();
          size_t particle2 = std::numeric_limits<size_t>::max();
          unsigned int type2 = std::numeric_limits<unsigned int>::max();
          decode_cell_particle(b[1], cell2, particle2, type2);
          cell_locks[cell2].lock();
          particle_to_bond[cell2].push_back( { particle2 , i*4+1 } ); // position 1 in bond i => i*4+1        
          cell_locks[cell2].unlock();
        }

#       pragma omp for schedule(static)
        for(size_t i=0;i<n_bends;i++)
        {
          auto& b = (*chemical_angles)[i];
          b[0] = atom_from_idmap( b[0] , *id_map , *id_map_ghosts );
          b[1] = atom_from_idmap( b[1] , *id_map , *id_map_ghosts );
          b[2] = atom_from_idmap( b[2] , *id_map , *id_map_ghosts );

#         ifdef XSTAMP_USE_BEND_FORCE_BUFFER
          // bond atom 0
          size_t cell0 = std::numeric_limits<size_t>::max();
          size_t particle0 = std::numeric_limits<size_t>::max();
          unsigned int type0 = std::numeric_limits<unsigned int>::max();
          decode_cell_particle(b[0], cell0, particle0, type0);
          cell_locks[cell0].lock();
          particle_to_angle[cell0].push_back( { particle0 , i*4+0 } ); // position 0 in bond i => i*4+0
          cell_locks[cell0].unlock();

          // bond atom 1
          size_t cell1 = std::numeric_limits<size_t>::max();
          size_t particle1 = std::numeric_limits<size_t>::max();
          unsigned int type1 = std::numeric_limits<unsigned int>::max();
          decode_cell_particle(b[1], cell1, particle1, type1);
          cell_locks[cell1].lock();
          particle_to_angle[cell1].push_back( { particle1 , i*4+1 } ); // position 1 in bond i => i*4+1
          cell_locks[cell1].unlock();

          // bond atom 2
          size_t cell2 = std::numeric_limits<size_t>::max();
          size_t particle2 = std::numeric_limits<size_t>::max();
          unsigned int type2 = std::numeric_limits<unsigned int>::max();
          decode_cell_particle(b[2], cell2, particle2, type2);
          cell_locks[cell2].lock();
          particle_to_angle[cell2].push_back( { particle2 , i*4+2 } ); // position 2 in bond i => i*4+2
          cell_locks[cell2].unlock();
#         endif
        }

#       pragma omp for schedule(static)
        for(size_t i=0;i<n_torsions;i++)
        {
          auto& t = (*chemical_torsions)[i];
          t[0] = atom_from_idmap( t[0] , *id_map , *id_map_ghosts );
          t[1] = atom_from_idmap( t[1] , *id_map , *id_map_ghosts );
          t[2] = atom_from_idmap( t[2] , *id_map , *id_map_ghosts );
          t[3] = atom_from_idmap( t[3] , *id_map , *id_map_ghosts );
        }

#       pragma omp for schedule(static)
        for(size_t i=0;i<n_impropers;i++)
        {
          auto& t = (*chemical_impropers)[i];
          t[0] = atom_from_idmap( t[0] , *id_map , *id_map_ghosts );
          t[1] = atom_from_idmap( t[1] , *id_map , *id_map_ghosts );
          t[2] = atom_from_idmap( t[2] , *id_map , *id_map_ghosts );
          t[3] = atom_from_idmap( t[3] , *id_map , *id_map_ghosts );
        }

      }

      atoms_to_bonds->resize( n_cells );

#     ifdef XSTAMP_USE_BEND_FORCE_BUFFER
      atoms_to_angles->resize( n_cells );
#     endif

#     pragma omp parallel for schedule(dynamic)
      for(size_t i=0;i<n_cells;i++)
      {
        std::sort( particle_to_bond[i].begin() , particle_to_bond[i].end() );
        size_t n_particles = cells[i].size();
        
        auto& cell_ab = (*atoms_to_bonds)[i];
        cell_ab.clear();
        size_t k = 0;
        size_t nbonds = particle_to_bond[i].size();
        for(size_t j=0;j<n_particles;j++)
        {
#         ifndef NDEBUG
          cell_ab.push_back(j); // for debug only
#         endif
          size_t count_idx = cell_ab.size();
          cell_ab.push_back(0);
          while( k<nbonds && particle_to_bond[i][k].first==j )
          {
            ++ cell_ab[count_idx];
            cell_ab.push_back( particle_to_bond[i][k].second );
            ++ k;
          }
        }
        assert( k == nbonds );

#       ifdef XSTAMP_USE_BEND_FORCE_BUFFER
        std::sort( particle_to_angle[i].begin() , particle_to_angle[i].end() );
        auto& cell_aa = (*atoms_to_angles)[i];
        cell_aa.clear();
        k = 0;
        size_t nangles = particle_to_angle[i].size();
        for(size_t j=0;j<n_particles;j++)
        {
#         ifndef NDEBUG
          cell_aa.push_back(j); // for debug only
#         endif
          size_t count_idx = cell_aa.size();
          cell_aa.push_back(0);
          while( k<nangles && particle_to_angle[i][k].first==j )
          {
            ++ cell_aa[count_idx];
            cell_aa.push_back( particle_to_angle[i][k].second );
            ++ k;
          }
        }
        assert( k == nangles );
#       endif
      }

    }
  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "mol_optimize_connectivity", make_grid_variant_operator< MolOptimizeConnectivity > );
  }

}
