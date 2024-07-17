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
    ADD_SLOT( ChemicalAngles    , chemical_angles    , INPUT_OUTPUT );    
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
          b[0] = atom_from_idmap( b[0] , *id_map , *id_map_ghosts );
          b[1] = atom_from_idmap( b[1] , *id_map , *id_map_ghosts );
        }

#       pragma omp for schedule(static)
        for(size_t i=0;i<n_bends;i++)
        {
          auto& b = (*chemical_angles)[i];
          b[0] = atom_from_idmap( b[0] , *id_map , *id_map_ghosts );
          b[1] = atom_from_idmap( b[1] , *id_map , *id_map_ghosts );
          b[2] = atom_from_idmap( b[2] , *id_map , *id_map_ghosts );
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

    }
  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "mol_optimize_connectivity", make_grid_variant_operator< MolOptimizeConnectivity > );
  }

}
