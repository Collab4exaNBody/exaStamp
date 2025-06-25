#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <exanb/core/profiling_tools.h>
#include <exaStamp/molecule/id_map.h>
#include <exanb/core/particle_id_codec.h>

#include <exaStamp/molecule/mol_connectivity.h>
#include <exaStamp/molecule/molecule_compute_param.h>
#include <exaStamp/molecule/intramolecular_pair_weight.h>

#include <mpi.h>

namespace exaStamp
{

  class IntramolecularPairPotentialList : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies             , species       , INPUT , REQUIRED );
    ADD_SLOT( IdMap                       , id_map        , INPUT  );
    ADD_SLOT( IdMapGhosts                 , id_map_ghosts , INPUT  );
    ADD_SLOT( IntramolecularPairWeighting , weight        , INPUT , IntramolecularPairWeighting{} );
    ADD_SLOT( MoleculeComputeParameterSet      , molecule_compute_parameters , INPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( IntramolecularParameterIndexLists, intramolecular_parameters , INPUT_OUTPUT, DocString{"Intramolecular functional parmater index lists"} );

    ADD_SLOT( ChemicalBonds      , chemical_bonds        , INPUT_OUTPUT );
    ADD_SLOT( ChemicalAngles     , chemical_angles       , INPUT_OUTPUT );    
    ADD_SLOT( ChemicalTorsions   , chemical_torsions     , INPUT_OUTPUT );
    ADD_SLOT( ChemicalPairs      , chemical_pairs        , INPUT_OUTPUT );
    ADD_SLOT( ChemicalPairPotMap , chemical_pair_pot_map , INPUT_OUTPUT , DocString{"Map of intramolecular pair potential parameter set and associated weights"} );

  public:
    inline void execute () override final
    {
      const size_t n_bonds = chemical_bonds->size();
      const size_t n_bends = chemical_angles->size();
      const size_t n_torsions = chemical_torsions->size();

      // index chemical links by their extrema to rapidly find if an atom pair belongs to a chemical link
      // ChemicalPairPotMap chemical_pair_pot_map;

      ProfilingTimer T;
      profiling_timer_start( T );

      // first build list of interaction pair to be excluded form standard neighbor list
      // and for which a pair weight is applied
      // this map is built using original ids, not localized ones transformed in the second parallel section
#     pragma omp parallel
      {
        // parallel construction of chemical link map : extrema pair -> link type (0 bond, 1 angle, 2 torsion)
#       pragma omp for schedule(static)
        for(size_t i=0;i<n_torsions;i++)
        {
          const auto& torsion = (*chemical_torsions)[i];
          chemical_pair_pot_map->insert( { { torsion[0] , torsion[3] } , 2 } );
        }
        
#       pragma omp for schedule(static)
        for(size_t i=0;i<n_bends;i++)
        {
          const auto & angle = (*chemical_angles)[i];
          chemical_pair_pot_map->insert( { { angle[0] , angle[2] } , 1 } );
        }
        
#       pragma omp for schedule(static)
        for(size_t i=0;i<n_bonds;i++)
        {
          const auto & bond = (*chemical_bonds)[i];
          chemical_pair_pot_map->insert( { { bond[0] , bond[1] } , 0 } );
        }
      }

      ldbg << "map build time = " << profiling_timer_elapsed_restart(T) << std::endl;

      const size_t n_chemical_pairs = chemical_pair_pot_map->size();
      chemical_pairs->assign( n_chemical_pairs , {0,0} );
      intramolecular_parameters->m_pair_param_idx.assign( n_chemical_pairs , -1 );

      const int n_type_pairs = unique_pair_count( species->size() );
      size_t pair_idx = 0;
      for(auto & m:chemical_pair_pot_map->m_meta_bucket)
      {
        for(auto & x:m)
        {
          auto pair = x.first;
          const int chemicalLinkType = x.second;
          auto & chemical_pair = (*chemical_pairs)[pair_idx];
          chemical_pair[0] = atom_from_idmap( pair[0] , *id_map , *id_map_ghosts );
          chemical_pair[1] = atom_from_idmap( pair[1] , *id_map , *id_map_ghosts );
          size_t c,p; // unused
          unsigned int ta=0, tb=0;
          decode_cell_particle( chemical_pair[0], c, p, ta );
          decode_cell_particle( chemical_pair[1], c, p, tb );
          const int pair_id = unique_pair_id(ta,tb);
          const auto & potparm = molecule_compute_parameters->m_pair_params[ n_type_pairs * chemicalLinkType + pair_id ];
          int param_idx = -1;
          if( potparm.m_pair_weight != 0.0f || potparm.m_rf_weight != 0.0f )
          {
             param_idx = n_type_pairs * chemicalLinkType + pair_id ;
          }
          intramolecular_parameters->m_pair_param_idx[ pair_idx ] = param_idx;
          //x.second = param_idx;
          ++ pair_idx;
        }
      }
      assert( pair_idx == n_chemical_pairs );

      ldbg << "pairs build time = " << profiling_timer_elapsed_restart(T) << std::endl;
      ldbg << "n_chemical_pairs = "<<n_chemical_pairs<<std::endl;
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(intramolecular_pair_list)
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "intramolecular_pair_list", make_simple_operator< IntramolecularPairPotentialList > );
  }

}
