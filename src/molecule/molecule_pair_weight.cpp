#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/std_array_hash.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/compact_grid_pair_weights.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/mt_concurrent_map.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/particle_neighbors/chunk_neighbors_apply.h>

#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/intramolecular_pair_weight.h>
#include <exaStamp/molecule/molecule_compute_param.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_type, field::_id>
    >
  class MoleculePairWeightChunk : public OperatorNode
  {
    ADD_SLOT( GridT              , grid              , INPUT , REQUIRED);    
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );    

    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT , DocString{"Molecule descriptions"} );

    ADD_SLOT( IntramolecularPairWeighting , weight   , INPUT , IntramolecularPairWeighting{} );
    ADD_SLOT( MoleculeComputeParameterSet      , molecule_compute_parameters , INPUT, DocString{"Intramolecular functionals' parameters"} );

    ADD_SLOT( CompactGridPairWeights , compact_nbh_weight , INPUT_OUTPUT );
    
    ADD_SLOT( ParticleSpecies    , species           , INPUT , REQUIRED );
    ADD_SLOT( ChemicalBonds      , chemical_bonds    , INPUT , OPTIONAL );
    ADD_SLOT( ChemicalAngles     , chemical_angles   , INPUT , OPTIONAL );
    ADD_SLOT( ChemicalTorsions   , chemical_torsions , INPUT , OPTIONAL );

    ADD_SLOT( bool               , ignore_intramolecular , INPUT , false );

    ADD_SLOT( ChemicalPairPotMap , chemical_pair_pot_map , INPUT , DocString{"Map of intramolecular pair potential parameter set and associated weights"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      const bool no_intramolecular = *ignore_intramolecular;
      
      std::map< std::string , size_t > molecule_name_map;
      size_t n_molecules = molecules->m_molecules.size();
      for(size_t i=0;i<n_molecules;i++)
      {
        molecule_name_map[ molecules->m_molecules[i].name() ] = i;
      }
      std::unordered_map< uint64_t , MolecularPairWeight > molid_weight_map;

      ldbg << "Intramolecular pair weighting map :"<<std::endl;      
      for( const auto & p : weight->m_molecule_weight )
      {
        molid_weight_map[ molecule_name_map[p.first] ] = p.second;
        ldbg << p.first << " : bond="<<p.second.m_bond_weight<<" , bond_rf="<<p.second.m_rf_bond_weight
                        <<" , bend="<<p.second.m_bend_weight<<" , bend_rf="<<p.second.m_rf_bend_weight
                        <<" , torsion="<<p.second.m_torsion_weight<<" , torsion_rf="<<p.second.m_rf_torsion_weight<<std::endl;
      }

      auto cells = grid->cells_accessor();
      auto field_idmol = grid->field_accessor( field::idmol );
      size_t n_cells = grid->number_of_cells(); // nbh.size();

      CompactGridPairWeights& cgpw = *compact_nbh_weight;
      cgpw.m_cell_weights.resize( n_cells );

      const IJK dims = grid->dimension();
      int gl = grid->ghost_layers();

      size_t n_total_nbh=0, n_positive_weights=0;

      // index chemical links by their extrema to rapidly find if an atom pair belongs to a chemical link
      const auto & chemPotMap = *chemical_pair_pot_map;
      const auto & mcp = *molecule_compute_parameters;

#     pragma omp parallel
      {

#       pragma omp for
        for(size_t i=0;i<n_cells;i++) cgpw.m_cell_weights[i].clear();

#       pragma omp barrier

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,block_loc, schedule(dynamic) reduction(+:n_total_nbh,n_positive_weights) )
        {
          IJK loc_a = block_loc + gl;
          size_t cell_a = grid_ijk_to_index( dims , loc_a );

          const uint64_t* idmol_a = cells[cell_a].field_pointer_or_null(field_idmol);
          const uint64_t* id_a    = cells[cell_a][field::id];
          const uint8_t* type_a   = cells[cell_a][field::type];

          cgpw.m_cell_weights[cell_a].set_number_of_particles( cells[cell_a].size() );
          const auto& sp = *species;
  
          apply_cell_particle_neighbors(*grid, *chunk_neighbors, cell_a, loc_a, std::false_type() /* not symetric */,
            [cells,cell_a,no_intramolecular,field_idmol,&cgpw,idmol_a,id_a,type_a,&molid_weight_map,&sp ,&n_total_nbh, &n_positive_weights, &chemPotMap, &mcp]
            ( unsigned int p_a, size_t cell_b, unsigned int p_b , size_t p_nbh_index )
            {
#             ifndef NDEBUG            
              const uint64_t* idmol_b = cells[cell_b].field_pointer_or_null(field_idmol);
#             endif
              const uint64_t* id_b    = cells[cell_b][field::id];

              // assert( cell_a>=0 && cell_a<n_cells );
              //assert( p_nbh_index == nbh_weight[cell_a].size() );
              //assert( p_nbh_index == cgpw.m_cell_weights[cell_a].size() );

              onika::oarray_t<uint64_t,2> pair = { id_a[p_a] , id_b[p_b] };
              if( pair[0] > pair[1] )
              {
                pair = { id_b[p_b] , id_a[p_a] };
              }
              assert( pair[0] != pair[1] );

              const auto it = chemPotMap.find( pair );
              double pair_weight = 1.0;
              
              if( no_intramolecular )
              {
                pair_weight = ( it != chemPotMap.end() ) ? 0.0 : 1.0;
              }
              else
              {
                const int moltype_a = molecule_type_from_id( idmol_a[p_a] );
                // Third case : we check weights----------------------------------------------
                if( it != chemPotMap.end() )
                {
                  assert( moltype_a == molecule_type_from_id( idmol_b[p_b] ) );
                
                  if( it->second == 0 ) pair_weight = molid_weight_map[moltype_a].m_bond_weight;
                  else if( it->second == 1 ) pair_weight = molid_weight_map[moltype_a].m_bend_weight;
                  else if( it->second == 2 ) pair_weight = molid_weight_map[moltype_a].m_torsion_weight;
                  else
                  {
                    fatal_error() << "Internal error: inconsistent chemical link type "<<it->second<<" in chemicalMap"<<std::endl;
                  }
                }
                else
                {
                  pair_weight = 1.0 ;
                }                
              }
                            
              ++ n_total_nbh;
              if(pair_weight>0.0) ++n_positive_weights;
              
              cgpw.m_cell_weights[cell_a].set_neighbor_weight( p_a , p_nbh_index , pair_weight );
            }
          );

        }
        GRID_OMP_FOR_END
      }

//      T2 = std::chrono::high_resolution_clock::now();
//      lout << "T1-T0 = "<< (T1-T0).count()/1000.0 << ", T2-T1 = "<< (T2-T1).count()/1000.0 << std::endl;
      int ratio_precent = n_total_nbh>0 ? (n_positive_weights*100)/n_total_nbh : 0 ;
      ldbg << "nbh="<<n_total_nbh<<", n_positive_weights="<<n_positive_weights<<", ratio="<< ratio_precent << "%" << std::endl;

    }
  };

  template<class GridT> using MoleculePairWeightChunkTmpl = MoleculePairWeightChunk<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_pair_weight)
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "molecule_pair_weight", make_grid_variant_operator< MoleculePairWeightChunkTmpl > );
  }

}
