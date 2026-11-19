/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

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
#include <exaStamp/molecule/intramolecular_image.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <mpi.h>

namespace exaStamp
{

  template<class GridT>
  class IntramolecularPairPotentialList : public OperatorNode
  {
    ADD_SLOT( GridT                       , grid          , INPUT , REQUIRED );
    ADD_SLOT( Domain                      , domain        , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies             , species       , INPUT , REQUIRED );
    ADD_SLOT( IdMap                       , id_map        , INPUT , DocString{"Owned atom ids to grid locations"} );
    ADD_SLOT( IdMapGhosts                 , id_map_ghosts , INPUT , DocString{"Ghost atom ids to grid locations"} );
    ADD_SLOT( IntramolecularPairWeighting , weight        , INPUT , IntramolecularPairWeighting{} , DocString{"Not used"} );
    ADD_SLOT( MoleculeComputeParameterSet      , molecule_compute_parameters , INPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( IntramolecularParameterIndexLists, intramolecular_parameters , INPUT_OUTPUT, DocString{"Intramolecular functional parameter index lists"} );

    ADD_SLOT( ChemicalBonds      , chemical_bonds        , INPUT_OUTPUT , DocString{"Bonds (pairs of bonded atoms)"} );
    ADD_SLOT( ChemicalAngles     , chemical_angles       , INPUT_OUTPUT , DocString{"Bending angles (atom triplets, center atom second)"} );    
    ADD_SLOT( ChemicalTorsions   , chemical_torsions     , INPUT_OUTPUT , DocString{"Torsions (4 atoms chained by bonds)"} );
    ADD_SLOT( ChemicalPairs      , chemical_pairs        , INPUT_OUTPUT , DocString{"Intramolecular pairs (atoms 1-2, 1-3 and 1-4 apart)"} );
    ADD_SLOT( ChemicalPairPotMap , chemical_pair_pot_map , INPUT_OUTPUT , DocString{"Map of intramolecular pair potential parameter set and associated weights"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Builds the list of intramolecular pairs (chemical_pairs) : atoms 1-2, 1-3 and 1-4 apart in the bond graph, from the bond,
angle and torsion lists. A pair found in several lists is counted once, with the shortest path. chemical_pair_pot_map
holds these pairs, to exclude them from the neighbor lists (extramolecular_neighbors).
)EOF";
    }

    inline void execute () override final
    {
      const size_t n_bonds = chemical_bonds->size();
      const size_t n_bends = chemical_angles->size();
      const size_t n_torsions = chemical_torsions->size();

      // index chemical links by their extrema to rapidly find if an atom pair belongs to a chemical link.
      // The map is a persistent slot : it must be rebuilt from scratch, otherwise pairs whose atoms have
      // migrated away from this MPI process accumulate and later can't be found in id_map / id_map_ghosts.
      chemical_pair_pot_map->clear();

      ProfilingTimer T;
      profiling_timer_start( T );

      // pairs excluded from the neighbor lists and weighted, keyed by original ids (min id, max id). A pair found in several
      // lists keeps the shortest path (0 = 1-2 bond, 1 = 1-3 angle, 2 = 1-4 torsion). Concurrent insertions are safe.

      // Helper: insert with canonical key; if key already exists, keep the lower priority value (shorter path).
      // We lock the meta-bucket manually (same pattern as MultiThreadedConcurrentMap::insert()) to make
      // the find+insert/update atomic and thread-safe.
      auto try_insert = [&](uint64_t a, uint64_t b, int priority)
      {
        if (a > b) std::swap(a, b);  // canonical: smaller id first, deduplicates (i,j)/(j,i)
        const ChemicalPairPotMap::key_type key = {a, b};
        const size_t mb = ChemicalPairPotMap::meta_bucket( key );
        chemical_pair_pot_map->m_meta_bucket_locks[ mb ].lock();
        auto& bucket = chemical_pair_pot_map->m_meta_bucket[ mb ];
        auto [it, inserted] = bucket.insert( { key, priority } );
        if ( !inserted && priority < it->second )
        {
          it->second = priority;  // shorter path takes precedence
        }
        chemical_pair_pot_map->m_meta_bucket_locks[ mb ].unlock();
      };

#     pragma omp parallel
      {
        // 1-4 torsion pairs  (priority 2)
#       pragma omp for schedule(static) nowait
        for (size_t i = 0; i < n_torsions; i++)
        {
          const auto& torsion = (*chemical_torsions)[i];
          try_insert(torsion[0], torsion[3], 2);
        }

        // 1-3 angle pairs  (priority 1)
#       pragma omp for schedule(static) nowait
        for (size_t i = 0; i < n_bends; i++)
        {
          const auto& angle = (*chemical_angles)[i];
          try_insert(angle[0], angle[2], 1);
        }

        // 1-2 bond pairs  (priority 0)
#       pragma omp for schedule(static) nowait
        for (size_t i = 0; i < n_bonds; i++)
        {
          const auto& bond = (*chemical_bonds)[i];
          try_insert(bond[0], bond[1], 0);
        }
      }

      ldbg << "map build time = " << profiling_timer_elapsed_restart(T) << std::endl;

      const size_t n_chemical_pairs = chemical_pair_pot_map->size();
      chemical_pairs->assign( n_chemical_pairs , ChemicalPair{0,0} );
      intramolecular_parameters->m_pair_param_idx.assign( n_chemical_pairs , -1 );
      intramolecular_parameters->m_pair_image.assign( n_chemical_pairs , 0 );
      static constexpr int pair_anchor[2] = { 0, 0 };
      const IntramolecularImageBuilder image_builder( *domain );
      auto cells = grid->cells();

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
          const Vec3d ra = { cells[c][field::rx][p] , cells[c][field::ry][p] , cells[c][field::rz][p] };
          decode_cell_particle( chemical_pair[1], c, p, tb );
          const Vec3d r[2] = { ra , Vec3d{ cells[c][field::rx][p] , cells[c][field::ry][p] , cells[c][field::rz][p] } };
          if( ! image_builder.build( r, pair_anchor, 2, intramolecular_parameters->m_pair_image[pair_idx] ) )
          {
            fatal_error() << "intramolecular pair atoms are more than "<<IntramolecularImageCodec::MAX_SHIFT<<" periodic images apart" << std::endl;
          }
          const int pair_id = unique_pair_id(ta,tb);
          int param_idx = -1;
          // Since weight is only applied on the Coulombic part of RF, we still need that parameters set to be considered
          param_idx = n_type_pairs * chemicalLinkType + pair_id ;
          intramolecular_parameters->m_pair_param_idx[ pair_idx ] = param_idx;
          ++ pair_idx;
        }
      }
      assert( pair_idx == n_chemical_pairs );

      ldbg << "pairs build time = " << profiling_timer_elapsed_restart(T) << std::endl;
      ldbg << "n_chemical_pairs = "<<n_chemical_pairs<<std::endl;
    }
  };

  template<class GridT> using IntramolecularPairPotentialListTmpl = IntramolecularPairPotentialList<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(intramolecular_pair_list)
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "intramolecular_pair_list", make_grid_variant_operator< IntramolecularPairPotentialListTmpl > );
  }

}
