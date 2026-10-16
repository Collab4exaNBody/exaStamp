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
#include <exaStamp/molecule/id_map.h>
#include <exanb/core/particle_id_codec.h>

#include <exaStamp/molecule/mol_connectivity.h>
#include <exaStamp/molecule/molecule_compute_param.h>
#include <exaStamp/molecule/intramolecular_pair_weight.h>

#include <mpi.h>

namespace exaStamp
{

  class MolOptimizeConnectivity : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                    , mpi           , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( ParticleSpecies             , species       , INPUT , REQUIRED );
    ADD_SLOT( IdMap                       , id_map        , INPUT  );
    ADD_SLOT( IdMapGhosts                 , id_map_ghosts , INPUT  );
    ADD_SLOT( IntramolecularPairWeighting , weight        , INPUT , IntramolecularPairWeighting{} );

    ADD_SLOT( ChemicalBonds     , chemical_bonds     , INPUT_OUTPUT );
    ADD_SLOT( ChemicalAngles    , chemical_angles    , INPUT_OUTPUT );    
    ADD_SLOT( ChemicalTorsions  , chemical_torsions  , INPUT_OUTPUT );
    ADD_SLOT( ChemicalImpropers , chemical_impropers , INPUT_OUTPUT );

    ADD_SLOT( MoleculeComputeParameterSet      , molecule_compute_parameters , INPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( IntramolecularParameterIndexLists, intramolecular_parameters , INPUT_OUTPUT, DocString{"Intramolecular functional parmater index lists"} );

  public:
    inline void execute () override final
    {
      const size_t n_bonds = chemical_bonds->size();
      const size_t n_bends = chemical_angles->size();
      const size_t n_torsions = chemical_torsions->size();
      const size_t n_impropers = chemical_impropers->size();
      
      intramolecular_parameters->m_bond_param_idx.assign( n_bonds , -1 );
      intramolecular_parameters->m_angle_param_idx.assign( n_bends , -1 );
      intramolecular_parameters->m_torsion_param_idx.assign( n_torsions , -1 );
      intramolecular_parameters->m_improper_param_idx.assign( n_impropers , -1 );
      
      const auto & param_map = molecule_compute_parameters->m_intramol_param_map;
      const auto & null_keys = molecule_compute_parameters->m_null_param_keys;
      unsigned long long n_missing_bond=0, n_missing_angle=0, n_missing_torsion=0, n_missing_improper=0;

#     pragma omp parallel
      {
        size_t c,p; // unused
        unsigned int ta, tb, tc, td;

#       pragma omp for schedule(static) reduction(+:n_missing_bond) nowait
        for(size_t i=0;i<n_bonds;i++)
        {
          auto& b = (*chemical_bonds)[i];
          b[0] = atom_from_idmap( b[0] , *id_map , *id_map_ghosts );
          b[1] = atom_from_idmap( b[1] , *id_map , *id_map_ghosts );
          // to get the id of the particles types ta and tb
          decode_cell_particle( /*IN*/ b[0], /*OUT*/ c,p, ta );
          decode_cell_particle( /*IN*/ b[1], /*OUT*/ c,p, tb );
	  //ldbg << "In mol_optimize_connectivity.cpp bond full atome types id ta = " << ta << " tb = " << tb << std::endl;
          // to convert the ids in particles FFtypes, same name ta and tb
	  ta = species->at(ta).m_FFtypeId;
	  tb = species->at(tb).m_FFtypeId;
	  //ldb << "In mol_optimize_connectivity.cpp bond FF types id ta = " << ta << " tb = " << tb << std::endl;
          const uint64_t key = bond_key(ta,tb);
          const auto it = param_map.find( key );
          if( it != param_map.end() ) intramolecular_parameters->m_bond_param_idx[i] = it->second;
          else if( null_keys.find( key ) == null_keys.end() ) ++ n_missing_bond;
        }

#       pragma omp for schedule(static) reduction(+:n_missing_angle) nowait
        for(size_t i=0;i<n_bends;i++)
        {
          auto& b = (*chemical_angles)[i];
          b[0] = atom_from_idmap( b[0] , *id_map , *id_map_ghosts );
          b[1] = atom_from_idmap( b[1] , *id_map , *id_map_ghosts );
          b[2] = atom_from_idmap( b[2] , *id_map , *id_map_ghosts );
	  // to get the id of the particles types ta, tb, and tc
          decode_cell_particle( b[0], c,p, ta );
          decode_cell_particle( b[1], c,p, tb );
          decode_cell_particle( b[2], c,p, tc );
	  // to convert the ids in particles FFtypes
          ta = species->at(ta).m_FFtypeId;
          tb = species->at(tb).m_FFtypeId;
          tc = species->at(tc).m_FFtypeId;
          const uint64_t key = angle_key(ta,tb,tc);
          const auto it = param_map.find( key );
          if( it != param_map.end() ) intramolecular_parameters->m_angle_param_idx[i] = it->second;
          else if( null_keys.find( key ) == null_keys.end() ) ++ n_missing_angle;
        }

#       pragma omp for schedule(static) reduction(+:n_missing_torsion) nowait
        for(size_t i=0;i<n_torsions;i++)
        {
          auto& t = (*chemical_torsions)[i];
          t[0] = atom_from_idmap( t[0] , *id_map , *id_map_ghosts );
          t[1] = atom_from_idmap( t[1] , *id_map , *id_map_ghosts );
          t[2] = atom_from_idmap( t[2] , *id_map , *id_map_ghosts );
          t[3] = atom_from_idmap( t[3] , *id_map , *id_map_ghosts );
	  // to get the id of the particles types ta, tb, tc and td
          decode_cell_particle( t[0], c,p, ta );
          decode_cell_particle( t[1], c,p, tb );
          decode_cell_particle( t[2], c,p, tc );
          decode_cell_particle( t[3], c,p, td );
	  // to convert the ids in particles FFtypes
          ta = species->at(ta).m_FFtypeId;
          tb = species->at(tb).m_FFtypeId;
          tc = species->at(tc).m_FFtypeId;
          td = species->at(td).m_FFtypeId;
          const uint64_t key = torsion_key(ta,tb,tc,td);
          const auto it = param_map.find( key );
          if( it != param_map.end() ) intramolecular_parameters->m_torsion_param_idx[i] = it->second;
          else if( null_keys.find( key ) == null_keys.end() ) ++ n_missing_torsion;
        }

#       pragma omp for schedule(static) reduction(+:n_missing_improper) nowait
        for(size_t i=0;i<n_impropers;i++)
        {
          auto& t = (*chemical_impropers)[i];
          t[0] = atom_from_idmap( t[0] , *id_map , *id_map_ghosts );
          t[1] = atom_from_idmap( t[1] , *id_map , *id_map_ghosts );
          t[2] = atom_from_idmap( t[2] , *id_map , *id_map_ghosts );
          t[3] = atom_from_idmap( t[3] , *id_map , *id_map_ghosts );
	  // to get the id of the particles types ta, tb, tc and td
          decode_cell_particle( t[0], c,p, ta );
          decode_cell_particle( t[1], c,p, tb );
          decode_cell_particle( t[2], c,p, tc );
          decode_cell_particle( t[3], c,p, td );
	  // to convert the ids in particles FFtypes
          ta = species->at(ta).m_FFtypeId;
          tb = species->at(tb).m_FFtypeId;
          tc = species->at(tc).m_FFtypeId;
          td = species->at(td).m_FFtypeId;
          const uint64_t key = improper_key(ta,tb,tc,td);
          const auto it = param_map.find( key );
          if( it != param_map.end() ) intramolecular_parameters->m_improper_param_idx[i] = it->second;
          else if( null_keys.find( key ) == null_keys.end() ) ++ n_missing_improper;
        }
      }

      // terms without parameters are silently skipped by force computations : warn once
      // (flag has the same value on all processes, so the collective call stays consistent)
      if( ! m_missing_params_reported )
      {
        unsigned long long n_missing[4] = { n_missing_bond, n_missing_angle, n_missing_torsion, n_missing_improper };
        MPI_Allreduce( MPI_IN_PLACE, n_missing, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, *mpi );
        if( ( n_missing[0] + n_missing[1] + n_missing[2] + n_missing[3] ) > 0 )
        {
          lout << "WARNING: intramolecular terms without potential parameters are ignored : bonds="<<n_missing[0]
               <<" , angles="<<n_missing[1]<<" , torsions="<<n_missing[2]<<" , impropers="<<n_missing[3]
               <<" (terms involving ghost atoms may be counted on several MPI processes)" << std::endl;
          m_missing_params_reported = true;
        }
      }
    }

  private:
    bool m_missing_params_reported = false;
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(mol_optimize_connectivity)
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "mol_optimize_connectivity", make_simple_operator< MolOptimizeConnectivity > );
  }

}
