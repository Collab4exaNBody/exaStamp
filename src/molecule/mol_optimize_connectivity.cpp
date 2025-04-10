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
      
#     pragma omp parallel
      {
        size_t c,p; // unused
        unsigned int ta, tb, tc, td;

#       pragma omp for schedule(static) nowait
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
          const auto it = molecule_compute_parameters->m_intramol_param_map.find( bond_key(ta,tb) );
          intramolecular_parameters->m_bond_param_idx[i] = ( it != molecule_compute_parameters->m_intramol_param_map.end() ) ? it->second : -1;
        }

#       pragma omp for schedule(static) nowait
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
          const auto it = molecule_compute_parameters->m_intramol_param_map.find( angle_key(ta,tb,tc) );
          intramolecular_parameters->m_angle_param_idx[i] = ( it != molecule_compute_parameters->m_intramol_param_map.end() ) ? it->second : -1;
        }

#       pragma omp for schedule(static) nowait
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
          const auto it = molecule_compute_parameters->m_intramol_param_map.find( torsion_key(ta,tb,tc,td) );
          intramolecular_parameters->m_torsion_param_idx[i] = ( it != molecule_compute_parameters->m_intramol_param_map.end() ) ? it->second : -1;
        }

#       pragma omp for schedule(static) nowait
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
          const auto it = molecule_compute_parameters->m_intramol_param_map.find( improper_key(ta,tb,tc,td) );
          intramolecular_parameters->m_improper_param_idx[i] = ( it != molecule_compute_parameters->m_intramol_param_map.end() ) ? it->second : -1;
        }
      }

    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(mol_optimize_connectivity)
  {
    /* ', field::_idmol' : this ensures that only grids with idmol field will be accepted to instantiate this operator */
    OperatorNodeFactory::instance()->register_factory( "mol_optimize_connectivity", make_simple_operator< MolOptimizeConnectivity > );
  }

}
