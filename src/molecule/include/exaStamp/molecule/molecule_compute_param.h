#pragma once

#include <exaStamp/molecule/molecule_species.h>
#include <exanb/core/basic_types_def.h>
#include <onika/oarray.h>
#include <exaStamp/molecule/pair_potential_parameters.h>

namespace exaStamp
{
  using MoleculeGenericFuncParam = onika::oarray_t<double,4>;

  struct MoleculeComputeParams
  {
    unsigned int m_nb_pairs = 0;
    unsigned int m_nb_bonds = 0;
    unsigned int m_nb_bends = 0;
    unsigned int m_nb_torsions = 0;
    unsigned int m_nb_impropers = 0;

    uint64_t m_pairs[MAX_MOLECULE_PAIRS];
    uint64_t m_torsions[MAX_MOLECULE_TORSIONS]; // contains atom places and func parameter index
    uint64_t m_impropers[MAX_MOLECULE_IMPROPERS];
    uint32_t m_bonds[MAX_MOLECULE_BONDS];
    uint32_t m_bends[MAX_MOLECULE_BENDS];
  };
  
  struct MoleculeSetComputeParams
  {
    std::vector<MoleculeComputeParams> m_molecules;
    std::vector<MoleculeGenericFuncParam> m_func_params;
    std::vector<IntramolecularRFParam> m_rf_params;
    std::vector<IntramolecularLJExp6Param> m_ljexp6_params;
  };

}

