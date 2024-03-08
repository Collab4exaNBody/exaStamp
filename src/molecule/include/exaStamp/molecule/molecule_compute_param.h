#pragma once

#include <exaStamp/molecule/molecule_species.h>
#include <exanb/core/basic_types_def.h>
#include <onika/oarray.h>
#include <exaStamp/molecule/pair_potential_parameters.h>

#include <onika/memory/allocator.h>
#include <onika/cuda/ro_shallow_copy.h>

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
    onika::memory::CudaMMVector<MoleculeComputeParams> m_molecules;
    onika::memory::CudaMMVector<MoleculeGenericFuncParam> m_func_params;
    onika::memory::CudaMMVector<IntramolecularRFParam> m_rf_params;
    onika::memory::CudaMMVector<IntramolecularLJExp6Param> m_ljexp6_params;
  };

  struct MoleculeSetComputeParamsRO
  {
    const MoleculeComputeParams * __restrict__ m_molecules = nullptr;
    const MoleculeGenericFuncParam * __restrict__ m_func_params = nullptr;
    const IntramolecularRFParam * __restrict__ m_rf_params = nullptr;
    const IntramolecularLJExp6Param * __restrict__ m_ljexp6_params = nullptr;
    
    MoleculeSetComputeParamsRO() = default;
    MoleculeSetComputeParamsRO(const MoleculeSetComputeParamsRO&) = default;
    MoleculeSetComputeParamsRO(MoleculeSetComputeParamsRO&&) = default;
    MoleculeSetComputeParamsRO& operator = (const MoleculeSetComputeParamsRO&) = default;
    MoleculeSetComputeParamsRO& operator = (MoleculeSetComputeParamsRO&&) = default;

    inline MoleculeSetComputeParamsRO( const MoleculeSetComputeParams& ml )
      : m_molecules( ml.m_molecules.data() )
      , m_func_params( ml.m_func_params.data() )
      , m_rf_params( ml.m_rf_params.data() )
      , m_ljexp6_params( ml.m_ljexp6_params.data() )
      {}

  };

}

// specialize ReadOnlyShallowCopyType so MoleculeListsRO is the read only type for MoleculeLists
namespace onika { namespace cuda { template<> struct ReadOnlyShallowCopyType< exaStamp::MoleculeSetComputeParams > { using type = exaStamp::MoleculeSetComputeParamsRO; }; } }

