#pragma once

#include <cstdint>
#include <onika/memory/allocator.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>

namespace exaStamp
{  
  struct MoleculeLists
  {
    // start of ith molecule data in m_molecule_data;
    onika::memory::CudaMMVector<uint64_t> m_offset;

    // serialized molecules' data
    // [ MolType , AtomId... ]
    onika::memory::CudaMMVector<int64_t> m_data;
    
    inline size_t number_of_molecules() const { return m_offset.size(); }
    
    inline unsigned int type(size_t m) const { return m_data[ m_offset[m] ]; }

    inline uint64_t atom_data(size_t m, size_t a) const { return m_data[ m_offset[m] + 1 + a]; }
  };

  struct MoleculeListsRO
  {
    const uint64_t * __restrict__ m_offset = nullptr;
    const int64_t * __restrict__ m_data = nullptr;
    
    MoleculeListsRO() = default;
    MoleculeListsRO(const MoleculeListsRO&) = default;
    MoleculeListsRO(MoleculeListsRO&&) = default;
    MoleculeListsRO& operator = (const MoleculeListsRO&) = default;
    MoleculeListsRO& operator = (MoleculeListsRO&&) = default;

    inline MoleculeListsRO( const MoleculeLists& ml )
      : m_offset( ml.m_offset.data() )
      , m_data( ml.m_data.data() )
      {}
   
    ONIKA_HOST_DEVICE_FUNC  
    inline unsigned int type(size_t m) const { return m_data[ m_offset[m] ]; }
    
    ONIKA_HOST_DEVICE_FUNC
    inline uint64_t atom_data(size_t m, size_t a) const { return m_data[ m_offset[m] + 1 + a]; }
  };

}

// specialize ReadOnlyShallowCopyType so MoleculeListsRO is the read only type for MoleculeLists
namespace onika { namespace cuda { template<> struct ReadOnlyShallowCopyType< exaStamp::MoleculeLists > { using type = exaStamp::MoleculeListsRO; }; } }

