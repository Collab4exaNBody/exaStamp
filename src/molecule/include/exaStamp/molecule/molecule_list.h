#pragma once

#include <cstdint>
#include <onika/memory/allocator.h>

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

}

