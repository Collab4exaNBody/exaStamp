#pragma once

#include <onika/math/basic_types.h>
#include <exanb/compute/compute_pair_buffer.h>

#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <vector>
#include <algorithm>

#include <omp.h>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;  

  struct alignas(DEFAULT_ALIGNMENT) NumNeighborsOp
    {
      const double m_rcut;
      GridParticleLocalStructuralMetrics & m_local_structural_data;
      
      template<class CellParticlesT>
    	inline void operator ()
        (
    	 size_t n,
    	 ComputePairBuffer2<false,false>& buf,
    	 CellParticlesT
    	 ) const	
      {
    	GridParticleLocalStructuralMetrics & local_structural_data = m_local_structural_data;

      local_structural_data[ buf.cell ].numneighbors[ buf.part ] = n;

      }
    };
  
}  
