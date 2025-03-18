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

  struct alignas(DEFAULT_ALIGNMENT) CentroSymmetryOp
    {
      const double m_rcut;
      const int m_nnn;
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

      double avg_dist = 0.;
      for (size_t i=0;i<n;i++)
        {
          avg_dist += sqrt( buf.drx[i]*buf.drx[i] + buf.dry[i]*buf.dry[i] + buf.drz[i]*buf.drz[i] );
        }
      avg_dist /= n;
      local_structural_data[ buf.cell ].csp[ buf.part ] = avg_dist;

      }
    };
  
}  
