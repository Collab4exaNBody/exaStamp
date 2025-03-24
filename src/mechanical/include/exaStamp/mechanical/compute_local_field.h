#pragma once

#include <onika/math/basic_types.h>
#include <exanb/compute/compute_pair_buffer.h>

#include <vector>
#include <algorithm>

#include <omp.h>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;  
  
  struct alignas(DEFAULT_ALIGNMENT) LocalFieldOp
    {
      const double m_rcut;
      template<class CellParticlesT>
    	inline void operator ()
        (
         size_t n,
         ComputePairBuffer2<false, false>& tab,
         double &local_field,
         
         CellParticlesT
         ) const	
      {
        
        std::cout << "My local field calculator over " << n << " neighbors" << std::endl;
       
        //        std::cout << "local field = " << n << std::endl;
        local_field = 3.14*n;
        
      }
    };
  
}  
