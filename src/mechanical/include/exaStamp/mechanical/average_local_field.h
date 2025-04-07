#pragma once

#include <onika/math/basic_types.h>
#include <exanb/compute/compute_pair_buffer.h>
#include <omp.h>

#include <exaStamp/mechanical/compute_local_field.h>

namespace exaStamp
{
  using namespace exanb;

  struct AverageOp
    {
      const double m_rcut;
      template<class ComputePairBufferT, class CellParticlesT>
    	ONIKA_HOST_DEVICE_FUNC inline void operator ()
        (
         size_t n,
         ComputePairBufferT& tab,
         double &local_field,
         
         CellParticlesT
         ) const	
      {
        
        lout << "My average field calculator over " << n << " neighbors" << std::endl;

        double sum_f = 0.;
        for(size_t i=0;i<n;i++)
          {
            sum_f += tab.nbh_pt[i][field::local_field];
            std::cout << "\tneigh val = " << tab.nbh_pt[i][field::local_field] << std::endl;            
          }
        sum_f /= (2.*n);
        lout << "\n\t\t average final value = " << sum_f << std::endl;
        
      }
    };
  
}  
