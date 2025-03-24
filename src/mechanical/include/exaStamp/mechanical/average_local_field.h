#pragma once

#include <onika/math/basic_types.h>
#include <exanb/compute/compute_pair_buffer.h>

#include <vector>
#include <algorithm>

#include <omp.h>

XNB_DECLARE_FIELD(double          , local_field , "mylocalfield");

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
        
        std::cout << "My average field calculator over " << n << " neighbors" << std::endl;

        double sum_f = 0.;
        for(size_t i=0;i<n;i++)
          {
            sum_f += tab.nbh_pt[i][field::local_field];
            std::cout << "\tneigh val = " << tab.nbh_pt[i][field::local_field] << std::endl;            
          }
        sum_f /= (2.*n);
        std::cout << "\n\t\t average final value = " << sum_f << std::endl;
        
      }
    };
  
}  
