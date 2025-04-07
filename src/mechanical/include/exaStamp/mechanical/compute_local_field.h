#pragma once

#include <onika/math/basic_types.h>
#include <exanb/compute/compute_pair_buffer.h>
#include <omp.h>

XNB_DECLARE_FIELD(double, local_field, "local_field");

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;  
  
  struct alignas(DEFAULT_ALIGNMENT) LocalFieldOp
    {
      const double m_rcut;

      template <class CellParticlesT>
      inline void operator()(size_t n, ComputePairBuffer2<false, false>& tab, double& local_field,
                             CellParticlesT) const {
        lout << "My local field calculator over " << n << " neighbors" << std::endl;
        local_field = 3.14 * n;
      }
  };
}  
