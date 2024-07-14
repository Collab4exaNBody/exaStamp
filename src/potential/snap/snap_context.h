#pragma once

#include <onika/memory/allocator.h>
#include "sna.h"

namespace exaStamp
{

  using namespace exanb;

  struct SnapLMPContext
  {
    SnapExt::SnapConfig m_config;
    LAMMPS_NS::SNA * sna = nullptr;
    onika::memory::CudaMMVector<LAMMPS_NS::SnapScratchBuffers> m_scratch;
    onika::memory::CudaMMVector<double> m_coefs;    // size = number of bispectrum coefficients
    onika::memory::CudaMMVector<double> m_factor;   // size = number of elements
    onika::memory::CudaMMVector<double> m_radelem;  // size = number of elements
    onika::memory::CudaMMVector<double> m_beta;       // size = number of particles
    onika::memory::CudaMMVector<double> m_bispectrum; // size = number of particles
    double m_rcut = 0.0;
    
    inline ~SnapLMPContext() { for(auto& s:m_scratch) s.finalize(); }
  };

}


