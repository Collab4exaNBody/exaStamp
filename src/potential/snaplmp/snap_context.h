#pragma once

#include <onika/memory/allocator.h>

namespace exaStamp
{

  using namespace exanb;
//  using onika::memory::DEFAULT_ALIGNMENT;
//  using namespace SnapExt;

  struct SnapLMPThreadContext
  {
    LAMMPS_NS::SnapScratchBuffers * scratch = nullptr;
    double * beta = nullptr;
    double * bispectrum = nullptr;
  };

  struct SnapLMPContext
  {
    SnapExt::SnapConfig m_config;
    LAMMPS_NS::SNA * sna = nullptr;
    onika::memory::CudaMMVector<SnapLMPThreadContext> m_thread_ctx;
    onika::memory::CudaMMVector<double> m_coefs;    // size = number of bispectrum coefficients
    onika::memory::CudaMMVector<double> m_factor;   // size = number of elements
    onika::memory::CudaMMVector<double> m_radelem;  // size = number of elements
    onika::memory::CudaMMVector<double> m_beta;       // size = number of particles
    onika::memory::CudaMMVector<double> m_bispectrum; // size = number of particles
    double m_rcut = 0.0;
  };

}


