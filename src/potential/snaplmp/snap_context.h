#pragma once

namespace exaStamp
{

  using namespace exanb;
//  using onika::memory::DEFAULT_ALIGNMENT;
//  using namespace SnapExt;

  struct SnapLMPThreadContext
  {
    LAMMPS_NS::SNA * sna = nullptr;
    double * beta = nullptr;
    double * bispectrum = nullptr;
  };

  struct SnapLMPContext
  {
    LAMMPS_NS::LAMMPS * ptr = nullptr;
    std::vector<SnapLMPThreadContext> m_thread_ctx;
    SnapExt::SnapConfig m_config;
    std::vector<double> m_coefs;    // size = number of bispectrum coefficients
    std::vector<double> m_factor;   // size = number of elements
    std::vector<double> m_radelem;  // size = number of elements
    
    std::vector<double> m_beta;       // size = number of particles
    std::vector<double> m_bispectrum; // size = number of particles
    double m_rcut = 0.0;
  };

}


