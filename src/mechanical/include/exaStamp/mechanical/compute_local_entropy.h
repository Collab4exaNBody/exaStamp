#pragma once

#include <onika/math/basic_types.h>
#include <exanb/compute/compute_pair_buffer.h>

#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <vector>
#include <algorithm>
#include <omp.h>

XNB_DECLARE_FIELD(double, local_entropy, "jesaispascequejefais");

namespace exaStamp {

using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

struct alignas(DEFAULT_ALIGNMENT) LocalEntropyOp {
  const double m_rcut;
  const ThermodynamicState& m_thermo_state;

  double m_sigma = 0;

  double m_dr = 0;
  size_t m_nbins = 0;
  double m_2_sigma_sq = 0;

  double m_volume = 0;
  size_t m_particle_count = 0;
  double m_density = 0;

  std::vector<double> m_rbin = {0};
  std::vector<double> m_rbinsq = {0};
  std::vector<double> m_prefactor = {0};

  bool m_allocated = false;

  LocalEntropyOp(const double rcut, const ThermodynamicState& thermo) : m_rcut(rcut), m_thermo_state(thermo) {
    if (!m_allocated)
      allocate();
  };

private:
  void allocate() {
    if (m_allocated)
      return;

    m_sigma = 0.15;
    m_nbins = static_cast<size_t>(m_rcut / m_sigma) + 1;
    m_dr = m_sigma;

    m_volume = m_thermo_state.volume();
    m_particle_count = m_thermo_state.particle_count();
    m_density = static_cast<double>(m_particle_count) / m_volume;

    m_rbin.resize(m_nbins);
    m_rbinsq.resize(m_nbins);
    m_prefactor.resize(m_nbins);

    std::memset(m_rbin.data(), 0., sizeof(double));
    std::memset(m_rbinsq.data(), 0., sizeof(double));
    std::memset(m_prefactor.data(), 0., sizeof(double));

    for (size_t i = 0; i < m_nbins; ++i) {
      m_rbin[i] = i * m_dr;
      m_rbinsq[i] = m_rbin[i] * m_rbin[i];
    }

    for (size_t i = 1; i < m_nbins; ++i) {
      m_prefactor[i] = 1. / ((4. * M_PI * m_density * sqrt(2. * M_PI) * m_sigma) * m_rbinsq[i]);
    }
    m_prefactor[0] = m_prefactor[1];

    // precompute some constants
    m_2_sigma_sq = 2 * m_sigma * m_sigma;

    m_allocated = true;
  }

public:
  template <class CellParticlesT>
  inline void operator()(size_t jnum, ComputePairBuffer2<false, false>& buf, double& local_entropy,
                         CellParticlesT /* cells */) const {
    
    lout << "------" << std::endl;




    std::vector<double> gm(m_nbins, 0.0);
    std::vector<double> integral(m_nbins, 0.0);

    

    for (size_t ii = 0; ii < jnum; ++ii) {
      double r = sqrt(buf.d2[ii]);
      size_t bin = static_cast<size_t>(std::floor(r / m_dr));
      double delta = r - m_rbin[bin];
      gm[bin] += std::min(std::exp(-1. * (delta * delta) / m_2_sigma_sq) * m_prefactor[bin], 1.0e-10);
    }

    for (size_t k = 0; k < m_nbins; ++k) {
      integral[k] = (gm[k] * std::log(gm[k]) - gm[k] + 1.) * m_rbinsq[k];
    }

    double tmp = 0.;
    for (size_t k = 0; k < m_nbins; ++k) {
      lout << k << " " << m_rbinsq[k] << " " << integral[k] << std::endl;
      tmp += integral[k];
    }

    // tmp += 0.5 * integral[0];
    // tmp += 0.5 * integral[m_nbins - 1];
    tmp *= m_dr;

    lout << -2. * M_PI * m_density * tmp << std::endl;
    local_entropy = 0.0;

    lout << "------" << std::endl;
  }
};
} // namespace exaStamp
