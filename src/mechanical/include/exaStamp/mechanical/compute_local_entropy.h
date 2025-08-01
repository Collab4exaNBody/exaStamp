#pragma once

#include <exanb/compute/compute_pair_buffer.h>
#include <onika/math/basic_types.h>
#include <onika/string_utils.h>

#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <omp.h>
#include <vector>

namespace exaStamp {

using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

struct alignas(DEFAULT_ALIGNMENT) LocalEntropyOp {

  bool verbose = true;
  double rcut;
  double sigma;
  size_t nbins;
  bool local;

  double inv_sigma = 0.;
  double drm = 0.;
  double inv_drm = 0.;

  size_t particle_count = 0;
  double volume = 0;
  double rho_global = 0;
  double local_volume = 0;
  double inv_local_volume = 0;

  double sigmasq2 = 0;
  double inv_sigmasq2 = 0;
  double deltabin = 0;

  constexpr static double gtol = 1.0e-10;
  constexpr static double eps = 1.0e-8;

  std::vector<double> rbin = {0};
  std::vector<double> rbinsq = {0};
  std::vector<double> prefactor = {0};

  void initialize(const ThermodynamicState& thermo_state) {

    if (!(nbins > 0)) {
      nbins = static_cast<size_t>(rcut / sigma) + 1;
    }

    drm = rcut / (nbins - 1);
    inv_drm = 1. / drm;

    volume = thermo_state.volume();
    particle_count = thermo_state.particle_count();
    rho_global = static_cast<double>(particle_count) / volume;
    local_volume = (4. / 3.) * M_PI * rcut * rcut * rcut;
    inv_local_volume = 1. / local_volume;

    // precompute some variables
    sigmasq2 = 2. * sigma * sigma;
    inv_sigmasq2 = 1. / sigmasq2;

    rbin.resize(nbins);
    rbinsq.resize(nbins);
    prefactor.resize(nbins);

    std::memset(rbin.data(), 0., sizeof(double));
    std::memset(rbinsq.data(), 0., sizeof(double));
    std::memset(prefactor.data(), 0., sizeof(double));

    // at k = 0, r = r^2 = 0.0
    for (size_t k = 1; k < nbins; ++k) {
      rbin[k] = k * drm;
      rbinsq[k] = k * k * drm * drm;
      prefactor[k] = 1. / ((4. * M_PI * rho_global * sqrt(2. * M_PI) * sigma) * rbinsq[k]);
    }
    prefactor[0] = prefactor[1]; // ensure no division by zero

    // evaluate exponential term to compute detltabin
    for (size_t k = 0; k < nbins; ++k) {
      double gk = std::exp(-rbinsq[k] * inv_sigmasq2);
      if (gk < 1.0e-5) {
        break;
      }
      deltabin++;
    }
  }

  template <class CellParticlesT>
  inline void operator()(size_t n, ComputePairBuffer2<false, false>& buf, double& entropy,
                         CellParticlesT /* cells */) const {

    double rho = local ? (static_cast<double>(n) / local_volume) : rho_global;
    double rhofac = rho_global / rho;

    std::vector<double> gm(nbins, 0.0);

    for (size_t i = 0; i < n; ++i) {
      double rij = std::sqrt(buf.d2[i]);
      size_t bin = static_cast<size_t>(rij * inv_drm + eps);

      size_t minbin = (bin > deltabin) ? bin - deltabin : 0;
      size_t maxbin = (nbins < bin + deltabin) ? nbins : bin + deltabin;

      for (size_t k = minbin; k < maxbin; ++k) {
        double dr = rbin[k] - rij;
        gm[k] += std::exp(-dr * dr * inv_sigmasq2) * prefactor[k] * rhofac;
      }
    }

    // compute the intergral with trapzeoid rule
    double tmp = 0.;
    for (size_t k = 0; k < nbins; ++k) {
      double g = gm[k];
      double f = (g >= gtol) ? (g * std::log(g) - g + 1.0) : 1.0;
      double w = (k == 0 || k == nbins - 1) ? 0.5 : 1.0;
      tmp += f * rbinsq[k] * w;
    }

    // don't forget to mulitply the integral by dr
    entropy = -2. * M_PI * rho * tmp * drm;
  }
};

} // namespace exaStamp
