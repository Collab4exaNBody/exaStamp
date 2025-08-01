#pragma once

#include <onika/math/basic_types.h>
#include <onika/yaml/yaml_enum.h>

#include <exanb/core/grid_fields.h>
#include <exanb/compute/compute_pair_buffer.h>

#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <vector>
#include <algorithm>
#include <omp.h>

namespace exaStamp {

using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;  


struct alignas(DEFAULT_ALIGNMENT) CentroSymmetryOp {

  using NeighInfo = std::pair<double, size_t>;

  const size_t nnn;

  template <class CellParticlesT>
  inline void operator()(size_t n, ComputePairBuffer2<false, false>& buf, double& csp, CellParticlesT) const {

    csp = 0.0;

    if (n < nnn) {
      return;
    }

    size_t npairs = nnn * (nnn - 1) / 2;
    size_t nhalf = nnn / 2;

    std::vector<NeighInfo> neigh_infos;
    neigh_infos.resize(n);

    // Find the N-Nearest neighbours
    // sort neighbours by squared distance
    for (size_t i = 0; i < n; ++i) {
      neigh_infos[i].first = buf.d2[i];
      neigh_infos[i].second = i;
    }

    std::partial_sort(neigh_infos.begin(), neigh_infos.begin() + nnn, neigh_infos.end(),
                      [](const NeighInfo& a, const NeighInfo& b) { return a.first < b.first; });

    std::vector<double> neigh_pairs;
    neigh_pairs.reserve(npairs);

    // loop over each npairs (i,j) among nearest neighbours
    for (size_t j = 0; j < nnn; ++j) {
      for (size_t k = j + 1; k < nnn; ++k) {

        size_t ij = neigh_infos[j].second;
        size_t ik = neigh_infos[k].second;

        double dx = buf.drx[ij] + buf.drx[ik];
        double dy = buf.dry[ij] + buf.dry[ik];
        double dz = buf.drz[ij] + buf.drz[ik];

        neigh_pairs.push_back((dx * dx) + (dy * dy) + (dz * dz));
      }
    }

    // Find N/2 smallest pair distances
    std::partial_sort(neigh_pairs.begin(), neigh_pairs.begin() + nhalf, neigh_pairs.end());

    // compute CSP
    for (size_t i = 0; i < nhalf; ++i) {
      csp += neigh_pairs[i];
    }
  }
};

} // namespace exaStamp
