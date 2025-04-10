#pragma once

#include <onika/math/basic_types.h>
#include <onika/yaml/yaml_enum.h>

#include <exanb/core/grid_fields.h>
#include <exanb/compute/compute_pair_buffer.h>

#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <vector>
#include <algorithm>
#include <omp.h>

EXANB_YAML_ENUM(exaStamp, LocalCentroMethod, GES, MWM);
XNB_DECLARE_FIELD(double, csp, "csp")

namespace exaStamp
{

using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;  

// structure used to sort nearest neighbours
struct neigh_info_t {
  double dist;
  size_t index;
};

struct alignas(DEFAULT_ALIGNMENT) CentroSymmetryOp {

  const double m_rcut;
  const size_t m_nnn;
  const LocalCentroMethod m_method;

  template <class CellParticlesT>
  inline void operator()(size_t n,                              // number of neighbours of atom i
                         ComputePairBuffer2<false, false>& buf, // buffer that contains information about neighbours
                         double& csp,
                         CellParticlesT                         // ??????
  ) const {



    // size_t nn = static_cast<size_t>(this->m_nnn);
    // size_t npairs = nn * (nn  - 1) / 2;
    // size_t nhalf = nn / 2;

    // std::vector<neigh_info_t> neigh;
    // neigh.reserve(n);

    // // Find the N-Nearest neighbours
    // for (size_t i = 0; i < n; ++i) {
    //   neigh.push_back((neigh_info_t) { 
    //     .dist = (buf.drx[i] * buf.drx[i]) + (buf.dry[i] * buf.dry[i]) + (buf.drz[i] * buf.drz[i]),
    //     .index  = i 
    //   });
    // }

    // // sort neighbours by squared distance
    // std::partial_sort(
    //   neigh.begin(), neigh.begin() + nn, neigh.end(),
    //   [] (const neigh_info_t  &a, const neigh_info_t &b) { return a.dist < b.dist; }
    // );

    // if ( m_csp_type == csp_t::ges ) {

    //   calculate_csp_ges(nn, nhalf, npairs, neigh, buf);

    // // } else if ( m_csp_type == csp_t::mwm ) {

    // //   calculate_csp_mwm(nn, nhalf, npairs, neigh, buf);

    // } else {

    //   local_structural_data[ buf.cell ].csp[ buf.part ] = 0.0;

    // }
  }

  // inline void calculate_csp_mwm(
  //   size_t &nn,
  //   size_t &nhalf,
  //   size_t &npairs,
  //   std::vector<neigh_info_t> & neigh,
  //   ComputePairBuffer2<false,false> & buf
  // ) const
  // {
  //   GridParticleLocalStructuralMetrics & local_structural_data = m_local_structural_data;

  //   std::vector<double> cost_matrix(nn * nn);

  //   lout << cost_matrix.size() << std::endl;

  //   // Build the cost matrix
  //   for (size_t j = 0; j < nn; ++j) {
  //     
  //     cost_matrix[j * nn + j] = 100000000000.0;

  //     for (size_t k = j + 1; k < nn; ++k) {

  //       size_t ij = neigh[j].index;
  //       size_t ik = neigh[k].index;

  //       double dx = buf.drx[ij] + buf.drx[ik];
  //       double dy = buf.dry[ij] + buf.dry[ik];
  //       double dz = buf.drz[ij] + buf.drz[ik];

  //       double rjk2 = (dx * dx) + (dy * dy) + (dz * dz);

  //       // adjency matrix is symmetric
  //       cost_matrix[ j * nn + k ] = rjk2;
  //       cost_matrix[ k + nn + j ] = rjk2;

  //     }
  //   }

  //   bool perfect_match = false;
  //   double min_bound = greedy_edge_assignement(nn, cost_matrix, perfect_match);

  //   lout << perfect_match << std::endl;

  //   local_structural_data[ buf.cell ].csp[ buf.part ] = min_bound;

  // }

  // inline double greedy_edge_assignement(size_t nn, const std::vector<double> & weights, bool & perfect_match) const {
  //   std::vector<int> pairs(nn);
  //   double csp = 0.0;

  //   for (int i = 0; i < nn; i++) {
  //     
  //     int k = -1;
  //     double min = 100000.0;

  //     for (int j = 0; j < nn; j++) {
  //       double w = weights[i * nn + j];
  //       if (( i != j ) && ( w < min )) {
  //         min = w;
  //         k = j;
  //       }
  //     }

  //     pairs[i] = k;
  //     csp += min;
  //   }

  //   for (int i = 0; i < nn; ++i)
  //     if (pairs[pairs[i]] != i)
  //       perfect_match = false;

  //   return csp / 2;
  // }

//   inline void calculate_csp_ges(
//     size_t &nn,
//     size_t &nhalf,
//     size_t &npairs,
//     std::vector<neigh_info_t> & neigh,
//     ComputePairBuffer2<false,false> & buf
//   ) const
//   {

//     GridParticleLocalStructuralMetrics & local_structural_data = m_local_structural_data;

//     double csp = 0.0;

//     std::vector<double> pairs;
//     pairs.reserve(npairs);

//     // loop over each npairs (i,j) among nearest neighbours
//     for (size_t j = 0; j < nn; ++j) {
//       for (size_t k = j + 1; k < nn; ++k) {

//         size_t ij = neigh[j].index;
//         size_t ik = neigh[k].index;

//         double dx = buf.drx[ij] + buf.drx[ik];
//         double dy = buf.dry[ij] + buf.dry[ik];
//         double dz = buf.drz[ij] + buf.drz[ik];

//         pairs.push_back((dx * dx) + (dy * dy) + (dz * dz));
//       }
//     }

//     // Find N/2 smallest pair distances
//     std::partial_sort(pairs.begin(), pairs.begin() + nhalf, pairs.end());

//     // compute CSP
//     for (size_t i = 0; i < nhalf; ++i) {
//       csp += pairs[i];
//     }

//     local_structural_data[ buf.cell ].csp[ buf.part ] = csp;
//   }

};

}  
