/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/
#pragma once

// MiLaDy's 2-body kernel patch (lib/milady-cpp compute_kernel_nbody.cpp) as a pair potential,
// one instance per element pair (compute_force_pair_multimat or k2b_v2_multi_force):
//
//   V_ab(r) = F * delta^2 * f(r) * sum_q w_q f(z_q) exp( -(r - z_q)^2 / (2 sigma^2) )
//
//   w_q  : MiLaDy weights of this element pair (one per radial point z_q)
//   f    : MiLaDy fcut_rij of type fcut_type (1, 2 or 3) with r_cut, r_cut_width
//   F    : pair_factor. MiLaDy adds the kernel from both atoms, and counts a same-element pair
//          twice more (tau = 2), so F = 4 for a-a and 2 for a-b. Default (0) = auto: 4 when both
//          species have the same z, else 2.
//
// exaStamp shifts pair energies by V(rcut): use MiLaDy's descriptor cutoff (where f = 0 for the usual
// r_cut <= rcut) to keep MiLaDy's energies.
// Same layout as pair potential "k2b": owning K2bV2Parms (CudaMMVector arrays) and the device-copyable
// read-only view K2bV2ParmsRO. k2b_v2_multi_force and k2b_v2_force run on GPU; compute_force_pair_multimat on CPU.
//
// Unlike pair potential "k2b": per element-pair weights, arbitrary centres z_q, smooth cutoff
// f(r) and weighting f(z_q), amplitude delta^2. Get the parameters from a MiLaDy .pot with
// milady_lammps/exastamp's milady_pair_yaml tool.
//
// YAML (lengths in angstrom unless given as quantities):
//   parameters: { sigma: 0.2 ang, delta: 1.0, fcut_type: 3, r_cut: 5.0 ang, r_cut_width: 0.5 ang,
//                 zr: [z_1, ..., z_np], w: [w_1, ..., w_np] }        # optional: pair_factor: 4

#include <cmath>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <onika/memory/allocator.h>

namespace exaStamp
{
  struct K2bV2Parms
  {
    double sigma = 0.2;
    double delta = 1.0;
    int fcut_type = 3;
    double r_cut = 0.0;
    double r_cut_width = 0.0;
    double pair_factor = 0.0;                   // 0: auto (4 same z, 2 otherwise)
    onika::memory::CudaMMVector<double> zr;     // centres z_q (angstrom)
    onika::memory::CudaMMVector<double> cw;     // delta^2 f(z_q) w_q
  };

  // device-copyable read-only view of K2bV2Parms
  struct K2bV2ParmsRO
  {
    const double* __restrict__ zr = nullptr;
    const double* __restrict__ cw = nullptr;
    int np = 0;
    double sigma = 0.2;
    double delta = 1.0;
    int fcut_type = 3;
    double r_cut = 0.0;
    double r_cut_width = 0.0;
    double pair_factor = 0.0;

    K2bV2ParmsRO() = default;
    K2bV2ParmsRO(const K2bV2ParmsRO&) = default;
    K2bV2ParmsRO& operator = (const K2bV2ParmsRO&) = default;

    inline K2bV2ParmsRO(const K2bV2Parms& p)
      : zr(p.zr.data()), cw(p.cw.data()), np(static_cast<int>(p.zr.size())), sigma(p.sigma), delta(p.delta)
      , fcut_type(p.fcut_type), r_cut(p.r_cut), r_cut_width(p.r_cut_width), pair_factor(p.pair_factor)
    {}
  };

  // MiLaDy neighbours.cpp fcut_rij, types 1..3 (verbatim, including type 1's derivative)
  ONIKA_HOST_DEVICE_FUNC inline void k2b_v2_fcut(int type, double r, double rc, double rw, double& f, double& df)
  {
    const double pi = 3.14159265358979323846;
    f = 0.; df = 0.;
    if( r > rc ) return;
    if( type == 1 ) { const double xx = (r / rc) * (r / rc) - 1.0; f = xx * xx; df = 4.0 * xx / (rc * rc); }
    else if( type == 2 ) { f = 0.5 * (cos(pi * r / rc) + 1.0); df = -0.5 * pi * sin(pi * r / rc) / rc; }
    else if( type == 3 )
    {
      if( r < rc - rw ) f = 1.;
      else { const double xx = pi * (r - rc + rw) / rw; f = 0.5 * (cos(xx) + 1.0); df = -0.5 * sin(xx) * pi / rw; }
    }
  }
}

namespace onika { namespace cuda {
  template<> struct ReadOnlyShallowCopyType<exaStamp::K2bV2Parms> { using type = exaStamp::K2bV2ParmsRO; };
} }

namespace YAML
{
  template<> struct convert<exaStamp::K2bV2Parms>
  {
    static bool decode(const Node& node, exaStamp::K2bV2Parms& v)
    {
      using onika::physics::Quantity;
      v = exaStamp::K2bV2Parms{};
      if( !node.IsMap() || !node["sigma"] || !node["r_cut"] || !node["zr"] || !node["w"] ) { return false; }
      v.sigma = node["sigma"].as<Quantity>().convert();
      if( node["delta"] ) v.delta = node["delta"].as<double>();
      if( node["fcut_type"] ) v.fcut_type = node["fcut_type"].as<int>();
      if( v.fcut_type < 1 || v.fcut_type > 3 ) { return false; }
      v.r_cut = node["r_cut"].as<Quantity>().convert();
      if( node["r_cut_width"] ) v.r_cut_width = node["r_cut_width"].as<Quantity>().convert();
      if( v.fcut_type == 3 && !( v.r_cut_width > 0. ) ) { return false; }
      if( node["pair_factor"] ) v.pair_factor = node["pair_factor"].as<double>();
      const auto zr = node["zr"].as< std::vector<double> >();
      const auto w = node["w"].as< std::vector<double> >();
      if( zr.size() != w.size() || zr.empty() ) { return false; }
      v.zr.assign( zr.begin(), zr.end() );
      v.cw.resize( w.size() );
      for( size_t q = 0; q < w.size(); ++q )
      {
        double fz, dfz;
        exaStamp::k2b_v2_fcut( v.fcut_type, zr[q], v.r_cut, v.r_cut_width, fz, dfz );
        v.cw[q] = v.delta * v.delta * fz * w[q];
      }
      return true;
    }
  };
}

namespace exaStamp
{
  // Pair energy e and de/dr (eV, eV/angstrom before conversion to internal units).
  ONIKA_HOST_DEVICE_FUNC inline void k2b_v2_compute_energy(const K2bV2ParmsRO& p, const PairPotentialMinimalParameters& pair, double r, double& e, double& de)
  {
    e = 0.; de = 0.;
    double f, df;
    k2b_v2_fcut( p.fcut_type, r, p.r_cut, p.r_cut_width, f, df );
    if( f == 0. && df == 0. ) return;
    const double inv_s2 = 1. / ( p.sigma * p.sigma );
    double g = 0., dg = 0.;
    for( int q = 0; q < p.np; ++q )
    {
      const double dr = r - p.zr[q];
      const double t = p.cw[q] * exp( -0.5 * dr * dr * inv_s2 );
      g += t;
      dg -= t * dr * inv_s2;
    }
    const double F = ( p.pair_factor != 0. ) ? p.pair_factor : ( pair.m_atom_a.m_z == pair.m_atom_b.m_z ? 4.0 : 2.0 );
    e = F * f * g;
    de = F * ( df * g + f * dg );
    static const double conv_energy_inv = 1e-4 * onika::physics::elementaryCharge / onika::physics::atomicMass;   // eV -> internal
    e *= conv_energy_inv;
    de *= conv_energy_inv;
  }

  // host overload on the owning type (singlemat / multimat operators)
  inline void k2b_v2_compute_energy(const K2bV2Parms& p, const PairPotentialMinimalParameters& pair, double r, double& e, double& de)
  {
    k2b_v2_compute_energy( K2bV2ParmsRO(p), pair, r, e, de );
  }
}

#define USTAMP_POTENTIAL_NAME     k2b_v2
#define USTAMP_POTENTIAL_PARAMS   K2bV2Parms
#define USTAMP_POTENTIAL_COMPUTE  k2b_v2_compute_energy
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) onika::make_flat_tuple(p.m_atom_a.m_z,p.m_atom_b.m_z)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1
