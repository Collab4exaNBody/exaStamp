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

// ZBL as in MiLaDy (lib/milady-cpp compute_zbl.cpp), as a pair potential.
//
//   phi(r) = Z1 Z2 / (4 pi eps0) * sum_{k=1..3} p_k exp(e_k r / a) / r,  1/a = (Z1^zz + Z2^zz) / pz
//
// Defaults are MiLaDy's screening constants (3 exponentials, not LAMMPS' universal 4-exponential ZBL,
// see pair potential "zbl"). Two ways to reach zero, MiLaDy's zbl_type:
//   mode: switched  (zbl_type 2)  phi(r) * S(r), S = 1 - X^3 (6X^2 - 15X + 10), X = (r - r1) / (r2 - r1),
//                                 i.e. S = 1 below r1, 0 above r2 (C2 at both ends, no energy shift)
//   mode: bridge    (zbl_type 1)  phi(r) up to r1, then the bridge fitted to the ML short-range part
//                                 up to rr, 0 beyond:
//                                   bridge_type exp  : exp(q0 + q1 r + q2 r^2 + q3 r^3)
//                                   bridge_type poly : q0 + q1 r + ... + q5 r^5
//                                 (q in eV and angstrom, as in MiLaDy's params_k2b_to_zbl)
//
// YAML:
//   parameters: { mode: switched, r1: 1.3 ang, r2: 2.4 ang }
//   parameters: { mode: bridge, r1: 1.3 ang, rr: 1.7 ang, bridge_type: exp, bridge: [q0, q1, q2, q3] }
//   optional screening overrides: pz, zz, p: [p1, p2, p3], e: [e1, e2, e3]
// Z1, Z2 come from the species (z). One pair counted once (MiLaDy adds phi/2 from each atom).
// exaStamp shifts pair energies by V(rcut): to keep MiLaDy's energies use rcut >= r2 (switched),
// or rcut slightly above rr (bridge: MiLaDy's bridge is nonzero up to rr and 0 beyond).

#include <cmath>
#include <string>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <onika/cuda/cuda.h>
#include <onika/flat_tuple.h>

namespace exaStamp
{
  using namespace exanb;

  struct ZBLv2Parms
  {
    int mode = 2;                    // 2 switched, 1 bridge
    double r1 = 0.0;                 // ZBL alone below r1 (angstrom)
    double r2 = 0.0;                 // mode 2: zero above r2
    double rr = 0.0;                 // mode 1: bridge up to rr
    int bridge_type = 1;             // mode 1: 1 exp, 2 poly
    double q[6] = {0., 0., 0., 0., 0., 0.};
    double pz = 0.46848;             // MiLaDy ZBLData defaults
    double zz = 0.23;
    double p[3] = {0.32825, 0.09219, 0.58110};
    double e[3] = {-2.54931, -0.29182, -0.59231};
    double eps0 = 55.26349406e-4;    // e^2 / (eV angstrom)
  };
}

namespace YAML
{
  template<> struct convert<exaStamp::ZBLv2Parms>
  {
    static bool decode(const Node& node, exaStamp::ZBLv2Parms& v)
    {
      using onika::physics::Quantity;
      v = exaStamp::ZBLv2Parms{};
      if( !node.IsMap() || !node["r1"] ) { return false; }
      const std::string mode = node["mode"] ? node["mode"].as<std::string>() : std::string("switched");
      if( mode == "switched" ) v.mode = 2;
      else if( mode == "bridge" ) v.mode = 1;
      else { return false; }
      v.r1 = node["r1"].as<Quantity>().convert();
      if( v.mode == 2 )
      {
        if( !node["r2"] ) { return false; }
        v.r2 = node["r2"].as<Quantity>().convert();
        if( !( v.r2 > v.r1 ) ) { return false; }
      }
      else
      {
        if( !node["rr"] || !node["bridge"] ) { return false; }
        v.rr = node["rr"].as<Quantity>().convert();
        const std::string bt = node["bridge_type"] ? node["bridge_type"].as<std::string>() : std::string("exp");
        if( bt == "exp" ) v.bridge_type = 1;
        else if( bt == "poly" ) v.bridge_type = 2;
        else { return false; }
        const auto q = node["bridge"].as< std::vector<double> >();
        if( q.size() != ( v.bridge_type == 1 ? 4u : 6u ) ) { return false; }
        for( size_t i = 0; i < q.size(); ++i ) v.q[i] = q[i];
      }
      if( node["pz"] ) v.pz = node["pz"].as<double>();
      if( node["zz"] ) v.zz = node["zz"].as<double>();
      if( node["p"] ) { const auto a = node["p"].as< std::vector<double> >(); if( a.size() != 3 ) return false; for( int k = 0; k < 3; ++k ) v.p[k] = a[k]; }
      if( node["e"] ) { const auto a = node["e"].as< std::vector<double> >(); if( a.size() != 3 ) return false; for( int k = 0; k < 3; ++k ) v.e[k] = a[k]; }
      return true;
    }
  };
}

namespace exaStamp
{
  // phi(r) * S(r) between r1 and r2 (S = 1 when r2 <= r1 is not used), energy e and de/dr in eV, eV/angstrom.
  ONIKA_HOST_DEVICE_FUNC inline void zbl_v2_switched(const ZBLv2Parms& p, double z1, double z2, double r, double r1, double r2, double& e, double& de)
  {
    const double pi = 3.14159265358979323846;
    const double cst = z1 * z2 / ( 4.0 * pi * p.eps0 );
    const double a = ( pow( z1, p.zz ) + pow( z2, p.zz ) ) / p.pz;
    const double x = r * a;
    double phi = 0., dphi = 0.;
    for( int k = 0; k < 3; ++k ) { const double t = p.p[k] * exp( p.e[k] * x ); phi += t; dphi += p.e[k] * t; }
    dphi *= a;
    double f = 0., df = 0.;
    if( r <= r1 ) f = 1.;
    else if( r < r2 )
    {
      const double X = ( r - r1 ) / ( r2 - r1 );
      f = 1.0 - X * X * X * ( 6.0 * X * X - 15.0 * X + 10.0 );
      df = ( -30.0 * X * X * X * X + 60.0 * X * X * X - 30.0 * X * X ) / ( r2 - r1 );
    }
    e = cst * phi * f / r;
    de = cst * ( ( phi * df + dphi * f ) / r - phi * f / ( r * r ) );
  }

  ONIKA_HOST_DEVICE_FUNC inline void zbl_v2_compute_energy(const ZBLv2Parms& p, const PairPotentialMinimalParameters& pair, double r, double& e, double& de)
  {
    const double z1 = pair.m_atom_a.m_z, z2 = pair.m_atom_b.m_z;
    e = 0.; de = 0.;
    if( p.mode == 2 ) zbl_v2_switched( p, z1, z2, r, p.r1, p.r2, e, de );
    else if( r <= p.r1 ) zbl_v2_switched( p, z1, z2, r, p.r1, p.r1 + 1.0, e, de );   // S = 1 below r1
    else if( r <= p.rr )
    {
      const double* q = p.q;
      if( p.bridge_type == 1 )
      {
        e = exp( q[0] + q[1] * r + q[2] * r * r + q[3] * r * r * r );
        de = ( q[1] + 2.0 * q[2] * r + 3.0 * q[3] * r * r ) * e;
      }
      else
      {
        e = q[0] + r * ( q[1] + r * ( q[2] + r * ( q[3] + r * ( q[4] + r * q[5] ) ) ) );
        de = q[1] + r * ( 2.0 * q[2] + r * ( 3.0 * q[3] + r * ( 4.0 * q[4] + r * 5.0 * q[5] ) ) );
      }
    }
    static const double conv_energy_inv = 1e-4 * onika::physics::elementaryCharge / onika::physics::atomicMass;   // eV -> internal
    e *= conv_energy_inv;
    de *= conv_energy_inv;
  }
}

#define USTAMP_POTENTIAL_NAME     zbl_v2
#define USTAMP_POTENTIAL_PARAMS   ZBLv2Parms
#define USTAMP_POTENTIAL_COMPUTE  zbl_v2_compute_energy
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) onika::make_flat_tuple(p.m_atom_a.m_z,p.m_atom_b.m_z)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1
