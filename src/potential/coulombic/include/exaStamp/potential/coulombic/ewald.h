/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/


#pragma once

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/log.h>
#include <exaStamp/unit_system.h>
#include <exaStamp/coulomb_constant.h>
#include <cmath>
#include <algorithm>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
inline namespace coulombic_ewald // distinct symbols from the legacy ewald plugin (plugins are loaded RTLD_GLOBAL)
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  namespace ewald_constants
  {
    // Coulomb constant 1/(4.pi.epsilon0) in internal units (LAMMPS metal units value, see exaStamp/coulomb_constant.h)
    static constexpr double qqr2e = COULOMB_CONSTANT;
    static constexpr double fpe0 = 1.0 / qqr2e;                 // 4.pi.epsilon0
    static constexpr double epsilonZero = fpe0 / ( 4.0 * M_PI ); // epsilon0

    // Abramowitz & Stegun approximation of erfc(x) (same as LAMMPS pair coul/long without tables)
    static constexpr double EWALD_P = 0.3275911;
    static constexpr double A1 = 0.254829592;
    static constexpr double A2 = -0.284496736;
    static constexpr double A3 = 1.421413741;
    static constexpr double A4 = -1.453152027;
    static constexpr double A5 = 1.061405429;
  }

  struct EwaldRho
  {
    size_t nk = 0;
    onika::memory::CudaMMVector<Complexd> rho;
  };

  struct EwaldCoeffs
  {
    double Gx;
    double Gy;
    double Gz;
    double Gc; // 2.pi/(4.pi.epsilon0.V) exp(-G^2/(4g^2))/G^2
    double Gv; // 2.(1+G^2/(4g^2))/G^2 , used for reciprocal virial
  };
  
  struct alignas(DEFAULT_ALIGNMENT) EwaldParameters
  {
    double g_ewald = 0.0;
    double radius = 0.0; 
    double accuracy_relative = 0.0;

    ssize_t kmax = 0;
    ssize_t kxmax = 0;
    ssize_t kymax = 0;
    ssize_t kzmax = 0;
    ssize_t nk = 0;
    ssize_t nknz = 0;

    double gm = 0.0;
    double gm_sr = 0.0; // 1/(4.pi.epsilon0), used in short range computation
    double qqr2e = 0.0;
    double bt_sr = 0.0; // 2.g/sqrt(pi), used in short range computation
    double qsum = 0.0;  // total charge, for neutralizing background energy
    
    double volume = 0.0;
    Vec3d box = { 0.0 , 0.0 , 0.0 }; // box size used to build k vectors
    Vec3d unitk = { 0.0 , 0.0 , 0.0 };
    
    onika::memory::CudaMMVector<EwaldCoeffs> Gdata;
  };

  // trivially copyable view of EwaldParameters, suitable for GPU functors
  struct ReadOnlyEwaldParameters
  {
    double g_ewald = 0.0;
    ssize_t nknz = 0;
    double gm_sr = 0.0;
    double bt_sr = 0.0;
    double qsum = 0.0;
    double volume = 0.0;
    
    const EwaldCoeffs* __restrict__ Gdata = nullptr;
    
    ReadOnlyEwaldParameters() = default;
    ReadOnlyEwaldParameters(const ReadOnlyEwaldParameters&) = default;
    ReadOnlyEwaldParameters(ReadOnlyEwaldParameters&&) = default;
    ReadOnlyEwaldParameters& operator = (const ReadOnlyEwaldParameters&) = default;
    ReadOnlyEwaldParameters& operator = (ReadOnlyEwaldParameters&&) = default;
    
    inline ReadOnlyEwaldParameters( const EwaldParameters & p )
      : g_ewald( p.g_ewald )
      , nknz( p.nknz )
      , gm_sr( p.gm_sr )
      , bt_sr( p.bt_sr )
      , qsum( p.qsum )
      , volume( p.volume )
      , Gdata( p.Gdata.data() )
    {}
  };

  // real space part : e = qi.qj/(4.pi.epsilon0) erfc(g.r)/r , de = de/dr
  template<class EwaldParametersT>
  ONIKA_HOST_DEVICE_FUNC static inline void ewald_compute_energy(const EwaldParametersT& p, double c, double r, double& e, double& de)
  {
    using namespace ewald_constants;
    const double cf = p.gm_sr * c;
    const double grij = p.g_ewald * r;
    const double expm2 = exp(-grij * grij);
    const double t = 1.0 / (1.0 + EWALD_P * grij);
    const double erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
    e = cf * erfc / r;
    de = - (cf * p.bt_sr * expm2 + e) / r;
  }

  // self energy + neutralizing background energy of one particle with charge q
  template<class EwaldParametersT>
  ONIKA_HOST_DEVICE_FUNC static inline double ewald_self_energy(const EwaldParametersT& p, double q)
  {
    using namespace ewald_constants;
    return - qqr2e * ( p.g_ewald / sqrt(M_PI) * q * q + 0.5 * M_PI * q * p.qsum / ( p.g_ewald * p.g_ewald * p.volume ) );
  }

  // rms force error estimate of the reciprocal part (same as LAMMPS Ewald::rms), q2 = sum of squared charges
  inline double ewald_error_accuracy(double g_ewald, int km, double length, uint64_t natoms, double q2)
  {
    if (natoms == 0) natoms = 1;
    double value = 2.0*q2*g_ewald/length * sqrt(1.0/(M_PI*km*natoms)) * std::exp(-M_PI*M_PI*km*km/(g_ewald*g_ewald*length*length));
    return value;
  }
  
  inline void ewald_init_parameters(double g_ewald, double radius, double accuracy_relative, long in_kmax, const Vec3d& domainSize, const uint64_t natoms, double qsq, double qsum, EwaldParameters& p , ::exanb::LogStreamWrapper& ldbg )
  {
    using ewald_constants::fpe0;
    
    p.g_ewald = g_ewald;
    p.radius = radius;
    p.accuracy_relative = accuracy_relative;
    p.qsum = qsum;
    p.kmax = in_kmax;
    p.kxmax = p.kymax = p.kzmax = 0;
    p.box = domainSize;
    p.volume = domainSize.x * domainSize.y * domainSize.z ;

    const double xL = domainSize.x;
    const double yL = domainSize.y;
    const double zL = domainSize.z;
    
    if( p.volume == 0.0 ) return;

    // ------------------------------------------------------------------- //
    // 1st step : g_ewald calculation (LAMMPS Ewald::init). accuracy_relative is relative to the
    // force between two unit charges at 1 ang ; it cancels out with qsq except in the log() branch below
    const double accuracy = accuracy_relative;
    if(p.g_ewald <= 0.)
    {
      double g = accuracy_relative*sqrt(natoms*radius*xL*yL*zL) / (2.0*qsq);
      const double accuracy_abs = accuracy_relative * COULOMB_CONSTANT_EV_ANG; // eV/ang, as in LAMMPS metal units
      if (g >= 1.0) g = (1.35 - 0.15*std::log(accuracy_abs))/radius;
      else g = sqrt(-std::log(g)) / radius;
      p.g_ewald = g;
    }
    
    if( ! ( p.g_ewald > 0. ) )
    {
      ::onika::fatal_error() << "ewald_init_parameters : g_ewald=" << p.g_ewald << " - Decrease accuracy_relative of Ewald method" << std::endl;
    }
    // ------------------------------------------------------------------- //

    // ------------------------------------------------------------------- //
    // 2nd step : kmax calculation
    if(p.kmax <= 0)
    {
      p.kxmax = 1;
      while( ewald_error_accuracy(p.g_ewald,p.kxmax,xL,natoms,qsq) > accuracy ) ++ p.kxmax;
      p.kymax = 1;
      while( ewald_error_accuracy(p.g_ewald,p.kymax,yL,natoms,qsq) > accuracy ) ++ p.kymax;
      p.kzmax = 1;
      while( ewald_error_accuracy(p.g_ewald,p.kzmax,zL,natoms,qsq) > accuracy ) ++ p.kzmax;
      p.kmax = std::max( p.kxmax , std::max( p.kymax , p.kzmax ) );
    }
    else
    {
      // user defined kmax, same in all directions (LAMMPS kmax/ewald kmax kmax kmax)
      p.kxmax = p.kymax = p.kzmax = p.kmax;
    }
    
    if(p.kmax < 2)
    {
      ::onika::fatal_error() << "ewald_init_parameters : kmax=" << p.kmax << " - Decrease accuracy_relative of Ewald method" << std::endl;
    }
    // ------------------------------------------------------------------- //

    p.unitk = Vec3d{ 2.*M_PI/xL , 2.*M_PI/yL , 2.*M_PI/zL };
    const double GnMax_x = p.unitk.x * p.unitk.x * p.kxmax * p.kxmax;
    const double GnMax_y = p.unitk.y * p.unitk.y * p.kymax * p.kymax;
    const double GnMax_z = p.unitk.z * p.unitk.z * p.kzmax * p.kzmax;
    // 1.00001 margin as in LAMMPS, so that k vectors exactly on the sphere are not lost to rounding
    const double GnMax = std::max( GnMax_x , std::max( GnMax_y , GnMax_z ) ) * 1.00001;
    
    p.nk = (2 * p.kxmax + 1) * (2 * p.kymax + 1) * (2 * p.kzmax + 1) - 1;
    p.Gdata.resize( p.nk );

    const double bt = 2. * M_PI / fpe0 / p.volume;
    p.bt_sr = 2. * p.g_ewald / std::sqrt(M_PI);
    p.gm = 1. / (4. * p.g_ewald * p.g_ewald);
    p.gm_sr = ewald_constants::qqr2e;
    p.qqr2e = ewald_constants::qqr2e;

    size_t kk = 0;
    double gcmin = 1e30;
    double gcmax = 0.0;

    for (ssize_t kx=-p.kxmax; kx<=p.kxmax; ++kx )
    {
      for (ssize_t ky=-p.kymax; ky<=p.kymax; ++ky)
      {
        for (ssize_t kz=-p.kzmax; kz<=p.kzmax; ++kz)
        {
          if( kx*kx + ky*ky + kz*kz > 0 )
          {
            const Vec3d G_kk = { kx * p.unitk.x, ky * p.unitk.y, kz * p.unitk.z };
            const double Gn_kk = norm2(G_kk);
            if ( Gn_kk <= GnMax)
            {
              assert( kk < static_cast<size_t>(p.nk) );
              double Gc_kk = std::exp(-p.gm * Gn_kk ) / Gn_kk;
              gcmin = std::min( gcmin , Gc_kk );
              gcmax = std::max( gcmax , Gc_kk );
              Gc_kk *= bt;
              const double Gv_kk = 2.0 * ( 1.0 + Gn_kk * p.gm ) / Gn_kk;
              p.Gdata[kk] = EwaldCoeffs{ G_kk.x , G_kk.y , G_kk.z , Gc_kk , Gv_kk };
              ++kk;
            }
          }
        }
      }
    }

    ldbg<<"   exp(-G^2/4g_ewald^2)/G^2 : minimum value :"<<gcmin<<std::endl;
    ldbg<<"                            : maximum value :"<<gcmax<<std::endl;

    // number of non zero values
    p.nknz = kk;    
    ldbg << "   number of k points="<< p.nknz <<std::endl;

    // adjust coeffs size
    p.Gdata.resize( p.nknz );
    p.Gdata.shrink_to_fit();
  }

  inline void ewald_init_parameters(double g_ewald, double radius, double accuracy_relative, long in_kmax, const Vec3d& domainSize, const uint64_t natoms, double qsq, double qsum, EwaldParameters& p )
  {
    ewald_init_parameters(g_ewald,radius,accuracy_relative,in_kmax,domainSize,natoms,qsq,qsum,p , ::exanb::ldbg<<"" );
  }

}
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EwaldParameters;
  
  using onika::physics::Quantity;
  using exanb::Vec3d;
  using exaStamp::ewald_init_parameters;

  template<> struct convert<EwaldParameters>
  {
    static bool decode(const Node& node, EwaldParameters& v)
    {
      if( !node.IsMap() ) { return false; }
      double g_ewald = 0.0;
      if( node["g_ewald"] )
      {
        g_ewald = node["g_ewald"].as<Quantity>().convert();
      }
      Vec3d domSize = node["size"].as<Vec3d>();
      ewald_init_parameters(
        g_ewald,
        node["radius"].as<Quantity>().convert(),
        node["accuracy_relative"].as<Quantity>().convert(),
        node["kmax"].as<long>(), domSize, 1, 1.,0.,v);
      return true;
    }
  };
}
