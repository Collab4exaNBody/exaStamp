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
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_yaml.h>
#include <cmath>

#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  namespace ewald_constants
  {
    static constexpr double epsilonZero = 0x1.337f14f782bb5p-21; //0x1.2c4a1a79b5aafp-11 * 0.001;
    static constexpr double fpe0 = 0x1.e303a50b0d1ap-18; // 0x1.d7b18f2ccacb6p-8 * 0.001 ;
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
    double Gc;
  };
  
  struct alignas(DEFAULT_ALIGNMENT) EwaldParameters
  {
    // Approximation of erfc() function
    double EWALD_F = 1.12837917;
    double EWALD_P = 0.3275911;
    double A1 = 0.254829592;
    double A2 = -0.284496736;
    double A3 = 1.421413741;
    double A4 = -1.453152027;
    double A5 = 1.061405429;
    
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
    double gm_sr = 0.0; // the gm used in short range computation
    double qqr2e = 0.0;
    double bt_sr = 0.0; // the bt value used in short range computation
    
    double volume = 0.0;
    Vec3d unitk = { 0.0 , 0.0 , 0.0 };
    
    onika::memory::CudaMMVector<EwaldCoeffs> Gdata;
  };

  struct ReadOnlyEwaldParameters
  {
    double g_ewald = 0.0;
    ssize_t nknz = 0;
    double gm_sr = 0.0; // the gm used in short range computation
    double bt_sr = 0.0; // the bt value used in short range computation
    
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
      , Gdata( p.Gdata.data() )
    {}
  };

  template<class EwaldParametersT>
  ONIKA_HOST_DEVICE_FUNC static inline void ewald_compute_energy(const EwaldParametersT& p, double c, double r, double& e, double& de)
  {
    const double cf = p.gm_sr * c;
    const double grij = p.g_ewald * r;
    const double expm2 = exp(-grij * grij);
    const double t = 1.0 / (1.0 + p.EWALD_P * grij);
    const double erfc = t * (p.A1 + t * (p.A2 + t * (p.A3 + t * (p.A4 + t * p.A5)))) * expm2;
    e = cf * erfc / r;
    de = - (cf * p.bt_sr * expm2 + e) / r;
  }

  inline double ewald_error_accuracy(double g_ewald, int km, double length, uint64_t natoms, double q2){
    if (natoms == 0) natoms = 1;
    double value = 2.0*q2*g_ewald/length * sqrt(1.0/(M_PI*km*natoms)) * std::exp(-M_PI*M_PI*km*km/(g_ewald*g_ewald*length*length));
    return value;
  }
  
  inline void ewald_init_parameters(double g_ewald, double radius, double accuracy_relative, unsigned int in_kmax, const Vec3d& domainSize, const uint64_t natoms, double qsq, double qsum, EwaldParameters& p , ::exanb::LogStreamWrapper& ldbg )
  {
    
    using ewald_constants::fpe0;
    using ewald_constants::epsilonZero;
    
    p.g_ewald = g_ewald;
    p.radius = radius;
    p.accuracy_relative = accuracy_relative;
    double accuracy = accuracy_relative;
    p.kmax = in_kmax;

    p.volume = domainSize.x * domainSize.y * domainSize.z ;

    double xL = domainSize.x;
    double yL = domainSize.y;
    double zL = domainSize.z;
    
    if( p.volume == 0.0 ) return;

    // ------------------------------------------------------------------- //
    // 1st step : g_ewald calculation
    if(p.g_ewald <= 0.)
    {
      double g_ewald = 0.0;
      g_ewald = accuracy_relative*sqrt(natoms*radius*xL*yL*zL) / (2.0*qsq);
      if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*std::log(accuracy))/radius;
      else g_ewald = sqrt(-std::log(g_ewald)) / radius;
      p.g_ewald = g_ewald;
    }
    
    if(p.g_ewald == 0.)
    {
        lerr<<"######## ERROR in ewald_init_parameters : g_ewald ="<<p.g_ewald<<" - Decrease accuracy_relative of Ewald method"<<std::endl<<std::flush;
        std::abort();
    }
    // ------------------------------------------------------------------- //

    // ------------------------------------------------------------------- //
    // 2nd step : kmax calculation
    if(p.kmax <= 0)
    {
        if(p.g_ewald == 0.)
        {
          lerr<<"######## ERROR in ewald_init_parameters : g_ewald ="<<p.g_ewald<<" in kmax calculation"<<std::endl<<std::flush;
          std::abort();
        }
        
        double err;
        int kmax = 0;
        int kxmax = 1;
        int kymax = 1;
        int kzmax = 1;

        err = ewald_error_accuracy(p.g_ewald,kxmax,xL,natoms,qsq);
        while (err > accuracy) {
          kxmax++;
          err = ewald_error_accuracy(p.g_ewald,kxmax,xL,natoms,qsq);
        }
      
        err = ewald_error_accuracy(p.g_ewald,kymax,yL,natoms,qsq);
        while (err > accuracy) {
          kymax++;
          err = ewald_error_accuracy(p.g_ewald,kymax,yL,natoms,qsq);
        }
      
        err = ewald_error_accuracy(p.g_ewald,kzmax,zL,natoms,qsq);
        while (err > accuracy) {
          kzmax++;
          err = ewald_error_accuracy(p.g_ewald,kzmax,zL,natoms,qsq);
        }
      
        kmax = std::max(kxmax,kymax);
        kmax = std::max(kmax,kzmax);
        p.kmax = kmax;
        p.kxmax = kxmax;
        p.kymax = kymax;
        p.kzmax = kzmax;
    }
    
    if(p.kmax < 2)
      {
        lout<<"######## ERROR in ewald_init_parameters : kmax ="<<p.kmax<<" - Decrease epsilon of Ewald method"<<std::endl;
        std::abort();
      }
    // ------------------------------------------------------------------- //

    p.unitk = Vec3d{ 2.*M_PI/xL , 2.*M_PI/yL , 2.*M_PI/zL };
    double GnMax_x = (p.unitk).x * (p.unitk).x * p.kxmax * p.kxmax;
    double GnMax_y = (p.unitk).y * (p.unitk).y * p.kymax * p.kymax;
    double GnMax_z = (p.unitk).z * (p.unitk).z * p.kzmax * p.kzmax;
    double GnMax = std::max(GnMax_x , GnMax_y);
    GnMax = std::max(GnMax , GnMax_z);
    
    p.nk = (2 * p.kmax + 1) * (2 * p.kmax + 1) * (2 * p.kmax + 1) - 1;
    p.Gdata.resize( p.nk );

    double bt = 2. * M_PI / fpe0 / p.volume;
    p.bt_sr = 2. * p.g_ewald / std::sqrt(M_PI);
    p.gm = 1. / (4. * p.g_ewald * p.g_ewald);
    p.gm_sr = 0.25 / M_PI / epsilonZero;
    p.qqr2e = 0.25 / M_PI / epsilonZero;

    // note : this triple-loop is parallelizable
    const ssize_t kxmin = - static_cast<ssize_t>(p.kxmax);
    const ssize_t kxmax =   static_cast<ssize_t>(p.kxmax);
    const ssize_t kymin = - static_cast<ssize_t>(p.kymax);
    const ssize_t kymax =   static_cast<ssize_t>(p.kymax);
    const ssize_t kzmin = - static_cast<ssize_t>(p.kzmax);
    const ssize_t kzmax =   static_cast<ssize_t>(p.kzmax);
    
    size_t kk = 0;
    double gcmin = 1e30;
    double gcmax = 0.0;

    for (ssize_t kx=kxmin; kx<=kxmax; ++kx )
    {
      for (ssize_t ky=kymin; ky<=kymax; ++ky)
      {
        for (ssize_t kz=kzmin; kz<=kzmax; ++kz)
        {
	        const ssize_t nk = kx*kx + ky*ky + kz*kz;
	        if( nk > 0 )
          {
            const Vec3d G_kk = { kx * p.unitk.x, ky * p.unitk.y, kz * p.unitk.z };
            const double Gn_kk = norm2(G_kk);
            if ( Gn_kk <= GnMax)
            {
              assert( kk < static_cast<size_t>(p.nk) );
              double Gc_kk = std::exp(-p.gm * Gn_kk ) / Gn_kk;
              if(Gc_kk<gcmin) gcmin=Gc_kk;
              if(Gc_kk>gcmax) gcmax=Gc_kk;
              ldbg<<"G^2="<<Gn_kk<<" exp(-G^2/(4 g_ewald^2)/G^2)="<<Gc_kk<<std::endl<<std::flush;
              Gc_kk *= bt;
              p.Gdata[kk] = EwaldCoeffs{ G_kk.x , G_kk.y , G_kk.z , Gc_kk };
              ++kk;
            }
          }
        }
      }
    }

    ldbg<<"   exp(-G^2/4g_ewald^2)/G^2 : minimum value :"<<gcmin<<std::endl<<std::flush;
    ldbg<<"                          : maximum value :"<<gcmax<<std::endl<<std::flush;

    // number of non zero values
    p.nknz = kk;    
    ldbg << "   number of k points="<< p.nknz       <<std::endl<<std::flush;

    // adjust coeffs size
    p.Gdata.resize( p.nknz );
    p.Gdata.shrink_to_fit();
  }

  inline void ewald_init_parameters(double g_ewald, double radius, double accuracy_relative, unsigned int in_kmax, const Vec3d& domainSize, const uint64_t natoms, double qsq, double qsum, EwaldParameters& p )
  {
    ewald_init_parameters(g_ewald,radius,accuracy_relative,in_kmax,domainSize,natoms,qsq,qsum,p , ::exanb::ldbg<<"" );
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
        node["kmax"].as<long>(), domSize, 1, 1.,1.,v);
      return true;
    }
  };
}
