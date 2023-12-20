#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>

#include <exanb/core/quantity_yaml.h>
#include <exanb/core/unityConverterHelper.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exanb/core/physics_constants.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // CoulWolf Parameters
  struct CoulWolfParms
  {
    double alpha = 0.0;
    double rc = 0.0;
    double qqrd2e = 14.399645;
    double e_shift = 0.0;
    double f_shift = 0.0;
  };
  
  // core computation kernel for coul wolf potential
  ONIKA_HOST_DEVICE_FUNC inline void coul_wolf_compute_energy(const CoulWolfParms& p, double c1, double c2, double r, double& e, double& de)
  {
    assert( r > 0. );
    
    // LAMMPS
    double prefactor = p.qqrd2e * c1 * c2 / r;
    double erfcc = erfc(p.alpha * r);
    double erfcd = exp(-p.alpha * p.alpha * r * r);
    double v_sh = (erfcc - p.e_shift * r) * prefactor;
    double dvdrr = (erfcc / ( r * r ) + 2.0 * p.alpha / sqrt(M_PI) * erfcd / r) + p.f_shift;
    double forcecoul = dvdrr * r * r * prefactor;
    double fpair = -forcecoul / r;
    
#   ifdef EXANB_UNITS_V2
    e = EXANB_QUANTITY( v_sh * eV ).convert();
    de = EXANB_QUANTITY( fpair * eV / ang ).convert();
#   else
    e = UnityConverterHelper::convert(v_sh, "eV"); // WARNING : not Cuda compatible
    de = UnityConverterHelper::convert(fpair, "eV/ang");
#   endif
    
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::CoulWolfParms;
  using exanb::UnityConverterHelper;
  using exanb::Quantity;

  template<> struct convert<CoulWolfParms>
  {
    static bool decode(const Node& node, CoulWolfParms& v)
    {
      if( !node.IsMap() ) { return false; }
      v.alpha = node["alpha"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      v.e_shift = erfc(v.alpha * v.rc) / v.rc;
      v.f_shift = -(v.e_shift + 2.0 * v.alpha / sqrt(M_PI) * exp(-v.alpha * v.alpha * v.rc * v.rc)) / v.rc;
      
      return true;
    }
  };
}

