#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>

#include <onika/physics/units.h>
#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <onika/physics/constants.h>

#include <onika/cuda/cuda.h>

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

namespace exaStamp
{
  using namespace exanb;

  // CoulDsf Parameters
  struct CoulDsfParms
  {
    double alpha = 0.0;
    double rc = 0.0;
    double qqrd2e = 14.399645;
  };
  // core computation kernel for coul dsf potential
  ONIKA_HOST_DEVICE_FUNC inline void coul_dsf_compute_energy(const CoulDsfParms& p, double c1, double c2, double r, double& e, double& de)
  {
    assert( r > 0. );
    double MY_PIS = sqrt(M_PI);
    double rsq = r*r;
    double cut_coul = p.rc;
    double cut_coulsq = cut_coul * cut_coul;    
    double qtmp = c1;
    double qj = c2;

    double erfcc = erfc(p.alpha * cut_coul);
    double erfcd = exp(-p.alpha * p.alpha * cut_coul * cut_coul);
    double f_shift = -(erfcc / cut_coulsq + 2.0 / MY_PIS * p.alpha * erfcd / cut_coul);
    double e_shift = erfcc / cut_coul - f_shift * cut_coul;

    //    double e_self = -(e_shift / 2.0 + p.alpha / MY_PIS) * qtmp * qtmp * p.qqrd2e;
    
    double prefactor = p.qqrd2e * qtmp * qj / r;
    erfcd = exp(-p.alpha * p.alpha * rsq);
    double t = 1.0 / (1.0 + EWALD_P * p.alpha * r);
    erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;

    double forcecoul = prefactor * (erfcc / r + 2.0 * p.alpha / MY_PIS * erfcd + r * f_shift) * r;
    double fpair = -forcecoul / r;
    double ecoul = prefactor * (erfcc - r * e_shift - rsq * f_shift);
    
#   ifdef EXANB_UNITS_V2
    e = EXANB_QUANTITY( ecoul * eV ).convert();
    de = EXANB_QUANTITY( fpair * eV / ang ).convert();
#   else
    e = UnityConverterHelper::convert(ecoul, "eV");
    de = UnityConverterHelper::convert(fpair, "eV/ang");;
#   endif

  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::CoulDsfParms;
  using exanb::UnityConverterHelper;
  using onika::physics::Quantity;

  template<> struct convert<CoulDsfParms>
  {
    static bool decode(const Node& node, CoulDsfParms& v)
    {
      if( !node.IsMap() ) { return false; }
      v.alpha = node["alpha"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      return true;
    }
  };
}

