#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>

#include <exanb/core/quantity_yaml.h>
#include <exaStamp/potential_factory/pair_potential.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;


  // Exp6 Parameters
  struct Exp6Parms
  {
    double A = 0.0;
    double B = 1.0;
    double C = 0.0;
    double D = 0.0;
  };

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::Exp6Parms;
  using exanb::Quantity;

  template<> struct convert<Exp6Parms>
  {
    static inline bool decode(const Node& node, Exp6Parms& v)
    {
      if( !node.IsMap() ) { return false; }

      v.A = node["A"].as<Quantity>().convert();
      v.B = node["B"].as<Quantity>().convert();
      v.C = node["C"].as<Quantity>().convert();
      v.D = node["D"].as<Quantity>().convert();

      return true;
    }
  };
}

namespace exaStamp
{
  using namespace exanb;

  ONIKA_HOST_DEVICE_FUNC inline void exp6_compute_energy(const Exp6Parms& p_exp6, const PairPotentialMinimalParameters&, double r, double& e, double& de)
  {
    assert(r>0);
    double one_rB = 1/(r*p_exp6.B);
    double Cr6 = p_exp6.C / pow(r,6);
    double Dtwelve_rB12 = p_exp6.D * pow(12*one_rB,12);
    double AexpmBr = p_exp6.A * exp(-p_exp6.B * r);

    e  =   AexpmBr       - Cr6   + Dtwelve_rB12;
    de = - p_exp6.B * AexpmBr + (6 * Cr6 - 12 * Dtwelve_rB12)/r;
  }

}


