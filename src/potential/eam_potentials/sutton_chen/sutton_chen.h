#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <exanb/core/unityConverterHelper.h>
#include <onika/math/basic_types.h>

namespace exaStamp
{
  using namespace exanb;

  struct EamSuttonChenParameters
  {
    double c;
    double epsilon;
    double a0;
    double n;
    double m;
  };

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  static inline void sutton_chen_phi(const EamSuttonChenParameters& p, double r, double& phiValue, double& dphi)
  {
    double ratio = p.a0/r;
    phiValue  = p.epsilon * std::pow(ratio, p.n);
    dphi = -1 * p.n * phiValue / r;
  }

  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  static inline void sutton_chen_rho(const EamSuttonChenParameters& p, double r, double& rhoValue, double& drho)
  {
      double ratio = p.a0/r;
      rhoValue  = std::pow(ratio, p.m);
      drho = -1 * p.m * rhoValue / r;
  }

  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  static inline void sutton_chen_fEmbed(const EamSuttonChenParameters& p, double rhoValue, double& f, double& df)
  {
      double sqrtRho = std::sqrt(rhoValue);
      f  = -1. * p.c * p.epsilon * sqrtRho;
      df = 0.5 * f / (rhoValue>0? rhoValue : 0);
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EamSuttonChenParameters;
  using exanb::UnityConverterHelper;
  using exanb::Quantity;

  template<> struct convert<EamSuttonChenParameters>
  {
    static bool decode(const Node& node, EamSuttonChenParameters& v)
    {
      if( !node.IsMap() ) { return false; }
      v.c = node["c"].as<Quantity>().convert();
      v.epsilon = node["epsilon"].as<Quantity>().convert();
      v.a0 = node["a0"].as<Quantity>().convert();
      v.n = node["n"].as<Quantity>().convert();
      v.m = node["m"].as<Quantity>().convert();
      return true;
    }
  };
}

