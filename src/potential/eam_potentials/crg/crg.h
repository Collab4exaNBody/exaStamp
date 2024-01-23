#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>
#include <exanb/core/unityConverterHelper.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/physics_constants.h>
#include <onika/cuda/cuda.h>

#include <exaStamp/potential/eam/eam_buffer.h>

namespace exaStamp
{
  using namespace exanb;
  using exanb::UnityConverterHelper;

  struct EamCRGParameters
  {
    //Pairwise
    double A      ;
    double rho_buck  ;
    double C      ;
    double D      ;
    double gamma  ;
    double r0     ;
    // EAM
    double G;
    double n;    
    // double G_a;
    // double n_a;
    // double G_b;
    // double n_b;
  };

  // Type to return the phi and dphi values of the pairwise terms
  typedef struct {
    double phi_pair;
    double dphi_pair;
  } phi_eval;

  // /// @brief Calculate the Coulomb term of the energy and its derivative
  // /// @param [in] r Interatomic distance
  // static inline void crg_coulomb(const EamCRGParameters& p, double r, phi_eval &phi)
  // {
  //   double qa=0.;
  //   double qb=0.;
  //   double dielectric = 1.;
  //   double qqr2e = 14.399645;
  //   double qqrd2e = UnityConverterHelper::convert(qqr2e, "eV.ang/e-^2") / dielectric;
  //   double C = qqrd2e;
  //   phi.phi_pair   =   C*qa*qb/r;
  //   phi.dphi_pair  =  -C*qa*qb/std::pow(r,2);
  // }
  // /// @brief Calculate the Coulomb term of the energy and its derivative
  // /// @param [in] r Interatomic distance
  // static inline void crg_coulomb_mm(const EamCRGParameters& p, double r, phi_eval &phi, double qa, double qb)
  // {
  //   const double FourPiEpsilon0 = 1.;//*UnityConverterHelper::convert(exanb::legacy_constant::epsilonZero, "C^2.s^2/m^3/kg^1");
  //   double dielectric = 1.;
  //   double qqr2e = 14.399645; // units = eV.ang/e-^2
  //   double qqrd2e = UnityConverterHelper::convert(qqr2e, "eV.ang/e-^2") / dielectric;    
  //   double C = qqrd2e;
  //   phi.phi_pair   =   C*qa*qb/r;
  //   phi.dphi_pair  =  -C*qa*qb/std::pow(r,2);
  // }
  
  /// @brief Calculate the Morse term of the energy and its derivative
  /// @param [in] r Interatomic distance
  ONIKA_HOST_DEVICE_FUNC
  static inline void crg_morse(const EamCRGParameters& p, double r, phi_eval &phi)
  {
    double k  = p.gamma * (r - p.r0);
    double dk = p.gamma;
    phi.phi_pair  = p.D * ( exp( -2. * k) - 2. * exp( -k ) );
    phi.dphi_pair = -2. * dk * p.D * ( exp( -2. * k) - exp( -k ) );
  }

  /// @brief Calculate the Buckingham term of the energy and its derivative
  /// @param [in] r Interatomic distance
  ONIKA_HOST_DEVICE_FUNC
  static inline void crg_buckingham(const EamCRGParameters& p, double r, phi_eval &phi)
  {
    double puis = -r / p.rho_buck;
    double dpuis = - 1. / p.rho_buck;
    double sec = - p.C / pow(r,6);
    double dsec = - 6. * sec / r;
    
    phi.phi_pair  = p.A * exp( puis ) + sec;
    phi.dphi_pair = p.A * dpuis * exp( puis ) + dsec;
    // phi.phi_pair  = p.A * std::exp( -r / p.rho_buck ) - (p.C / std::pow(r,6));
    // phi.dphi_pair = ( - p.A / p.rho_buck )  * std::exp( -r / p.rho_buck ) + 6. *(p.C / std::pow(r,7));
  }

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative, mono-material version
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) = \phi_c(r) + \phi_m(r) + \phi_b(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC
  static inline void crg_phi(const EamCRGParameters& p, double r, double& phiValue, double& dphi)
  {
    phi_eval phi_m ;
    phi_eval phi_b ;
    
    crg_morse(p,r,phi_m);
    crg_buckingham(p,r,phi_b);

    phiValue  = phi_m.phi_pair  + phi_b.phi_pair ;
    dphi      = phi_m.dphi_pair + phi_b.dphi_pair ;
  }

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative, multimaterial version
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) = \phi_c(r) + \phi_m(r) + \phi_b(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC
  static inline void crg_phi_mm(const EamCRGParameters& p, double r, double& phiValue, double& dphi, const EAMSpecyPairInfo& pair_info , bool pair_inversed = false )
  {
    // may use pair_info.m_charge_a , pair_info.m_charge_b
    
    phi_eval phi_m ;
    phi_eval phi_b ;
    
    crg_morse(p,r,phi_m);
    crg_buckingham(p,r,phi_b);

    phiValue  = phi_m.phi_pair  + phi_b.phi_pair;
    dphi      = phi_m.dphi_pair + phi_b.dphi_pair ;
  }

  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative, mono-material version
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC
  static inline void crg_rho(const EamCRGParameters& p, double r, double& rhoValue, double& drho, bool pair_inversed = false)
  {
    //    double n = pair_inversed ? p.n_b : p.n_a;    
    //    double sigma=n/std::pow(r,8);
    double sigma=p.n/pow(r,8);
    double ur = 20 * ( r - 1.5);
    rhoValue  = sigma * 0.5 * ( 1 + erf( ur ) );
    drho = -8. * rhoValue / r + sigma * 20. * exp(-pow(ur,2)) / sqrt(M_PI);
  }

  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative, multi-material version
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC
  static inline void crg_rho_mm(const EamCRGParameters& p, double r, double& rhoValue, double& drho, const EAMSpecyPairInfo& pair_info , bool pair_inversed = false )
  {
    // si on duplique n en n_a et n_b pour les atomes de la paire (a,b), alors :
    // pair_inversed=false => p.n_a = atome central, p.n_b = atome peripherique
    // pair_inversed=true => p.n_b = atome central, p.n_a = atome peripherique
    //    double n = pair_inversed ? p.n_a : p.n_b;
    
    // pour la charge
    // pair_inversed=false => pair_info.m_charge_a = atome central, pair_info.m_charge_b = atome peripherique
    // pair_inversed=true => pair_info.m_charge_b = atome central, pair_info.m_charge_a = atome peripherique
    
    // may use pair_info.m_charge_a , pair_info.m_charge_b

    double sigma=p.n/pow(r,8);
    double ur = 20 * ( r - 1.5);
    rhoValue = sigma * 0.5 * ( 1 + erf( ur ) );
    drho = -8. * sigma * 0.5 * ( 1 + erf( ur ) ) / r + 0.5 * sigma * 20. * exp(-pow(ur,2)) / sqrt(M_PI);

  }

  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  ONIKA_HOST_DEVICE_FUNC
  static inline void crg_fEmbed(const EamCRGParameters& p, double rhoValue, double& f, double& df , bool pair_inversed = false )
  {
    double sqrtRho = sqrt(rhoValue);
    double G = p.G;
    //      f  = -1. * p.G * sqrtRho;
    f  = -1. * G * sqrtRho;      
    if (sqrtRho>0.) {
      df = 0.5 * f / rhoValue;
    } else {
      df = 0.;
    }
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EamCRGParameters;
  using exanb::UnityConverterHelper;
  using exanb::Quantity;

  template<> struct convert<EamCRGParameters>
  {
    static bool decode(const Node& node, EamCRGParameters& v)
    {
      if( !node.IsMap() ) { return false; }
      v.A         = node["A"].as<Quantity>().convert();
      v.rho_buck  = node["rho"].as<Quantity>().convert();
      v.C         = node["C"].as<Quantity>().convert();
      v.D         = node["D"].as<Quantity>().convert();
      v.gamma     = node["gamma"].as<Quantity>().convert();
      v.r0        = node["r0"].as<Quantity>().convert();
      v.G        = node["G"].as<Quantity>().convert();
      v.n        = node["n"].as<Quantity>().convert();
      // v.G_a        = node["G_a"].as<Quantity>().convert();
      // v.n_a        = node["n_a"].as<Quantity>().convert();
      // v.G_b        = node["G_b"].as<Quantity>().convert();
      // v.n_b        = node["n_b"].as<Quantity>().convert();            
      //      std::cout<<"EamCRGParameters tout Ok\n";
      return true;
    }
  };
}

