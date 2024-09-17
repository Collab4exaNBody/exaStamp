#pragma once

#include <cmath>
#include <limits>
#include <utility>
#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exanb/core/physics_constants.h>

#include <onika/cuda/cuda.h>
#include <onika/flat_tuple.h>

namespace exaStamp
{
using namespace exanb;

// ZBL Parameters
struct ZBLParms
{
  // static constexpr double e = 1.6e-19;
  // static constexpr double ct1 = 9.e9; // 1 / (4*M_PI*e0)

  // static constexpr double a1   = 0.4685;
  // static constexpr double a2   = 0.23;
  // static constexpr double phi1 = 0.18175;
  // static constexpr double phi2 = 3.19980;
  // static constexpr double phi3 = 0.50986;
  // static constexpr double phi4 = 0.94229;
  // static constexpr double phi5 = 0.28022;
  // static constexpr double phi6 = 0.40290;
  // static constexpr double phi7 = 0.02817;
  // static constexpr double phi8 = 0.20162;

  double r1;
  double rc;

  static constexpr double pzbl = 0.23;
  static constexpr double a0 = 0.46850;
  static constexpr double c1 = 0.02817;
  static constexpr double c2 = 0.28022;
  static constexpr double c3 = 0.50986;
  static constexpr double c4 = 0.18175;
  static constexpr double d1 = 0.20162;
  static constexpr double d2 = 0.40290;
  static constexpr double d3 = 0.94229;
  static constexpr double d4 = 3.19980;
  
};
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::ZBLParms>
  {
    static bool decode(const Node& node, exaStamp::ZBLParms& v)
    {
      using exanb::Quantity;
      v = exaStamp::ZBLParms{};
      if( !node.IsMap() ) { return false; }
      if( ! node["r1"] ) { return false; }
      if( ! node["rc"] ) { return false; }
      v.r1 = node["r1"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      return true;
    }
  };
}

namespace exaStamp
{
using namespace exanb;

ONIKA_HOST_DEVICE_FUNC inline double e_zbl(const ZBLParms& p, double r, double d1a, double d2a, double d3a, double d4a, double zze)
{

  double d1aij = d1a;
  double d2aij = d2a;
  double d3aij = d3a;
  double d4aij = d4a;
  double zzeij = zze;
  double rinv = 1.0/r;

  double sum = p.c1*exp(-d1aij*r);
  sum += p.c2*exp(-d2aij*r);
  sum += p.c3*exp(-d3aij*r);
  sum += p.c4*exp(-d4aij*r);

  double result = zzeij*sum*rinv;

  return result;  

}

ONIKA_HOST_DEVICE_FUNC inline double dzbldr(const ZBLParms& p, double r, double d1a, double d2a, double d3a, double d4a, double zze)
{

  double d1aij = d1a;
  double d2aij = d2a;
  double d3aij = d3a;
  double d4aij = d4a;
  double zzeij = zze;
  double rinv = 1.0/r;

  double e1 = exp(-d1aij*r);
  double e2 = exp(-d2aij*r);
  double e3 = exp(-d3aij*r);
  double e4 = exp(-d4aij*r);

  double sum = p.c1*e1;
  sum += p.c2*e2;
  sum += p.c3*e3;
  sum += p.c4*e4;

  double sum_p = -p.c1*d1aij*e1;
  sum_p -= p.c2*d2aij*e2;
  sum_p -= p.c3*d3aij*e3;
  sum_p -= p.c4*d4aij*e4;

  double result = zzeij*(sum_p - sum*rinv)*rinv;

  return result;
}

ONIKA_HOST_DEVICE_FUNC inline double d2zbldr2(const ZBLParms& p, double r, double d1a, double d2a, double d3a, double d4a, double zze) {

  double d1aij = d1a;
  double d2aij = d2a;
  double d3aij = d3a;
  double d4aij = d4a;
  double zzeij = zze;
  double rinv = 1.0/r;

  double e1 = exp(-d1aij*r);
  double e2 = exp(-d2aij*r);
  double e3 = exp(-d3aij*r);
  double e4 = exp(-d4aij*r);

  double sum = p.c1*e1;
  sum += p.c2*e2;
  sum += p.c3*e3;
  sum += p.c4*e4;

  double sum_p = p.c1*e1*d1aij;
  sum_p += p.c2*e2*d2aij;
  sum_p += p.c3*e3*d3aij;
  sum_p += p.c4*e4*d4aij;

  double sum_pp = p.c1*e1*d1aij*d1aij;
  sum_pp += p.c2*e2*d2aij*d2aij;
  sum_pp += p.c3*e3*d3aij*d3aij;
  sum_pp += p.c4*e4*d4aij*d4aij;

  double result = zzeij*(sum_pp + 2.0*sum_p*rinv +
                         2.0*sum*rinv*rinv)*rinv;

  return result;
}

ONIKA_HOST_DEVICE_FUNC inline void zbl_compute_energy(const ZBLParms& p, const PairPotentialMinimalParameters& pair, double rij, double& e, double& de)
{
  assert( rij > 0. );
  
  double qqr2e = 14.399645;
  double qelectron = 1.0;
  
  double ainv = ( pow( pair.m_atom_a.m_z , p.pzbl ) + pow( pair.m_atom_b.m_z , p.pzbl ) ) / p.a0;
  
  double d1a = p.d1*ainv;
  double d2a = p.d2*ainv;
  double d3a = p.d3*ainv;
  double d4a = p.d4*ainv;
  double zze = pair.m_atom_a.m_z * pair.m_atom_b.m_z * qqr2e * qelectron * qelectron;

  // std::cout << "qqr2e     = " << qqr2e << std::endl;
  // std::cout << "qelectron = " << qelectron << std::endl;
  // std::cout << "ainv      = " << ainv << std::endl;
  // std::cout << "d1a       = " << d1a << std::endl;
  // std::cout << "d2a       = " << d2a << std::endl;
  // std::cout << "d3a       = " << d3a << std::endl;
  // std::cout << "d4a       = " << d4a << std::endl;
  // std::cout << "zze       = " << zze << std::endl << std::endl;
  
  // e =  t^3 (sw3 + sw4*t) + sw5
  //   = A/3*t^3 + B/4*t^4 + C
  // sw3 = A/3
  // sw4 = B/4
  // sw5 = C

  // dedr = t^2 (sw1 + sw2*t)
  //      = A*t^2 + B*t^3
  // sw1 = A
  // sw2 = B

  // de2dr2 = 2*A*t + 3*B*t^2

  // Require that at t = tc:
  // e = -Fc
  // dedr = -Fc'
  // d2edr2 = -Fc''

  // Hence:
  // A = (-3Fc' + tc*Fc'')/tc^2
  // B = ( 2Fc' - tc*Fc'')/tc^3
  // C = -Fc + tc/2*Fc' - tc^2/12*Fc''
  double cut_global=p.rc;
  double cut_inner=p.r1;
  double tc = cut_global-cut_inner;  
  double fc = e_zbl(p, cut_global, d1a, d2a, d3a, d4a, zze);
  double fcp = dzbldr(p, cut_global, d1a, d2a, d3a, d4a, zze);
  double fcpp = d2zbldr2(p, cut_global, d1a, d2a, d3a, d4a, zze);
  
  double swa = (-3.0*fcp + tc*fcpp)/(tc*tc);
  double swb = ( 2.0*fcp - tc*fcpp)/(tc*tc*tc);
  double swc = -fc + (tc/2.0)*fcp - (tc*tc/12.0)*fcpp;

  double sw1 = swa;
  double sw2 = swb;
  double sw3 = swa/3.0;
  double sw4 = swb/4.0;
  double sw5 = swc;

  // std::cout << "cut_global = " << cut_global << std::endl;
  // std::cout << "cut_inner  = " << cut_inner << std::endl;
  // std::cout << "tc         = " << tc << std::endl;
  // std::cout << "fc         = " << fc << std::endl;
  // std::cout << "fcp        = " << fcp << std::endl;
  // std::cout << "fcpp       = " << fcpp << std::endl;
  // std::cout << "sw1        = " << sw1 << std::endl;
  // std::cout << "sw2        = " << sw2 << std::endl;
  // std::cout << "sw3        = " << sw3 << std::endl;
  // std::cout << "sw4        = " << sw4 << std::endl;
  // std::cout << "sw5        = " << sw5 << std::endl << std::endl << std::endl;  
  
  double r,t,rsq,cut_innersq,cut_globalsq;
  double eswitch;
  double fswitch;
  
  rsq = rij * rij;
  cut_innersq = cut_inner * cut_inner;
  cut_globalsq = cut_global * cut_global;  

  // std::cout << "rsq          = " << rsq << std::endl;
  // std::cout << "cut_innersq  = " << cut_innersq << std::endl;
  // std::cout << "cut_globalsq = " << cut_globalsq << std::endl << std::endl << std::endl;
  
  if (rsq < cut_globalsq) {
    r = sqrt(rsq);
    de = dzbldr(p,r,d1a,d2a,d3a,d4a,zze);

    if (rsq > cut_innersq) {
      t = r - cut_inner;
      fswitch = t*t * (sw1 + sw2*t);

      // std::cout << "t       = " << t << std::endl;
      // std::cout << "de      = " << de << std::endl;
      
      de += fswitch;

      // std::cout << "fswitch = " << fswitch << std::endl;
      // std::cout << "de      = " << de << std::endl;
      
    }

    //    de *= -1.0/r;
    //    if (rsq < cut_innersq) {
    //    std::cout << "rsq     = " << rsq << std::endl;      
    //    std::cout << "de      = " << de << std::endl << std::endl;
    //    std::cout << "rsq , de      = " << rsq << " " << de << std::endl;    
      //    }
    e = e_zbl(p, r, d1a, d2a, d3a, d4a, zze);

    // if (r < 3.2) {
    //   std::cout << "rij   = " << r << std::endl;
    //   std::cout << "e_zbl = " << e << std::endl << std::endl;
    // }    
    e += sw5;


    if (rsq > cut_innersq) {
      t = r - cut_inner;
      eswitch = t*t*t * (sw3 + sw4*t);
      e += eswitch;
    }

  }

  static const double conv_energy_inv =  1e-4 * exanb::legacy_constant::elementaryCharge / exanb::legacy_constant::atomicMass;

  e *= conv_energy_inv;
  de *= conv_energy_inv;

}

}

#define USTAMP_POTENTIAL_NAME     zbl
#define USTAMP_POTENTIAL_PARAMS   ZBLParms
#define USTAMP_POTENTIAL_COMPUTE  zbl_compute_energy

// only atom z are meaningful for zbl
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) onika::make_flat_tuple(p.m_atom_a.m_z,p.m_atom_b.m_z)

#define USTAMP_POTENTIAL_ENABLE_CUDA 1

// static inline double zbl_V(const ZBLParms& p, const exanb::PairPotentialParameters& pair, double rij )
// {
//   return p.ct1 * ( pair.m_atom_a.m_z * pair.m_atom_b.m_z * p.e * p.e ) / rij ;
// }

// static inline double zbl_dV(const ZBLParms& p, const exanb::PairPotentialParameters& pair, double rij )
// {
//   return - p.ct1 * ( pair.m_atom_a.m_z * pair.m_atom_b.m_z * p.e * p.e ) / ( rij * rij );
// }

// static inline double zbl_ddV(const ZBLParms& p, const exanb::PairPotentialParameters& pair, double rij )
// {
//   return 2. * p.ct1 * ( pair.m_atom_a.m_z * pair.m_atom_b.m_z * p.e * p.e * p.e ) / ( rij * rij * rij );
// }

// static inline double zbl_Phi(const ZBLParms& p, double rij, double a )
// {
//   double x = rij / a;
//   return p.phi1 * std::exp( -p.phi2 * x)
//        + p.phi3 * std::exp( -p.phi4 * x)
//        + p.phi5 * std::exp( -p.phi6 * x)
//        + p.phi7 * std::exp( -p.phi8 * x) ;
// }

// static inline double zbl_dPhi(const ZBLParms& p, double rij, double a )
// {
//   double x = rij / a;
//   return ( - p.phi1 * p.phi2 * std::exp( -p.phi2 * x)
//            - p.phi3 * p.phi4 * std::exp( -p.phi4 * x)
//            - p.phi5 * p.phi6 * std::exp( -p.phi6 * x)
//            - p.phi7 * p.phi8 * std::exp( -p.phi8 * x) )
//          / a;
// }

// static inline double zbl_ddPhi(const ZBLParms& p, double rij, double a )
// {
//   double x = rij / a;
//   return ( - p.phi1 * p.phi2 * p.phi2 * std::exp( -p.phi2 * x)
//            - p.phi3 * p.phi4 * p.phi4 * std::exp( -p.phi4 * x)
//            - p.phi5 * p.phi6 * p.phi6 * std::exp( -p.phi6 * x)
//            - p.phi7 * p.phi8 * p.phi8 * std::exp( -p.phi8 * x) )
//          / ( a * a );
// }

// static inline double zbl_E(const ZBLParms& p, const exanb::PairPotentialParameters& pair, double rij, double a )
// {
//   return zbl_V(p,pair,rij) * zbl_Phi(p,rij,a) ;
// }

// static inline double zbl_dE(const ZBLParms& p, const exanb::PairPotentialParameters& pair, double rij, double a )
// {
//   return zbl_dV(p,pair,rij) * zbl_Phi(p,rij,a) + zbl_V(p,pair,rij) * zbl_dPhi(p,rij,a);
// }

// static inline double zbl_ddE(const ZBLParms& p, const exanb::PairPotentialParameters& pair, double rij, double a )
// {
//   return        zbl_ddV(p,pair,rij) * zbl_Phi(p,rij,a)
//         + 2.0 * zbl_dV(p,pair,rij)  * zbl_dPhi(p,rij,a)
//         +       zbl_V(p,pair,rij)   * zbl_ddPhi(p,rij,a);
// }

// static inline void initialize_ABC(const ZBLParms& p, const exanb::PairPotentialParameters& pair, double a, double& A, double& B, double& C)
// {
//   double rc_r1 = p.rc - p.r1;
//   double rc_r1_2 = rc_r1 * rc_r1;
//   double rc_r1_3 = rc_r1_2 * rc_r1;
  
//   double E_rc = zbl_E(p,pair,p.rc,a);
//   double dE_rc = zbl_dE(p,pair,p.rc,a);
//   double ddE_rc = zbl_ddE(p,pair,p.rc,a);
  
//   A = ( -3.0 * dE_rc + rc_r1 * ddE_rc ) / rc_r1_2;
//   B = ( 2.0 * dE_rc - rc_r1 * ddE_rc ) / rc_r1_3;
//   C = -E_rc + 0.5*rc_r1*dE_rc - (1.0/12.0) * rc_r1_2 * ddE_rc;
// }

// static inline double zbl_dS(const ZBLParms& p, double rij, double A, double B, double C)
// {
//   if( rij < p.r1 )
//   {
//     return 0.0;
//   }
//   else
//   {
//     if( rij > p.rc ) { rij = p.rc; }
//     double rij_r1 = rij - p.r1;
//     double rij_r1_2 = rij_r1 * rij_r1;
//     double rij_r1_3 = rij_r1_2 * rij_r1;
//     return A*rij_r1_2 + B*rij_r1_3;
//   }
// }

// static inline double zbl_S(const ZBLParms& p, double rij, double A, double B, double C)
// {
//   if( rij < p.r1 )
//   {
//     return C;
//   }
//   else
//   {
//     if( rij > p.rc ) { rij = p.rc; }
//     double rij_r1 = rij - p.r1;
//     double rij_r1_2 = rij_r1 * rij_r1;
//     double rij_r1_3 = rij_r1_2 * rij_r1;
//     double rij_r1_4 = rij_r1_2 * rij_r1_2;
//     return (A/3.0)*rij_r1_3 + (B/4.0)*rij_r1_4 + C;
//   }
// }
