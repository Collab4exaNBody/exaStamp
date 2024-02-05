#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>

namespace exaStamp
{

  struct EamAlloyParameters
  {
    std::vector<std::vector<std::vector<double>>> frho_spline;
    std::vector<std::vector<std::vector<double>>> rhor_spline;
    std::vector<std::vector<std::vector<double>>> z2r_spline;
    double rdr,rdrho,rc,rhomax;
    int nmats,nr,nrho;
  };

  static inline void eam_alloy_phi(const EamAlloyParameters& eam, double r, double& phi, double& dphi)
  {
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = std::min(m,eam.nr-1);
    p -= m;
    p = std::min(p,1.0);
    // rhoip = derivative of (density at atom j due to atom i)
    // rhojp = derivative of (density at atom i due to atom j)
    // phi = pair potential energy
    // phip = phi'
    // z2 = phi * r
    // z2p = (phi * r)' = (phi' r) + phi
    // psip needs both fp[i] and fp[j] terms since r_ij appears in two
    //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
    //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

    //    std::vector<double> coeff = rhor_spline[type2rhor[itype][jtype]][m];
    std::vector<double> coeff = eam.rhor_spline[0][m];
    double rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
    //    std::vector<double> coeff = rhor_spline[type2rhor[jtype][itype]][m];
    coeff.clear();
    coeff = eam.rhor_spline[0][m];    
    double rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
    //    std::vector<double> coeff = z2r_spline[type2z2r[itype][jtype]][m];
    coeff.clear();
    coeff = eam.z2r_spline[0][m];    
    double z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
    double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

    double recip = 1.0/r;
    phi = z2*recip;
    double phip = z2p*recip - phi*recip;
    //fp[i] fp [j] are computed elswhere. Need to find workaround ?
    double psip = /*fp[i]**/rhojp + /*fp[j]**/rhoip + phip;
    dphi = -psip*recip;
  }

  static inline void eam_alloy_rho(const EamAlloyParameters& eam, double r, double& rho, double& drho)
  {
    // Would need the types of central and neighbor atoms in this function (see below)
    
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = std::min(m,eam.nr-1);
    p -= m;
    p = std::min(p,1.0);
    // In Lammps type2rhor[jtype][itype] allows to define which rhor_spline should be used depending on atom types
    //    coeff = eam.rhor_spline[type2rhor[jtype][itype]]][m];
    std::vector<double> coeff = eam.rhor_spline[0][m];
    rho = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    drho=0.;
  }

  static inline void eam_alloy_fEmbed(const EamAlloyParameters& eam, double rho, double& f, double& df)
  {
    // Would need the type of central atom in this function (see below)
    
    double p = rho * eam.rdrho + 1.0;
    int m = static_cast<int> (p);
    m = std::max(1,std::min(m,eam.nrho-1));
    p -= m;
    p = std::min(p,1.0);
    // Again type2frho[type[i]] allows to choose the appropriate frho_spline function depending on central atom type
    //    coeff = eam.frho_spline[type2frho[type[i]]][m];
    std::vector<double> coeff = eam.frho_spline[0][m];
    // fp =derivative of embedding energy at each atom
    // phi = embedding energy at each atom
    df = (coeff[0]*p + coeff[1])*p + coeff[2];
    f = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    if (rho > eam.rhomax) f += df * (rho-eam.rhomax);
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EamAlloyParameters;
  using exanb::UnityConverterHelper;
  using exanb::Quantity;

  template<> struct convert<EamAlloyParameters>
    {
      static bool decode(const Node& node, EamAlloyParameters& v);
      static void interpolate(int n, double delta, std::vector<double> f, std::vector<std::vector<double>>& spline);
    };
}

