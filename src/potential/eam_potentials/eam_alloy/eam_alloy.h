#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>
#include <exaStamp/potential/eam/eam_buffer.h>

namespace exaStamp
{

  struct EamAlloyParameters
  {
    std::vector< std::vector< std::vector<double> > > frho_spline;
    std::vector< std::vector< std::vector<double> > > rhor_spline;
    std::vector< std::vector< std::vector<double> > > z2r_spline;
    double rdr = 0.0;
    double rdrho = 0.0;
    double rc = 0.0;
    double rhomax = 0.0;
    int nmats = 0;
    int nr = 0;
    int nrho = 0;
  };

  namespace EamAlloyTools
  {
  
    // ONIKA_HOST_DEVICE_FUNC
    static inline void interpolate(int n, double delta, const std::vector<double>& f, std::vector< std::vector<double> > & spline)
    {
      for (int m = 0; m < n; m++) spline[m][6] = f[m];

      spline[0][5] = spline[1][6] - spline[0][6];
      spline[0][5] = 0.5 * (spline[2][6]-spline[0][6]);
      spline[n-2][5] = 0.5 * (spline[n-1][6]-spline[n-3][6]);
      spline[n-1][5] = spline[n-1][6] - spline[n-2][6];

      for (int m = 2; m < n-2; m++)
        spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) + 8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

      for (int m = 0; m < n-1; m++) {
        spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) - 2.0*spline[m][5] - spline[m+1][5];
        spline[m][3] = spline[m][5] + spline[m+1][5] - 2.0*(spline[m+1][6]-spline[m][6]);
      }

      spline[n-1][4] = 0.0;
      spline[n-1][3] = 0.0;

      for (int m = 0; m < n; m++) {
        spline[m][2] = spline[m][5]/delta;
        spline[m][1] = 2.0*spline[m][4]/delta;
        spline[m][0] = 3.0*spline[m][3]/delta;
      }
    }
  }

  // ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_phi_mm(const EamAlloyParameters& eam, double r, double& phi, double& dphi, const EAMSpecyPairInfo& pair_info , bool pair_inversed = false)
  {
    // ONIKA_HOST_DEVICE_FUNC
    using onika::cuda::min;
  
    int itype=0, jtype=0;
    if( pair_inversed ) { itype = pair_info.m_type_b; jtype = pair_info.m_type_a; }
    else                { itype = pair_info.m_type_a; jtype = pair_info.m_type_b; }
  
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = min(m,eam.nr-1);
    p -= m;
    p = min(p,1.0);
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
    const int rhor_i_j_index = 0; // = type2rhor[itype][jtype]
    const auto& coeff_rhor_i_j = eam.rhor_spline[rhor_i_j_index][m];
    double rhoip = (coeff_rhor_i_j[0]*p + coeff_rhor_i_j[1])*p + coeff_rhor_i_j[2];

    //    std::vector<double> coeff = rhor_spline[type2rhor[jtype][itype]][m];
    const int rhor_j_i_index = 0; // = type2rhor[jtype][itype]
    const auto& coeff_rhor_j_i = eam.rhor_spline[rhor_j_i_index][m];
    const double rhojp = (coeff_rhor_j_i[0]*p + coeff_rhor_j_i[1])*p + coeff_rhor_j_i[2];

    //    std::vector<double> coeff = z2r_spline[type2z2r[itype][jtype]][m];
    const int z2r_i_j_index = 0; // = type2z2r[itype][jtype]
    const auto& coeff_z2r_i_j = eam.z2r_spline[z2r_i_j_index][m];    
    const double z2p = (coeff_z2r_i_j[0]*p + coeff_z2r_i_j[1])*p + coeff_z2r_i_j[2];
    const double z2 = ((coeff_z2r_i_j[3]*p + coeff_z2r_i_j[4])*p + coeff_z2r_i_j[5])*p + coeff_z2r_i_j[6];

    const double recip = 1.0/r;
    phi = z2*recip;
    const double phip = z2p*recip - phi*recip;
    
    //fp[i] fp [j] are computed elswhere. Need to find workaround ?
    const double psip = /*fp[i]**/rhojp + /*fp[j]**/rhoip + phip;
    dphi = -psip*recip;
  }

  static inline void eam_alloy_phi(const EamAlloyParameters& eam, double r, double& phi, double& dphi)
  {
    eam_alloy_phi_mm(eam,r,phi,dphi,EAMSpecyPairInfo{},false);
  }

  // ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_rho_mm(const EamAlloyParameters& eam, double r, double& rho, double& drho, const EAMSpecyPairInfo& pair_info , bool pair_inversed = false)
  {
    using onika::cuda::min;
    // Would need the types of central and neighbor atoms in this function (see below)

    int itype=0, jtype=0;
    if( pair_inversed ) { itype = pair_info.m_type_b; jtype = pair_info.m_type_a; }
    else                { itype = pair_info.m_type_a; jtype = pair_info.m_type_b; }
    
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = min(m,eam.nr-1);
    p -= m;
    p = min(p,1.0);
    
    // In Lammps type2rhor[jtype][itype] allows to define which rhor_spline should be used depending on atom types
    //    coeff = eam.rhor_spline[type2rhor[jtype][itype]]][m];
    const int rhor_i_j_index = 0; // = type2rhor[jtype][itype];
    const auto& coeff = eam.rhor_spline[rhor_i_j_index][m];
    rho = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    drho=0.;
  }

  static inline void eam_alloy_rho(const EamAlloyParameters& eam, double r, double& rho, double& drho)
  {
    eam_alloy_rho_mm(eam,r,rho,drho,EAMSpecyPairInfo{},false);
  }

  // ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_fEmbed(const EamAlloyParameters& eam, double rho, double& f, double& df)
  {
    using onika::cuda::min;
    using onika::cuda::max;
    // Would need the type of central atom in this function (see below)
    
    double p = rho * eam.rdrho + 1.0;
    int m = static_cast<int> (p);
    m = max(1,min(m,eam.nrho-1));
    p -= m;
    p = min(p,1.0);
    
    // Again type2frho[type[i]] allows to choose the appropriate frho_spline function depending on central atom type
    //    coeff = eam.frho_spline[type2frho[type[i]]][m];
    const int frho_i_j_index = 0; // type2frho[type[i]]
    const auto& coeff = eam.frho_spline[frho_i_j_index][m];
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
  template<> struct convert<EamAlloyParameters>
    {
      static bool decode(const Node& node, EamAlloyParameters& v);
    };
}

