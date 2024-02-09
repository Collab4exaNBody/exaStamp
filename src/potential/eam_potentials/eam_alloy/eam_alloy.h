#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>
#include <exanb/core/physics_constants.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>
#include <exaStamp/potential/eam/eam_buffer.h>
#include <iostream>

namespace exaStamp
{

  struct EamAlloyParameters
  {
    static inline constexpr size_t MAX_ARRAY_SIZE = 8;
    static inline constexpr size_t MAX_ELEMENTS = MAX_ARRAY_SIZE - 1;
    double conversion_z2r = UnityConverterHelper::convert(1., "eV*ang");
    double conversion_frho = UnityConverterHelper::convert(1., "eV");
    
    std::vector< std::vector< std::vector<double> > > frho_spline;
    std::vector< std::vector< std::vector<double> > > rhor_spline;
    std::vector< std::vector< std::vector<double> > > z2r_spline;
      
    double rdr = 0.0;
    double rdrho = 0.0;
    double rc = 0.0;
    double rhomax = 0.0;
    int nelements = 0;
    int nr = 0;
    int nrho = 0;

    int type2frho[MAX_ARRAY_SIZE]; // Max 7 materials
    int type2rhor[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
    int type2z2r[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
    int n_types = 0;
    
    EamAlloyParameters() = default;
    EamAlloyParameters(const EamAlloyParameters&) = default;

    inline EamAlloyParameters(const EamAlloyParameters& p, int ntypes , const uint8_t* pair_enabled)
    {
      if( static_cast<size_t>(ntypes) > MAX_ELEMENTS )
      {
        std::cerr << "Maximum number of elements is "<<MAX_ELEMENTS<<" , but ntypes="<<ntypes<<std::endl;
        std::abort();
      }
    
      *this = p;
    
      n_types = ntypes;
      for(size_t i=0;i<MAX_ARRAY_SIZE;i++)
      {
        type2frho[i] = -1;
        for(size_t j=0;j<MAX_ARRAY_SIZE;j++)
        {
          type2rhor[i][j] = -1;
          type2z2r[i][j] = -1;
        }
      }
    
      std::vector<int> map( ntypes + 1 , -1 );
      for(int i=1 ; i<=ntypes ; i++) {
        map[i] = i-1;
      }
      
      int nfrho = frho_spline.size();
    
      // type2frho
      for (int i = 1; i <= ntypes; i++) {
        if (map[i] >= 0)
          type2frho[i] = map[i];
        else
          type2frho[i] = nfrho - 1;
      }

      // type2rhor
      for (int i = 1; i <= ntypes; i++) {
        for (int j = 1; j <= ntypes; j++) {
          type2rhor[i][j] = map[i];
        }
      }
    
      // type2z2r
      int irow, icol;
      for (int i = 1; i <= ntypes; i++) {
        for (int j = 1; j <= ntypes; j++) {
          irow = map[i];
          icol = map[j];
          if (irow == -1 || icol == -1) {
            type2z2r[i][j] = 0;
            continue;
          }
          if (irow < icol) {
            irow = map[j];
            icol = map[i];
          }
          int n = 0;
          for (int m = 0; m < irow; m++) n += m + 1;
          n += icol;
          type2z2r[i][j] = n;
        }
      }
    }
  };

  namespace EamAlloyTools
  {  
    // ONIKA_HOST_DEVICE_FUNC
    static inline void interpolate(int n, double delta, const std::vector<double>& f, std::vector< std::vector<double> > & spline)
    {

      for (int m = 1; m <= n; m++) spline[m][6] = f[m];

      spline[1][5] = spline[2][6] - spline[1][6];
      spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
      spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
      spline[n][5] = spline[n][6] - spline[n-1][6];

      for (int m = 3; m <= n-2; m++)
        spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                        8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

      for (int m = 1; m <= n-1; m++) {
        spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
          2.0*spline[m][5] - spline[m+1][5];
        spline[m][3] = spline[m][5] + spline[m+1][5] -
          2.0*(spline[m+1][6]-spline[m][6]);
      }

      spline[n][4] = 0.0;
      spline[n][3] = 0.0;

      for (int m = 1; m <= n; m++) {
        spline[m][2] = spline[m][5]/delta;
        spline[m][1] = 2.0*spline[m][4]/delta;
        spline[m][0] = 3.0*spline[m][3]/delta;
      }
      
    }

  }

  // ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_phi_mm(const EamAlloyParameters& eam, double r, double& phi, double& dphi, int itype , int jtype )
  {
    // ONIKA_HOST_DEVICE_FUNC
    using onika::cuda::min;
  
    itype = itype + 1;
    jtype = jtype + 1;

    assert( static_cast<size_t>(eam.n_types) < EamAlloyParameters::MAX_ARRAY_SIZE );
    assert( itype >= 1 && itype <= eam.n_types );
    assert( jtype >= 1 && jtype <= eam.n_types );
  
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = min(m,eam.nr-1);
    p -= m;
    p = min(p,1.0);
        
    const int z2r_i_j_index = eam.type2z2r[itype][jtype];
    assert( z2r_i_j_index >= 0 && static_cast<size_t>(z2r_i_j_index) < eam.z2r_spline.size() );
    const auto& coeff = eam.z2r_spline[z2r_i_j_index][m];    
    const double z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
    const double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    const double recip = 1.0/r;
    
    phi = z2*recip;
    dphi = z2p*recip - phi*recip;

    phi *= eam.conversion_z2r;
    dphi *= eam.conversion_z2r;    
  }

  static inline void eam_alloy_phi(const EamAlloyParameters& eam, double r, double& phi, double& dphi)
  {
    eam_alloy_phi_mm(eam,r,phi,dphi,0,0);
  }

  // ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_rho_mm(const EamAlloyParameters& eam, double r, double& rho, double& drho, int itype , int jtype)
  {
    using onika::cuda::min;
    itype = itype + 1;
    jtype = jtype + 1;
    assert( static_cast<size_t>(eam.n_types) < EamAlloyParameters::MAX_ARRAY_SIZE );
    assert( itype >= 1 && itype <= eam.n_types );
    assert( jtype >= 1 && jtype <= eam.n_types );
    
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = min(m,eam.nr-1);
    p -= m;
    p = min(p,1.0);

    const int rhor_i_j_index = eam.type2rhor[itype][jtype];
    assert( rhor_i_j_index >= 0 && static_cast<size_t>(rhor_i_j_index) < eam.rhor_spline.size() );
    const auto& coeff = eam.rhor_spline[rhor_i_j_index][m];
    rho = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    drho = (coeff[0]*p + coeff[1])*p + coeff[2];    
  }

  static inline void eam_alloy_rho(const EamAlloyParameters& eam, double r, double& rho, double& drho)
  {
    eam_alloy_rho_mm(eam,r,rho,drho,0,0);
  }

  // ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_fEmbed_mm(const EamAlloyParameters& eam, double rho, double& phi, double& fp, int itype )
  {
    using onika::cuda::min;
    using onika::cuda::max;

    itype = itype + 1;
    assert( static_cast<size_t>(eam.n_types) < EamAlloyParameters::MAX_ARRAY_SIZE ); 
    assert( itype >= 1 && itype <= eam.n_types );
    
    double p = rho * eam.rdrho + 1.0;
    int m = static_cast<int> (p);
    m = max(1,min(m,eam.nrho-1));
    p -= m;
    p = min(p,1.0);
    
    const int frho_i_j_index = eam.type2frho[itype];
    assert( frho_i_j_index >= 0 && static_cast<size_t>(frho_i_j_index) < eam.frho_spline.size() );
    const auto& coeff = eam.frho_spline[frho_i_j_index][m];
    fp = (coeff[0]*p + coeff[1])*p + coeff[2];
    phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    if (rho > eam.rhomax) phi += fp * (rho-eam.rhomax);
    phi *= eam.conversion_frho;
    fp *= eam.conversion_frho;
  }

  static inline void eam_alloy_fEmbed(const EamAlloyParameters& eam, double rho, double& f, double& df)
  {
    eam_alloy_fEmbed_mm(eam,rho,f,df,0);
  }  
}

#define USTAMP_POTENTIAL_EAM_EMB_MM_CONSTRUCTOR \
inline EmbOp(const EamAlloyParameters& parms, size_t nt, const uint8_t * __restrict__ pair_enabled ) \
: p(parms,nt,pair_enabled) , n_types(nt) , m_pair_enabled(pair_enabled) { }

#define USTAMP_POTENTIAL_EAM_FORCE_MM_CONSTRUCTOR \
inline ForceOp(const EamAlloyParameters& parms, size_t nt, const uint8_t * __restrict__ pair_enabled ) \
: p(parms,nt,pair_enabled) , n_types(nt) , m_pair_enabled(pair_enabled) { }

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EamAlloyParameters;
  template<> struct convert<EamAlloyParameters>
    {
      static bool decode(const Node& node, EamAlloyParameters& v);
    };
}

