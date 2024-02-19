#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>
#include <exanb/core/physics_constants.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <onika/memory/allocator.h>
#include <exaStamp/potential/eam/eam_buffer.h>
#include <iostream>

namespace exaStamp
{

  struct alignas(sizeof(double)*8) SplineCoeffs
  {
    static constexpr size_t N_SPLINE_POINTS_STORAGE = 8;
    alignas(sizeof(double)*8) double coeffs[N_SPLINE_POINTS_STORAGE] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  };
  static_assert( sizeof(SplineCoeffs) == ( SplineCoeffs::N_SPLINE_POINTS_STORAGE * sizeof(double) ) );

  struct EamAlloyParameters
  {
    static constexpr size_t N_SPLINE_POINTS = 7;
    static_assert( N_SPLINE_POINTS <= SplineCoeffs::N_SPLINE_POINTS_STORAGE );
         
    std::vector< std::vector< std::vector<double> > > frho_spline;
    std::vector< std::vector< std::vector<double> > > rhor_spline;
    std::vector< std::vector< std::vector<double> > > z2r_spline;

    onika::memory::CudaMMVector< SplineCoeffs > frho_spline_data;
    onika::memory::CudaMMVector< SplineCoeffs > rhor_spline_data;
    onika::memory::CudaMMVector< SplineCoeffs > z2r_spline_data;

    double rdr = 0.0;
    double rdrho = 0.0;
    double rc = 0.0;
    double rhomax = 0.0;
    int nelements = 0;
    int nr = 0;
    int nrho = 0;
    int n_types = 0;
    
    std::vector< std::vector<int> > type2rhor;
    std::vector< std::vector<int> > type2z2r;
    std::vector<int> type2frho;

    inline void initialize_types_table(int nt, const uint8_t* pair_enabled)
    {
      assert( nt >= nelements );
      n_types = nt;

      type2frho.resize(n_types+1);
      type2rhor.resize(n_types+1);
      type2z2r.resize(n_types+1);
      for(int i=0 ; i <= n_types ; i++)
      {
        type2frho[i] = -1;
        type2rhor[i].resize(n_types+1);
        type2z2r[i].resize(n_types+1);
        for(int j=0 ; j <= n_types ; j++)
        {
          type2rhor[i][j] = -1;
          type2z2r[i][j] = -1;
        }
      }
    
      // TODO: change this to reflect pair_enabled
      std::vector<int> map( n_types + 1 , -1 );
      for(int i=1 ; i<= n_types ; i++) {
        map[i] = i-1;
      }
      
      const int nfrho = frho_spline.size();
    
      // type2frho
      for (int i = 1; i <= n_types; i++) {
        if (map[i] >= 0)
          type2frho[i] = map[i];
        else
          type2frho[i] = nfrho - 1;
      }

      // type2rhor
      for (int i = 1; i <= n_types; i++) {
        for (int j = 1; j <= n_types; j++) {
          type2rhor[i][j] = map[i];
        }
      }
    
      // type2z2r
      int irow, icol;
      for (int i = 1; i <= n_types; i++) {
        for (int j = 1; j <= n_types; j++) {
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

  struct EamAlloyParametersRO
  {
    static inline constexpr size_t MAX_ELEMENTS = 7;
    static inline constexpr size_t MAX_ARRAY_SIZE = MAX_ELEMENTS+1;
    static constexpr size_t N_SPLINE_POINTS = EamAlloyParameters::N_SPLINE_POINTS;
    static_assert( MAX_ARRAY_SIZE <= 127 );
        
    const SplineCoeffs * const __restrict__ frho_spline_data = nullptr;
    const SplineCoeffs * const __restrict__ rhor_spline_data = nullptr;
    const SplineCoeffs * const __restrict__ z2r_spline_data = nullptr;
      
    const double rdr = 0.0;
    const double rdrho = 0.0;
    const double rc = 0.0;
    const double rhomax = 0.0;
    const int nelements = 0;
    const int nr = 0;
    const int nrho = 0;
    const int n_types = 0;

    const double conversion_z2r = 0.0;
    const double conversion_frho = 0.0;

    int8_t type2rhor[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
    int8_t type2z2r[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
    int8_t type2frho[MAX_ARRAY_SIZE];

    EamAlloyParametersRO() = default;
    EamAlloyParametersRO(const EamAlloyParametersRO&) = default;

    inline EamAlloyParametersRO(const EamAlloyParameters& p)
    : frho_spline_data( p.frho_spline_data.data() )
    , rhor_spline_data( p.rhor_spline_data.data() )
    , z2r_spline_data( p.z2r_spline_data.data() )
    , rdr( p.rdr )
    , rdrho( p.rdrho )
    , rc( p.rc )
    , rhomax( p.rhomax )
    , nelements( p.nelements )
    , nr( p.nr )
    , nrho( p.nrho )
    , n_types( p.n_types )
    , conversion_z2r( EXANB_QUANTITY( 1. * eV * ang ) )
    , conversion_frho( EXANB_QUANTITY( 1.0 * eV ) )
    {
    
      if( static_cast<size_t>(n_types) > MAX_ELEMENTS )
      {
        std::cerr << "Maximum number of elements is "<<MAX_ELEMENTS<<" , but n_types="<<n_types<<std::endl;
        std::abort();
      }

      for(int i=0;i<=n_types;i++)
      {
        type2frho[i] = p.type2frho[i];
        for(int j=0;j<=n_types;j++)
        {
          type2rhor[i][j] = p.type2rhor[i][j];
          type2z2r[i][j] = p.type2z2r[i][j];
        }
      }
    }
  };

  ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_phi(const EamAlloyParametersRO& eam, double r, double& phi, double& dphi, int itype=0 , int jtype=0 )
  {  
    using onika::cuda::min;
  
    itype = itype + 1;
    jtype = jtype + 1;

    assert( static_cast<size_t>(eam.n_types) < EamAlloyParametersRO::MAX_ARRAY_SIZE );
    assert( itype >= 1 && itype <= eam.n_types );
    assert( jtype >= 1 && jtype <= eam.n_types );
  
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = min(m,eam.nr-1);
    p -= m;
    p = min(p,1.0);
        
    const int z2r_i_j_index = eam.type2z2r[itype][jtype];
    //const auto& coeff = eam.z2r_spline[z2r_i_j_index][m];
    const auto * __restrict__ coeff = eam.z2r_spline_data[ z2r_i_j_index * (eam.nr+1) + m ].coeffs;
    ONIKA_ASSUME_ALIGNED( coeff );
    const double z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
    const double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    const double recip = 1.0/r;
    
    phi = z2*recip;
    dphi = z2p*recip - phi*recip;

    phi *= eam.conversion_z2r;
    dphi *= eam.conversion_z2r;    
  }

  ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_rho(const EamAlloyParametersRO& eam, double r, double& rho, double& drho, int itype=0 , int jtype=0 )
  {
    using onika::cuda::min;
    itype = itype + 1;
    jtype = jtype + 1;
    
    assert( static_cast<size_t>(eam.n_types) < EamAlloyParametersRO::MAX_ARRAY_SIZE );
    assert( itype >= 1 && itype <= eam.n_types );
    assert( jtype >= 1 && jtype <= eam.n_types );
    
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = min(m,eam.nr-1);
    p -= m;
    p = min(p,1.0);

    const int rhor_i_j_index = eam.type2rhor[itype][jtype];
    //const auto& coeff = eam.rhor_spline[rhor_i_j_index][m];
    const auto * __restrict__ coeff = eam.rhor_spline_data[ rhor_i_j_index * (eam.nr+1) + m ].coeffs;
    ONIKA_ASSUME_ALIGNED( coeff );
    rho = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    drho = (coeff[0]*p + coeff[1])*p + coeff[2];    
  }

  ONIKA_HOST_DEVICE_FUNC
  ONIKA_ALWAYS_INLINE
  static Vec3d eam_alloy_mm_force(const EamAlloyParametersRO& eam, Vec3d dr, double r, double fpi, double fpj, int itype=0 , int jtype=0 )
  {
    // to ensure here-after use of builtin_assume_aligned is legit
    static_assert( SplineCoeffs::N_SPLINE_POINTS_STORAGE*sizeof(double) == 64 && alignof(SplineCoeffs) == 64 );
    
    using onika::cuda::min;
    itype = itype + 1;
    jtype = jtype + 1;
    
    assert( static_cast<size_t>(eam.n_types) < EamAlloyParametersRO::MAX_ARRAY_SIZE );
    assert( itype >= 1 && itype <= eam.n_types );
    assert( jtype >= 1 && jtype <= eam.n_types );
    
    double p = r * eam.rdr + 1.0;
    int m = static_cast<int> (p);
    m = min(m,eam.nr-1);
    p -= m;
    p = min(p,1.0);

    const int rhor_i_j_index = eam.type2rhor[itype][jtype];
    const double* __restrict__ coeff_i_j = (const double*) __builtin_assume_aligned( eam.rhor_spline_data[ rhor_i_j_index * (eam.nr+1) + m ].coeffs , 64 );
    const double rhoip = (coeff_i_j[0]*p + coeff_i_j[1])*p + coeff_i_j[2];    

    const int rhor_j_i_index = eam.type2rhor[jtype][itype];
    const double* __restrict__ coeff_j_i = (const double*) __builtin_assume_aligned( eam.rhor_spline_data[ rhor_j_i_index * (eam.nr+1) + m ].coeffs , 64 );
    const double rhojp = (coeff_j_i[0]*p + coeff_j_i[1])*p + coeff_j_i[2];  

    const int z2r_i_j_index = eam.type2z2r[itype][jtype];
    const double* __restrict__ coeff = (const double*) __builtin_assume_aligned( eam.z2r_spline_data[ z2r_i_j_index * (eam.nr+1) + m ].coeffs , 64 );
    const double z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
    const double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    const double recip = 1.0/r;
    
    const double phi = z2*recip;
    const double phip = ( z2p*recip - phi*recip ) * eam.conversion_z2r;

    const double psip = fpi * rhojp + fpj * rhoip + phip;
    const double fpair = psip*recip;
    
    return { dr.x * fpair , dr.y * fpair , dr.z * fpair };
  }


  ONIKA_HOST_DEVICE_FUNC
  static inline void eam_alloy_fEmbed(const EamAlloyParametersRO& eam, double rho, double& phi, double& fp, int itype=0 )
  {
    using onika::cuda::min;
    using onika::cuda::max;

    itype = itype + 1;
    assert( static_cast<size_t>(eam.n_types) < EamAlloyParametersRO::MAX_ARRAY_SIZE ); 
    assert( itype >= 1 && itype <= eam.n_types );
    
    double p = rho * eam.rdrho + 1.0;
    int m = static_cast<int> (p);
    m = max(1,min(m,eam.nrho-1));
    p -= m;
    p = min(p,1.0);
    
    const int frho_i_j_index = eam.type2frho[itype];
    const auto * __restrict__ coeff = eam.frho_spline_data[ frho_i_j_index * (eam.nrho+1) + m ].coeffs; //eam.frho_spline[frho_i_j_index][m];
    ONIKA_ASSUME_ALIGNED( coeff );
    fp = (coeff[0]*p + coeff[1])*p + coeff[2];
    phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    if (rho > eam.rhomax) phi += fp * (rho-eam.rhomax);
    phi *= eam.conversion_frho;
    fp *= eam.conversion_frho;
  }

}

namespace onika
{
  namespace cuda
  {
    template<> struct ReadOnlyShallowCopyType<exaStamp::EamAlloyParameters> { using type = exaStamp::EamAlloyParametersRO; };
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

