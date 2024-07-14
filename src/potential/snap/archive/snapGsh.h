#pragma once

#include <vector>                   
#include <cmath>                  
#include "snap_3Dtypes.h"
         
namespace SnapExt
{

                           
class snapGsh {

  public:
//    snapGsh(double const &jmax=0.) : m_jmax(jmax), ran(), rcut(0.), gsh() {};	//constructor
    snapGsh(double jmax=0.)
      : m_jmax(jmax)
      , m_two_jmax( std::floor(2*jmax) )
      {
        m_gsh.resize(size());
        m_dgsh.resize(size());
      };	//constructor
      
    int compute_gsh(double3d const &ran, double rcut, double rfac0, double rmin0);	//computes all gsh and dgsh coefficients; rcut could be transformed in const double input of the class?

    inline const Complexd * __restrict__ gsh_data() const { return m_gsh.data(); }
    inline const complex3d * __restrict__ dgsh_data() const { return m_dgsh.data(); }

    inline Complexd gsh_val(int j, int m1, int m2)
    {	// returns gsh value for (j,m1,m2)
       return m_gsh[idx(j,m1,m2)];
    }
    
    inline complex3d dgsh_val(int j, int m1, int m2)
    {	// returns dgsh value for (j,m1,m2)
       return m_dgsh[idx(j,m1,m2)];
    }
    
    static inline constexpr int idx(int j, int m1, int m2) 
    {	//returns index in gsh vector corresponding to (j,m1,m2)
        return j*(j+1)*(2*j+1)/6+(m2+j)/2*(j+1)+(m1+j)/2+2; //sum of j first square integers+offset
    }
    
    inline int size() const
    {
          size_t nb_j = m_two_jmax + 1; // Account for j = 0.
          size_t nb_m_per_j = m_two_jmax + 1; // -j <= m <= j with cnt 2.
          return nb_j*nb_m_per_j*nb_m_per_j; // Upper bound.
    }
    
    void set_jmax(double J)
    {
      m_jmax=J;
      m_two_jmax=std::floor(2*J);
      m_gsh.resize(size());
      m_dgsh.resize(size());
    }
    
  private:
    double m_jmax;	// input parameter jmax
    int m_two_jmax;
    std::vector< Complexd > m_gsh; //contains all gsh coefficients 
    std::vector< complex3d > m_dgsh; //contains all dgsh coefficients
};


}

