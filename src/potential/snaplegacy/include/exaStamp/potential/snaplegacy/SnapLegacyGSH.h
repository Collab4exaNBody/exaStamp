#pragma once

#include <vector>                   
#include <complex>                  
#include <exaStamp/potential/snaplegacy/snap_legacy_3Dtypes.h>
                                    
using namespace std;                

class SnapLegacyGSH {

  public:
//    SnapLegacyGSH(double const &jmax=0.) : m_jmax(jmax), ran(), rcut(0.), gsh() {};	//constructor
    SnapLegacyGSH(double const &jmax=0.) : m_jmax(jmax) {};	//constructor
    int compute_gsh(double3d const &ran, double rcut);	//computes all gsh and dgsh coefficients; rcut could be transformed in const double input of the class?
    complex<double> gsh_val(int j, int m1, int m2) {	// returns gsh value for (j,m1,m2)
       return m_gsh[idx(j,m1,m2)];
    }
    complex3d       dgsh_val(int j, int m1, int m2) {	// returns dgsh value for (j,m1,m2)
       return m_dgsh[idx(j,m1,m2)];
    }
    int idx(int j, int m1, int m2) {	//returns index in gsh vector corresponding to (j,m1,m2)
        int idx=j*(j+1)*(2*j+1)/6+(m2+j)/2*(j+1)+(m1+j)/2+2; //sum of j first square integers+offset
        return idx;
    }
    int size();         	           //returns gsh vector size
    void set_jmax(double J) {m_jmax=J;}
  private:
    double m_jmax;	// input parameter jmax
    vector<complex<double>> m_gsh;       //contains all gsh coefficients 
    vector<complex3d>       m_dgsh;      //contains all dgsh coefficients
};
