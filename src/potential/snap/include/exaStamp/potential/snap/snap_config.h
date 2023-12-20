#pragma once

#include <exaStamp/potential/snap/snap_material.h>
#include <vector>

namespace SnapExt
{

class SnapConfig
{
  public:
    inline void set_rfac0(double x) { m_rfac0=x; }
    inline double rfac0() const { return m_rfac0; }

    inline void set_rmin0(double x) { m_rmin0=x; }
    inline double rmin0() const { return m_rmin0; }

    inline void set_bzeroflag(bool x) { m_bzeroflag = x; }
    inline bool bzeroflag() const { return m_bzeroflag; }
  
    inline void set_rcutfac(double x) { m_rcutfac=x; }
    inline double rcutfac() const { return m_rcutfac; }    
 
    inline void set_twojmax(size_t x) { m_twojmax=x; }
    inline size_t twojmax() const { return m_twojmax; }     
    
    inline std::vector<SnapMaterial>& materials() { return m_materials; }
    inline const std::vector<SnapMaterial>& materials() const { return m_materials; }

    inline void set_nelements(int n) { m_nelements = n; }
    inline int nelements() const { return m_nelements; }
    
    inline void set_switchflag(bool b) { m_switchflag = b; }
    inline bool switchflag() const { return m_switchflag; }

    inline void set_chemflag(bool b) { m_chemflag = b; }
    inline bool chemflag() const { return m_chemflag; }

    inline void set_bnormflag(bool b) { m_bnormflag = b; }
    inline bool bnormflag() const { return m_bnormflag; }

    inline void set_wselfallflag(bool b) { m_wselfallflag = b; }
    inline bool wselfallflag() const { return m_wselfallflag; }

    inline void set_switchinnerflag(bool b) { m_switchinnerflag = b; }
    inline bool switchinnerflag() const { return m_switchinnerflag; }

    inline void set_quadraticflag(bool b) { m_quadraticflag = b; }
    inline bool quadraticflag() const { return m_quadraticflag; }

  private:
    std::vector<SnapMaterial> m_materials;
    double m_rfac0 = 0.99363;
    double m_rmin0 = 0.;
    double m_rcutfac = 0.;
    size_t m_twojmax = 6;
    int m_nelements = 1;
    bool m_bzeroflag = true;
    
    // extended parameters
    bool m_switchflag = true;
    bool m_chemflag = false;
    bool m_bnormflag = false;
    bool m_wselfallflag = false;
    bool m_switchinnerflag = false;
    bool m_quadraticflag = false;
};

}
