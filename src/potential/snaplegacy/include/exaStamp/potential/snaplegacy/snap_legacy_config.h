#pragma once

#include <exaStamp/potential/snaplegacy/snap_legacy_material.h>
#include <vector>

class SnapConfig
{
  public:
    inline void set_rfac0(double x) { m_rfac0=x; }
    inline double rfac0() const { return m_rfac0; }

    inline void set_rmin0(double x) { m_rmin0=x; }
    inline double rmin0() const { return m_rmin0; }
  
    inline void set_rcutfac(double x) { m_rcutfac=x; }
    inline double rcutfac() const { return m_rcutfac; }    
 
    inline void set_twojmax(size_t x) { m_twojmax=x; }
    inline size_t twojmax() const { return m_twojmax; }     
    
    inline std::vector<SnapMaterial>& materials() { return m_materials; }
    inline const std::vector<SnapMaterial>& materials() const { return m_materials; }

  private:
    double m_rfac0 = 1.;
    double m_rmin0 = 0.;
    double m_rcutfac = 0.;
    size_t m_twojmax = 6;
    std::vector<SnapMaterial> m_materials;
};

