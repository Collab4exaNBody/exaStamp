#pragma once

#include <vector>
#include <string>

namespace SnapExt
{

class SnapMaterial
{
  public:
  inline void set_radelem(double r) { m_radelem=r; }
  inline double radelem() const { return m_radelem; }

  inline void set_weight(double w) { m_weight=w; }
  inline double weight() const { return m_weight; }

  inline void set_name(const std::string& s) { m_name = s; }
  inline const std::string& name() const { return m_name; }

  inline size_t number_of_coefficients() const { return m_coeffs.size(); }
  inline void resize_coefficients(size_t s) { m_coeffs.assign( s , 0. ); }
  inline void set_coefficient(size_t i, double value) { m_coeffs[i] = value; }
  inline double coefficient(size_t i) const { return m_coeffs[i]; }
  inline double const * coefficient_data() const { return m_coeffs.data(); }

  private:
	double m_radelem = 0.5;
	double m_weight = 1.0;
	std::vector<double> m_coeffs;
	std::string m_name = "<unknown>";
};


}
