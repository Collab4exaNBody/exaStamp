/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <vector>
#include <string>

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

