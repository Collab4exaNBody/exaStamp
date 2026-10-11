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
#include <complex>                  
#include <cmath>

class SnapLegacyCG {

  public:
    SnapLegacyCG(double const &_jmax, int const &m_nt);		//constructor

    //returns cg(j,j1,j2,m,m1,m2)
    inline double val(int j, int j1, int j2, int m, int m1, int m2) const
    {  
      return m_cg_tab[index(j,j1,j2,m,m1,m2)];
    }

    inline const double* __restrict__ cg_tab_data() const { return m_cg_tab.data(); }

    //returns the number of cg coefficients (useless ?)
    inline int size()
    {
      return m_cg_tab.size();
    }

    //returns jmax
    inline double get_jmax() const
    {
       return m_jmax;
    }

    int compute(); 					 	//computes all the coefs

    // returns the index in the array of cg coefs
    inline int index(int j, int j1, int j2, int m, int m1, int m2) const
    {
      int idx_tabm= (m2+j2)/m_nt + (m1+j1)/m_nt*((2*j2)/m_nt+1)+(m+j)/m_nt*((2*j1)/m_nt+1)*((2*j2)/m_nt+1) ; //index in the subarray (m,m1,m2)
      return m_tab_idxj[idx_tabj(j,j1,j2)]+idx_tabm; // sum of the beginning of the subarray(offset)+index in the subarray
    }


    // returns the index in the array of offsets 
    inline int idx_tabj(int j, int j1, int j2) const
    {
      return j2+j1*(m_nt_jmax+1)+j*std::pow(m_nt_jmax+1,2);
    }

  private:

    // return f*n!
    int factorial(int const n, double & f);
    
    // computes wigner3j symbols
    int wigner3j(int const j, int const m, int const j1,int const m1, int const j2,int const m2, double & w3j);
		 
		// computes one coefficient
    int compute_cg(int const j,int const m, int const j1,int const m1, int const j2,int const m2, double & onecg);
		 
    std::vector<double> m_cg_tab;					// array of coefs
    std::vector<int> m_tab_idxj;					// array of offsets to find the subset(m,m1,m2) in the array
    int m_nt;							// input parameter <=> ntype
    int m_nt_jmax;						// jmax*ntype
    double m_jmax;						// input parameter jmax

    // tabulated n! values
    static double const m_factorial[168];
};
