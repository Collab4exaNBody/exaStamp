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

