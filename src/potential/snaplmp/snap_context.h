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

namespace exaStamp
{

  using namespace exanb;
//  using onika::memory::DEFAULT_ALIGNMENT;
//  using namespace SnapExt;

  struct SnapLMPThreadContext
  {
    LAMMPS_NS::SNA * sna = nullptr;
    double * beta = nullptr;
    double * bispectrum = nullptr;
  };

  struct SnapLMPContext
  {
    LAMMPS_NS::LAMMPS * ptr = nullptr;
    std::vector<SnapLMPThreadContext> m_thread_ctx;
    SnapExt::SnapConfig m_config;
    std::vector<double> m_coefs;    // size = number of bispectrum coefficients
    std::vector<double> m_factor;   // size = number of elements
    std::vector<double> m_radelem;  // size = number of elements
    
    std::vector<double> m_beta;       // size = number of particles
    std::vector<double> m_bispectrum; // size = number of particles
    double m_rcut = 0.0;
  };

}


