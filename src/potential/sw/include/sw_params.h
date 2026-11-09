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

#include <yaml-cpp/yaml.h>
#include <onika/memory/allocator.h>
#include <string>

namespace exaStamp
{

  // User parameters (YAML) of the sw_init operator
  struct SWParams
  {
    std::string file;       // .sw parameter file (element triplets, energies in eV, lengths in Angstrom)
    bool threebody = true;  // false: drop the three-body term (pure pair potential)
  };

  // One element-triplet record (i,j,k), with derived constants precomputed.
  // Energy-dimensioned constants are already converted to exaStamp internal units.
  // Two-body and radial three-body terms of pair (i,j) use record (i,j,j),
  // the angular term of triplet j-i-k uses record (i,j,k).
  struct SWTriplet
  {
    // raw parameters
    double epsilon = 0.0, sigma = 0.0, littlea = 0.0, lambda = 0.0, gamma = 0.0, costheta = 0.0;
    double biga = 0.0, bigb = 0.0, powerp = 0.0, powerq = 0.0, tol = 0.0;
    // derived
    double cut = 0.0, cutsq = 0.0;
    double sigma_gamma = 0.0, lambda_epsilon = 0.0, lambda_epsilon2 = 0.0;
    double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0, c6 = 0.0;
  };

  // Output of sw_init, consumed by sw_force.
  // Triplet table indexed by exaStamp species index : (ti*ntypes + tj)*ntypes + tk
  struct SWContext
  {
    onika::memory::CudaMMVector<SWTriplet> triplets;
    int ntypes = 0;
    double rcut = 0.0;      // max SW cutoff over all triplets
    bool threebody = true;
  };

}

namespace YAML
{
  template<> struct convert<exaStamp::SWParams>
  {
    static bool decode(const Node& node, exaStamp::SWParams& v)
    {
      if( !node.IsMap() )   { return false; }
      if( !node["file"] )   { return false; }
      v.file = node["file"].as<std::string>();
      if( node["threebody"] ) { v.threebody = node["threebody"].as<bool>(); }
      return true;
    }
  };
}
