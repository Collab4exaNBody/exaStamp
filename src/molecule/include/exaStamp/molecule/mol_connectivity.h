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
#include <vector>
#include <array>
#include <cstdint>

#include <onika/math/basic_types_def.h>
#include <exaStamp/atom_bond_connectivity.h>

#include <onika/memory/allocator.h>
#include <onika/oarray.h>


namespace exaStamp
{
  using namespace exanb;

  using ChemicalBond      = onika::oarray_t<uint64_t,2>;
  using ChemicalBonds     = onika::memory::CudaMMVector< ChemicalBond >;
  using ChemicalAngle     = onika::oarray_t<uint64_t,3>;
  using ChemicalAngles    = onika::memory::CudaMMVector< ChemicalAngle >;
  using ChemicalTorsion   = onika::oarray_t<uint64_t,4>;
  using ChemicalTorsions  = onika::memory::CudaMMVector< ChemicalTorsion >;
  using ChemicalImproper  = onika::oarray_t<uint64_t,4>;
  using ChemicalImpropers = onika::memory::CudaMMVector< ChemicalImproper >;
  using ChemicalPair      = onika::oarray_t<uint64_t,2>;
  using ChemicalPairs     = onika::memory::CudaMMVector< ChemicalPair >;

  using MoleculeConnectivity = exaStamp::AtomBondConnectivity;
} // namespace exaStamp

namespace YAML
{

  template<> struct convert<exaStamp::ChemicalBond>
  {
    static bool decode(const Node& node, exaStamp::ChemicalBond& value)
    {
      const char* chem = "Bond";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" YAML must be a sequence" << std::endl;
      const auto v = node.as< std::vector<long> >();
      if( v.size() != value.size() ) exanb::fatal_error() << chem << " sequence size must be "<<value.size() << std::endl;
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };

  template<> struct convert<exaStamp::ChemicalBonds>
  {
    static bool decode(const Node& node, exaStamp::ChemicalBonds& value)
    {
      const char* chem = "Bond";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" list must be a YAML sequence" << std::endl;
      const auto v = node.as< std::vector< exaStamp::ChemicalBond > >();
      value.resize( v.size() );
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };

  template<> struct convert<exaStamp::ChemicalAngle>
  {
    static bool decode(const Node& node, exaStamp::ChemicalAngle& value)
    {
      const char* chem = "Angle";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" YAML must be a sequence" << std::endl;
      const auto v = node.as< std::vector<long> >();
      if( v.size() != value.size() ) exanb::fatal_error() << chem << " sequence size must be "<<value.size() << std::endl;
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };

  template<> struct convert<exaStamp::ChemicalAngles>
  {
    static bool decode(const Node& node, exaStamp::ChemicalAngles& value)
    {
      const char* chem = "Angle";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" list must be a YAML sequence" << std::endl;
      const auto v = node.as< std::vector< exaStamp::ChemicalAngle > >();
      value.resize( v.size() );
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };

  template<> struct convert<exaStamp::ChemicalTorsion>
  {
    static bool decode(const Node& node, exaStamp::ChemicalTorsion& value)
    {
      const char* chem = "Torsion";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" YAML must be a sequence" << std::endl;
      const auto v = node.as< std::vector<long> >();
      if( v.size() != value.size() ) exanb::fatal_error() << chem << " sequence size must be "<<value.size() << std::endl;
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };

  template<> struct convert<exaStamp::ChemicalTorsions>
  {
    static bool decode(const Node& node, exaStamp::ChemicalTorsions& value)
    {
      const char* chem = "Torsion";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" list must be a YAML sequence" << std::endl;
      const auto v = node.as< std::vector< exaStamp::ChemicalTorsion > >();
      value.resize( v.size() );
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };

#if 0
  template<> struct convert<exaStamp::ChemicalImproper>
  {
    static bool decode(const Node& node, exaStamp::ChemicalImproper& value)
    {
      const char* chem = "Improper";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" YAML must be a sequence" << std::endl;
      const auto v = node.as< std::vector<long> >();
      if( v.size() != value.size() ) exanb::fatal_error() << chem << " sequence size must be "<<value.size() << std::endl;
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };

  template<> struct convert<exaStamp::ChemicalImpropers>
  {
    static bool decode(const Node& node, exaStamp::ChemicalImpropers& value)
    {
      const char* chem = "Improper";
      if( !node.IsSequence() ) exanb::fatal_error() << "Chemical "<< chem <<" list must be a YAML sequence" << std::endl;
      const auto v = node.as< std::vector< exaStamp::ChemicalImproper > >();
      value.resize( v.size() );
      for(size_t i=0;i<value.size();i++) value[i] = v[i];
      return true;
    }
  };
#endif

}

