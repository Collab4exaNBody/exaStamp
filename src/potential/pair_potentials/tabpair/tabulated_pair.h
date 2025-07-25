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

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exanb/core/ro_spline.h>

namespace exaStamp
{
  using namespace exanb;

  // Tabulated Parameters
  struct TabPairPotentialParms
  {
    exanb::Spline m_e;
    exanb::Spline m_de;
  };

  // read only Tabulated Parameters
  struct ReadOnlyTabPairPotentialParms
  {
    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_e;
    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_de;

    ReadOnlyTabPairPotentialParms() = default;
    ReadOnlyTabPairPotentialParms(const ReadOnlyTabPairPotentialParms&) = default;
    ReadOnlyTabPairPotentialParms(ReadOnlyTabPairPotentialParms&&) = default;
    ReadOnlyTabPairPotentialParms& operator = (const ReadOnlyTabPairPotentialParms&) = default;
    ReadOnlyTabPairPotentialParms& operator = (ReadOnlyTabPairPotentialParms&&) = default;
    
    inline ReadOnlyTabPairPotentialParms( const TabPairPotentialParms& p )
      : m_e( p.m_e )
      , m_de( p.m_de )
      {}
  };

}

// specialize ReadOnlyShallowCopyType so ReadOnlyEwaldParms is the read only type for EwaldParms
namespace onika { namespace cuda {
  template<> struct ReadOnlyShallowCopyType< exaStamp::TabPairPotentialParms > { using type = exaStamp::ReadOnlyTabPairPotentialParms; };
} }


namespace exaStamp
{

  ONIKA_HOST_DEVICE_FUNC inline void tab_pair_potential_energy(const ReadOnlyTabPairPotentialParms& p, const PairPotentialMinimalParameters&, double r, double& e, double& de)
  {
    e = p.m_e.eval(r);
    de = p.m_de.eval(r);
  }

  inline void tab_pair_potential_energy(const TabPairPotentialParms& p, const PairPotentialMinimalParameters&, double r, double& e, double& de)
  {
    e = p.m_e.eval(r);
    de = p.m_de.eval(r);
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::TabPairPotentialParms>
  {
    static bool decode(const Node& node, exaStamp::TabPairPotentialParms& v);
  };
}

