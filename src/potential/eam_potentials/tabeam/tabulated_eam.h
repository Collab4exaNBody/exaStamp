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
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <exanb/core/ro_spline.h>

namespace exaStamp
{
  using namespace exanb;

  // LennardJones Parameters
  struct TabEAMPotentialParms
  {
    exanb::Spline m_phi;
    exanb::Spline m_dphi;

    exanb::Spline m_rho;
    exanb::Spline m_drho;

    exanb::Spline m_f;
    exanb::Spline m_df;
  };

  // read only Tabulated Parameters
  struct ReadOnlyTabEAMPotentialParms
  {
    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_phi;
    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_dphi;

    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_rho;
    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_drho;

    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_f;
    onika::cuda::ro_shallow_copy_t<exanb::Spline> m_df;

    ReadOnlyTabEAMPotentialParms() = default;
    ReadOnlyTabEAMPotentialParms(const ReadOnlyTabEAMPotentialParms&) = default;
    ReadOnlyTabEAMPotentialParms(ReadOnlyTabEAMPotentialParms&&) = default;
    ReadOnlyTabEAMPotentialParms& operator = (const ReadOnlyTabEAMPotentialParms&) = default;
    ReadOnlyTabEAMPotentialParms& operator = (ReadOnlyTabEAMPotentialParms&&) = default;
    
    inline ReadOnlyTabEAMPotentialParms( const TabEAMPotentialParms& p )
      : m_phi( p.m_phi )
      , m_dphi( p.m_dphi )
      , m_rho( p.m_rho )
      , m_drho( p.m_drho )
      , m_f( p.m_f )
      , m_df( p.m_df )
      {}
  };
}

// specialize ReadOnlyShallowCopyType so ReadOnlyEwaldParms is the read only type for EwaldParms
namespace onika { namespace cuda {
  template<> struct ReadOnlyShallowCopyType< exaStamp::TabEAMPotentialParms > { using type = exaStamp::ReadOnlyTabEAMPotentialParms; };
} }

namespace exaStamp
{
  template<class TabEAMPotentialParmsT>
  static inline void tab_eam_phi(const TabEAMPotentialParmsT& p, double r, double& phi, double& dphi)
  {
    phi = p.m_phi.eval(r);
    dphi = p.m_dphi.eval(r);
  }

  template<class TabEAMPotentialParmsT>
  static inline void tab_eam_rho(const TabEAMPotentialParmsT& p, double r, double& rho, double& drho)
  {
    rho = p.m_rho.eval(r);
    drho = p.m_drho.eval(r);
  }
  
  template<class TabEAMPotentialParmsT>
  static inline void tab_eam_fEmbed(const TabEAMPotentialParmsT& p, double rho, double& f, double& df)
  {
    f = p.m_f.eval(rho);
    df = p.m_df.eval(rho);
  }
  
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::TabEAMPotentialParms;

  template<> struct convert<TabEAMPotentialParms>
  {
    static bool decode(const Node& node, TabEAMPotentialParms& v);
  };
}

