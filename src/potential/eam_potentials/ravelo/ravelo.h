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

#include <exanb/core/spline.h>

namespace exaStamp
{

  struct EamRaveloParameters
  {
    double Ec;
    double alpha;
    double a0;
    double f3;
    double f4;
    double U0;
    double r1;
    double alphap;
    double beta3;
    double beta4;
    double rs;
    double s;
    double rc;
    double a1;
    double a2;
    double a3;
    double a4;
    double p;
    double q;
    double rho0;

    exanb::Spline m_f;
    exanb::Spline m_df;
  };

  static inline void eam_ravelo_phi(const EamRaveloParameters& p, double r, double& phi, double& dphi)
  {
    // assume r >= 0
    if(r <= p.rs) // <=> 0 <= r <= rs
      {
	const double alphap_r1 = p.alphap / p.r1;

	const double rstar  = r * alphap_r1 - p.alphap ;
	const double rstar2 = rstar    * rstar;
	const double rstar3 = rstar2   * rstar;
	const double rstar4 = rstar2   * rstar2;
	
	const double tmp = -p.U0 * std::exp(-rstar);

	phi = tmp * (1.0 + rstar + p.beta3 * rstar3 + p.beta4 * rstar4);
	dphi = tmp * alphap_r1 * (-rstar + (3 * rstar2 - rstar3) * p.beta3  + (4 * rstar3 - rstar4) * p.beta4);
      }
    else if(r <= p.rc) // <=> rs <= r <= rc
      {
	const double rcr  = p.rc - r;
	const double rcr2 = rcr * rcr;
	const double rcr3 = rcr * rcr2;
	
	const double tmp2 = p.a2 * rcr;
	const double tmp3 = p.a3 * rcr2;
	const double tmp4 = p.a4 * rcr3;

	phi = p.U0 * std::pow(rcr, p.s) * (p.a1 + tmp2 + tmp3 + tmp4);
	dphi = -p.U0 * std::pow(rcr, p.s - 1) * (p.a1 * p.s + tmp2 * (p.s+1) + tmp3 * (p.s+2) + tmp4 * (p.s+3));
      }
    else // <=> r > rc
      {
	phi = 0;
	dphi = 0;  
      }
  }

  static inline void eam_ravelo_rho(const EamRaveloParameters& p, double r, double& rho, double& drho)
  {
    // assume r >= 0
    if(r <= p.rc)
      {
	const double coef = std::sqrt(3.0) / 2.0;
	const double r0 = coef * p.a0 ; 

	const double num = std::pow(p.rc, p.p) - std::pow(r, p.p);
	const double den = std::pow(p.rc, p.p) - std::pow(r0, p.p);	
	
	const double tmp = num / den;

	rho = p.rho0 * std::pow(tmp, p.q);
	drho = -p.rho0 * p.q * p.p * std::pow(r, p.p-1) * std::pow(tmp, p.q-1) / den;
      }
    else
      {
	rho = 0;
	drho = 0;
      }
  }

  static inline void eam_ravelo_fEmbed(const EamRaveloParameters& p, double rho, double& f, double& df)
  {
    f = p.m_f.eval(rho);
    df = p.m_df.eval(rho);
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EamRaveloParameters;
  using onika::physics::Quantity;

  template<> struct convert<EamRaveloParameters>
    {
      static bool decode(const Node& node, EamRaveloParameters& v);
	/* { */
	/* 	if( !node.IsMap() ) { return false; } */
	/* #     define CONVERT_EAMRAVELO_PARAM(x) v.x = node[#x].as<Quantity>().convert() //; std::cout << #x << " = " << v.x << std::endl  */
	/* 	CONVERT_EAMRAVELO_PARAM(Ec); */
	/* 	CONVERT_EAMRAVELO_PARAM(alpha); */
	/* 	CONVERT_EAMRAVELO_PARAM(a0); */
	/* 	CONVERT_EAMRAVELO_PARAM(f3); */
	/* 	CONVERT_EAMRAVELO_PARAM(f4); */
	/* 	CONVERT_EAMRAVELO_PARAM(U0); */
	/* 	CONVERT_EAMRAVELO_PARAM(r1); */
	/* 	CONVERT_EAMRAVELO_PARAM(alphap); */
	/* 	CONVERT_EAMRAVELO_PARAM(beta3); */
	/* 	CONVERT_EAMRAVELO_PARAM(beta4); */
	/* 	CONVERT_EAMRAVELO_PARAM(rs); */
	/* 	CONVERT_EAMRAVELO_PARAM(s); */
	/* 	CONVERT_EAMRAVELO_PARAM(rc); */
	/* 	CONVERT_EAMRAVELO_PARAM(a1); */
	/* 	CONVERT_EAMRAVELO_PARAM(a2); */
	/* 	CONVERT_EAMRAVELO_PARAM(a3); */
	/* 	CONVERT_EAMRAVELO_PARAM(a4); */
	/* 	CONVERT_EAMRAVELO_PARAM(p); */
	/* 	CONVERT_EAMRAVELO_PARAM(q); */
	/* 	CONVERT_EAMRAVELO_PARAM(rho0); */
	/* #     undef CONVERT_EAMRAVELO_PARAM */
	/* 	return true; */
	/*       } */
    };
}

