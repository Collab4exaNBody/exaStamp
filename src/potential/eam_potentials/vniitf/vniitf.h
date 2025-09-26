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

/// @file 
/// @brief Definition of EAM VNIITF potential
/// This is a EAM potential developed to calculate the properties of tin near melting curve
/// @see http://iopscience.iop.org/article/10.1088/1742-6596/500/3/032017/pdf


#include <cmath>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/units.h>
#include <onika/cuda/cuda.h>

namespace exaStamp
{

  /// @brief Structure to handle parameters of a EAM VNIITF potential
  struct EamVniitfParameters 
  {
    double rmax; ///< Maximal distance for density contribution
    double rmin; ///< Minimal distance for density contribution
    double rt0; ///< Characteristic radius
    double Ecoh; ///< Cohesive energy
    double E0; ///< Energy parameter
    double beta; ///< Attenuation rate
    double A; ///< Empirical parameter A
    double Z; ///< Number of neighbors in the reference structure
    double n; ///< Empirical parameter n
    double alpha; ///< Parameter alpha
    double D; ///< Parameter D
    double eta; ///< Parameter eta
    double mu; ///< Parameter mu
  };

  /// @brief Derivative of the smoothstep 3
  /// @param [in] x Input (only values between 0 and 1 are relevant)
  /// @return Output
  ONIKA_HOST_DEVICE_FUNC static inline double Spline_dS3(double x)
  {
    double x2 = x *x;
    double x3 = x2*x;
    return 140*x3 * ( -1*x3 + 3*x2 - 3*x + 1 );
  }

  /// @brief Smoothstep 3 : smootheststep
  /// @param [in] x Input (only values between 0 and 1 are relevant)
  /// @return Output
  ONIKA_HOST_DEVICE_FUNC static inline double Spline_S3(double x)
  {
    double x2 = x*x;
    return x2*x2 * ( -20*x2*x + 70*x2 - 84*x + 35 );
  }

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC static inline void eam_vniitf_phi(const EamVniitfParameters& p, double r, double& phi, double& dphi)
  {
    double ir   = 1/r;
    double irt0 = 1/p.rt0;
    double dr   = r*irt0 - 1.0;
    double dr2  = dr*dr;
    double a    = -2*p.Ecoh/p.Z;
    double b    = p.alpha*p.alpha*p.alpha*p.D*p.rt0;

    double f1  = a*( 1 + p.alpha*dr + p.eta*dr2 + (p.mu+b*ir)*dr2*dr );
    double df1 = a*( p.alpha*irt0 + 2*p.eta*irt0*dr + 3*p.mu*irt0*dr2 + b*(3*irt0-dr*ir)*dr2*ir );

    double f2  = std::exp(-p.alpha*dr);
    double df2 = -p.alpha*irt0*f2;

    double drS = (p.rmax-r)/(p.rmax-p.rmin);
    double S  = Spline_S3(drS);
    double dS = Spline_dS3(drS) / (p.rmin-p.rmax);

    // Maybe this could be made before
    if (drS<0) {
      S  = 0.;
      dS = 0.;
    }
    else if (drS>1) {
      S  = 1.;
      dS = 0.;
    }

    phi  = (p.E0 + f1*f2) *  S;
    dphi = (p.E0 + f1*f2) * dS + (f1*df2 + f2*df1) * S;
  }

  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC static inline void eam_vniitf_rho(const EamVniitfParameters& p, double r, double& rho, double& drho)
  {
    double irt0 = 1/p.rt0;

    double F  = exp(-p.beta * (r*irt0 - 1.0) ) / p.Z;
    double dF = -p.beta*F*irt0;

    double drS = (p.rmax-r)/(p.rmax-p.rmin);
    double S  = Spline_S3(drS);
    double dS = Spline_dS3(drS) / (p.rmin-p.rmax);

    // Maybe this could be made before
    if (drS<0) {
      S  = 0.;
      dS = 0.;
    }
    else if (drS>1) {
      S  = 1.;
      dS = 0.;
    }

    rho  = F*S;
    drho = F*dS + S*dF;
  }

  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  ONIKA_HOST_DEVICE_FUNC static inline void eam_vniitf_fEmbed(const EamVniitfParameters& p, double rho, double& f, double& df)
  {
    if (rho<=0.) {

      f  = 0.;
      df = 0.;

    }
    else {

      double a = pow(rho, p.n);
      double b = p.A*p.Ecoh * a;
      double c = log(a);
      
      f  = b * (c-1);
      df = p.n*b*c / rho; 

    }
  }

}


// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EamVniitfParameters;
  using onika::physics::Quantity;

  template<> struct convert<EamVniitfParameters>
  {
    static bool decode(const Node& node, EamVniitfParameters& v)
    {
      if( !node.IsMap() ) { return false; }
#     define CONVERT_EAMJVNIITF_PARAM(x) v.x = node[#x].as<Quantity>().convert()
      CONVERT_EAMJVNIITF_PARAM( rmax );
      CONVERT_EAMJVNIITF_PARAM( rmin );
      CONVERT_EAMJVNIITF_PARAM( rt0 );
      CONVERT_EAMJVNIITF_PARAM( Ecoh );
      CONVERT_EAMJVNIITF_PARAM( E0 );
      CONVERT_EAMJVNIITF_PARAM( beta );
      CONVERT_EAMJVNIITF_PARAM( A );
      CONVERT_EAMJVNIITF_PARAM( Z );
      CONVERT_EAMJVNIITF_PARAM( n );
      CONVERT_EAMJVNIITF_PARAM( alpha );
      CONVERT_EAMJVNIITF_PARAM( D );
      CONVERT_EAMJVNIITF_PARAM( eta );
      CONVERT_EAMJVNIITF_PARAM( mu );
#     undef CONVERT_EAMJVNIITF_PARAM
      return true;
    }
  };
}

