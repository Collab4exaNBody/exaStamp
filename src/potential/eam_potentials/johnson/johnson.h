#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/units.h>

namespace exaStamp
{
  using namespace exanb;

  struct EamJohnsonParameters
  {
    double re;
    double fe;
    double rhoe;
    double alpha;
    double beta;
    double A;
    double B;
    double kappa;
    double lambda;
    double Fn0;
    double Fn1;
    double Fn2;
    double Fn3;
    double F0;
    double F1;
    double F2;
    double F3;
    double Fo;
    double eta;
  };

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC inline void eam_johnson_phi(const EamJohnsonParameters& p, double r, double& phi, double& dphi)
  {
      static constexpr int n = 20;
      static constexpr int m = 20;
    
      double x = r / p.re;
      double ire = 1./p.re;
      
      /// Premier terme 
      double c1 = p.A * exp(-p.alpha*(x-1.));
      double c2 = x-p.kappa;
      double c3 = pow(c2,m);

      double num = c1;
      double den = 1. + c3;

      double dnum = -p.alpha * c1;
      double dden = m * c3 / c2;
          
      phi = num / den;
      dphi = ire * (dnum * den - num * dden)/(den * den);

      /// Second terme  
      c1 = -p.B * exp(-p.beta*(x-1.));
      c2 = x-p.lambda;
      c3 = pow(c2,n);

      num = c1;
      den = 1. + c3;

      dnum = -p.beta * c1;
      dden = n * c3 / c2;   

      phi += num / den;
      dphi += ire * (dnum * den - num * dden)/(den * den);  
  }

  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  ONIKA_HOST_DEVICE_FUNC inline void eam_johnson_rho(const EamJohnsonParameters& p, double r, double& rho, double& drho)
  {
      static constexpr int n = 20;
    	
      double x = r / p.re;
      double ire = 1./p.re;
	    double c1 = p.fe * exp(-p.beta*(x-1.));
      double c2 = x-p.lambda;
      double c3 = pow(c2,n);
      
      double num = c1;
      double den = 1. + c3;
      
      double dnum = -p.beta * c1 ;
      double dden = n * c3 / c2;
      
      rho = num / den;
      drho = ire * (dnum * den - num * dden)/(den * den);
  }

  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  ONIKA_HOST_DEVICE_FUNC inline void eam_johnson_fEmbed(const EamJohnsonParameters& p, double rho, double& f, double& df)
  {
      double rhon = 0.85 * p.rhoe;
      double rho0 = 1.15 * p.rhoe;
      double irhon = 1. / rhon;
      double irhoe = 1. / p.rhoe;
      double raprho = 0.;
      double puis[4];
      
      puis[0] = 1.;
         
      if(rho < rhon)
      {
		  raprho = rho / rhon; 
		  puis[1] = raprho - 1.;
		  puis[2] = puis[1] * puis [1];		
		  puis[3] = puis[1] * puis [2];	
		  
		  f = p.Fn0 + p.Fn1 * puis[1] + p.Fn2 * puis[2] + p.Fn3 * puis[3] ; 	
		  df = p.Fn1 + 2. * p.Fn2 * puis[1] + 3. * p.Fn3 * puis[2];
		  df *= irhon;
				   		   
      }
      else if(rho < rho0)
      {
		  raprho = rho * irhoe;     
		  puis[1] = raprho - 1.;
		  puis[2] = puis[1] * puis [1];		
		  puis[3] = puis[1] * puis [2];	
		  
		  f = p.F0 + p.F1 * puis[1] + p.F2 * puis[2] + p.F3 * puis[3] ; 	
		  df = p.F1 + 2. * p.F2 * puis[1] + 3. * p.F3 * puis[2];
		  df *= irhoe;		
		  
      }
      else
      {
		  raprho = rho * irhoe;     
		  double raprhopuiseta= pow(raprho,p.eta);
		  double lograprhopuiseta = log(raprhopuiseta);
		  
		  f = p.Fo * (1. - lograprhopuiseta) * raprhopuiseta;  
		  df = -p.eta * raprhopuiseta + (1. - lograprhopuiseta) * p.eta * raprhopuiseta;
		  df *= p.Fo * irhoe / raprho;
      }
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EamJohnsonParameters;
  using onika::physics::Quantity;

  template<> struct convert<EamJohnsonParameters>
  {
    static bool decode(const Node& node, EamJohnsonParameters& v)
    {
      if( !node.IsMap() ) { return false; }
#     define CONVERT_EAMJOHNSON_PARAM(x) v.x = node[#x].as<Quantity>().convert() //; std::cout << #x << " = " << v.x << std::endl 
      CONVERT_EAMJOHNSON_PARAM(re);
      CONVERT_EAMJOHNSON_PARAM(fe);
      CONVERT_EAMJOHNSON_PARAM(rhoe);
      CONVERT_EAMJOHNSON_PARAM(alpha);
      CONVERT_EAMJOHNSON_PARAM(beta);
      CONVERT_EAMJOHNSON_PARAM(A);
      CONVERT_EAMJOHNSON_PARAM(B);
      CONVERT_EAMJOHNSON_PARAM(kappa);
      CONVERT_EAMJOHNSON_PARAM(lambda);
      CONVERT_EAMJOHNSON_PARAM(Fn0);
      CONVERT_EAMJOHNSON_PARAM(Fn1);
      CONVERT_EAMJOHNSON_PARAM(Fn2);
      CONVERT_EAMJOHNSON_PARAM(Fn3);
      CONVERT_EAMJOHNSON_PARAM(F0);
      CONVERT_EAMJOHNSON_PARAM(F1);
      CONVERT_EAMJOHNSON_PARAM(F2);
      CONVERT_EAMJOHNSON_PARAM(F3);
      CONVERT_EAMJOHNSON_PARAM(Fo);
      CONVERT_EAMJOHNSON_PARAM(eta);
#     undef CONVERT_EAMJOHNSON_PARAM
      return true;
    }
  };
}

