#pragma once

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <exanb/core/unityConverterHelper.h>
//#include <exanb/pair_potential.h>
#include <exanb/core/physics_constants.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // ReactionField Parameters
  struct ReactionFieldParms
  {
    //double epsilon = 80.0;
    double rc = 12.5; 
    double RF0 = 0.0;
    double RF1 = 0.0;
    double RF2 = 0.0;
    double ecut = 0.0;
    
    ONIKA_HOST_DEVICE_FUNC
    inline bool is_null() const { return RF0==0.0 && RF1==0.0 && RF2==0.0 && ecut==0.0; }    
  };

  // core computation kernel for reaction field potential
  ONIKA_HOST_DEVICE_FUNC inline void reaction_field_compute_energy(const ReactionFieldParms& p_rc, double c /* =c1*c2 */, double r, double& e, double& de)
  {
    assert( r > 0. );
    const double r2 = r * r;
    // double c  = c1 * c2;
    e  = ( (                 p_rc.RF0/r  - p_rc.RF2                   + p_rc.RF1*r2 ) - p_rc.ecut ) * c;
    de =   ( - p_rc.RF0/r2                          + 2.*p_rc.RF1 * r               ) * c;
  }

  inline void init_rf(ReactionFieldParms& v, double rc, double epsilon)
  {
    const double one_FourPiEpsilon0 = 1./4./M_PI/UnityConverterHelper::convert(exanb::legacy_constant::epsilonZero, "C^2.s^2/m^3/kg^1");
    v.rc = rc;
    if( rc==0.0 || epsilon==0.0 )
    {
      v.rc = 1e-18;
      v.RF0 = 0.0;
      v.RF1 = 0.0;
      v.RF2 = 0.0;
      v.ecut = 0.0;
    }
    else
    {
      v.RF0     = one_FourPiEpsilon0;
      v.RF1     = one_FourPiEpsilon0 / pow(v.rc,3) * (epsilon - 1.)/(2.*epsilon + 1.);
      v.RF2     = one_FourPiEpsilon0 / v.rc * (3. * epsilon)/(2.*epsilon + 1.);
      v.ecut    = 0.0;
      double e=0.0, de=0.0;
      reaction_field_compute_energy( v , 1.0 , v.rc , e , de );
      v.ecut = e;
    }      
  }

}


// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::ReactionFieldParms;
  using exaStamp::reaction_field_compute_energy;
  using exanb::UnityConverterHelper;
  using exanb::Quantity;

  template<> struct convert<ReactionFieldParms>
  {
    static bool decode(const Node& node, ReactionFieldParms& v)
    {
      if( !node.IsMap() ) { return false; }
      init_rf( v, node["rc"].as<Quantity>().convert() , node["epsilon"].as<double>() );
      return true;
    }
  };
}

