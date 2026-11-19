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
#include <string>
#include <yaml-cpp/yaml.h>
#include <onika/log.h>
#include <exaStamp/yaml_check_keys.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/cuda/cuda.h>
#include <exaStamp/unit_system.h>

namespace exaStamp
{
  using namespace exanb;

  // Treatment of the force at the reaction field cutoff rc (the energy is always zero at rc): truncated (default),
  // shifted force E(r) - (r-rc) E'(rc), or switching S(x) E(r) with a quintic S between rs and rc.
  static constexpr int RF_CUTOFF_TRUNCATED = 0;
  static constexpr int RF_CUTOFF_SHIFTED_FORCE = 1;
  static constexpr int RF_CUTOFF_SWITCHING = 2;

  // Reaction field coulomb potential (dielectric continuum beyond rc). Single definition, used by the coul_rf pair style,
  // the ljrf, exp6rf and ljexp6rf pair styles, the ljexp6rf_pc molecular force field and coulombic_rf (per atom charges).
  struct ReactionFieldParms
  {
    //double epsilon = 80.0;
    double rc = 0.0; 
    double RF0 = 0.0;
    double RF1 = 0.0;
    double RF2 = 0.0;
    double ecut = 0.0;
    int cutoff_mode = RF_CUTOFF_TRUNCATED;
    double rs = 0.0;                // switching start (RF_CUTOFF_SWITCHING)
    double inv_rc2 = 0.0;           // 1/rc^2 (RF_CUTOFF_SHIFTED_FORCE)
    double inv_switch_width = 0.0;  // 1/(rc-rs) (RF_CUTOFF_SWITCHING)
    
    ONIKA_HOST_DEVICE_FUNC
    inline bool is_null() const { return RF0==0.0 && RF1==0.0 && RF2==0.0; }    
  };

  // Cutoff treatment, as given in input files
  struct ReactionFieldCutoffMode
  {
    int mode = RF_CUTOFF_TRUNCATED;
    double rs = 0.0;
  };

  inline const char* reaction_field_cutoff_mode_name(int mode)
  {
    switch( mode )
    {
      case RF_CUTOFF_SHIFTED_FORCE : return "shifted_force";
      case RF_CUTOFF_SWITCHING : return "switching";
      default : return "truncated";
    }
  }

  inline int reaction_field_cutoff_mode_from_name(const std::string& name)
  {
    if( name == "truncated" ) return RF_CUTOFF_TRUNCATED;
    if( name == "shifted_force" ) return RF_CUTOFF_SHIFTED_FORCE;
    if( name == "switching" ) return RF_CUTOFF_SWITCHING;
    onika::fatal_error() << "reaction field: unknown cutoff mode '" << name << "' (expected truncated, shifted_force or switching)" << std::endl;
    return RF_CUTOFF_TRUNCATED;
  }

  // Sets the cutoff treatment of initialized parameters (rc, RF0, RF1, RF2). Without reaction field (null parameters),
  // there is nothing to modify.
  inline void reaction_field_set_cutoff_mode(ReactionFieldParms& v, const ReactionFieldCutoffMode& m)
  {
    v.cutoff_mode = RF_CUTOFF_TRUNCATED;
    v.rs = 0.0;
    v.inv_rc2 = 0.0;
    v.inv_switch_width = 0.0;
    if( v.is_null() ) return;
    if( m.mode == RF_CUTOFF_SHIFTED_FORCE )
    {
      v.cutoff_mode = RF_CUTOFF_SHIFTED_FORCE;
      v.inv_rc2 = 1.0 / ( v.rc * v.rc );
    }
    else if( m.mode == RF_CUTOFF_SWITCHING )
    {
      if( !( m.rs > 0.0 && m.rs < v.rc ) )
      {
        onika::fatal_error() << "reaction field: switching start rs="<<m.rs<<" must be between 0 and rc="<<v.rc << std::endl;
      }
      v.cutoff_mode = RF_CUTOFF_SWITCHING;
      v.rs = m.rs;
      v.inv_switch_width = 1.0 / ( v.rc - m.rs );
    }
  }

  // Applies the cutoff treatment to the reaction field energy e and its derivative de = dE/dr, for r < rc.
  // weight scales the 1/r term, as in reaction_field_compute_energy_weighted.
  ONIKA_HOST_DEVICE_FUNC
  inline void reaction_field_apply_cutoff_mode(const ReactionFieldParms& p_rc, double c, double r, double& e, double& de, double weight = 1.0)
  {
    if( p_rc.cutoff_mode == RF_CUTOFF_SHIFTED_FORCE )
    {
      const double de_rc = c * ( 2. * p_rc.RF1 * p_rc.rc - weight * p_rc.RF0 * p_rc.inv_rc2 );
      e  -= ( r - p_rc.rc ) * de_rc;
      de -= de_rc;
    }
    else if( p_rc.cutoff_mode == RF_CUTOFF_SWITCHING && r > p_rc.rs )
    {
      const double x = ( r - p_rc.rs ) * p_rc.inv_switch_width;
      const double x2 = x * x;
      const double s = 1. - x * x2 * ( 10. - 15. * x + 6. * x2 );
      const double ds = -30. * x2 * ( 1. - x ) * ( 1. - x ) * p_rc.inv_switch_width;
      de = de * s + e * ds;
      e  = e * s;
    }
  }

  // core computation kernel for reaction field potential
  ONIKA_HOST_DEVICE_FUNC
  inline void reaction_field_compute_energy(const ReactionFieldParms& p_rc, double c /* =c1*c2 */, double r, double& e, double& de)
  {
    assert( r > 0. );
    const double r2 = r * r;
    // double c  = c1 * c2;

    double e_coul  =   c * p_rc.RF0/r;
    double de_coul = - c * p_rc.RF0/r2;

    double e_corr  =   c * ( p_rc.RF1*r2  - p_rc.RF2 );
    double de_corr =   c * 2. * p_rc.RF1 * r;
    
    e  =  e_coul +  e_corr;
    de = de_coul + de_corr;
    reaction_field_apply_cutoff_mode( p_rc, c, r, e, de );
  }

  // core computation kernel for reaction field potential
  ONIKA_HOST_DEVICE_FUNC
  inline void reaction_field_compute_energy_weighted(const ReactionFieldParms& p_rc, double c /* =c1*c2 */, double r, double& e, double& de, double weight)
  {
    assert( r > 0. );
    const double r2 = r * r;
    // double c  = c1 * c2;
    // e  = ( (                 p_rc.RF0/r  - p_rc.RF2                   + p_rc.RF1*r2 ) - p_rc.ecut ) * c;
    // de =   ( - p_rc.RF0/r2                          + 2.*p_rc.RF1 * r               ) * c;

    double e_coul  =   c * weight * p_rc.RF0/r;
    double de_coul = - c * weight * p_rc.RF0/r2;

    double e_corr  = c * ( p_rc.RF1*r2  - p_rc.RF2 );
    double de_corr = c * 2. * p_rc.RF1 * r;
    
    e  =  e_coul +  e_corr;
    de = de_coul + de_corr;
    reaction_field_apply_cutoff_mode( p_rc, c, r, e, de, weight );
  }

  // core computation kernel for reaction field potential
  ONIKA_HOST_DEVICE_FUNC
  inline void reaction_field_compute_self_energy(const ReactionFieldParms& p_rc, double c, double& e)
  {
    e  = - c * p_rc.RF2;
  }
  
  inline void init_rf(ReactionFieldParms& v, double rc, double epsilon)
  {
    static constexpr double Epsilon0 = EXASTAMP_CONST_QUANTITY( onika::physics::epsilonZero * ( C^2 ) * ( s^2 ) / ( m^3 ) / ( kg^1 ) );
    static constexpr double one_FourPiEpsilon0 = 1. / 4. / M_PI / Epsilon0;
    v.rc = rc;
    if( rc==0.0 || epsilon==0.0 )
    {
      v.rc = 1e-18;
      v.RF0 = 0.0;
      v.RF1 = 0.0;
      v.RF2 = 0.0;
    }
    else
    {
      v.RF0     = one_FourPiEpsilon0;
      v.RF1     = one_FourPiEpsilon0 / pow(v.rc,3) * (epsilon - 1.)/(2.*epsilon + 1.);
      v.RF2     = one_FourPiEpsilon0 / v.rc * (3. * epsilon)/(2.*epsilon + 1.);
    }
    // new rc and coefficients : back to the plain truncated reaction field
    reaction_field_set_cutoff_mode( v , ReactionFieldCutoffMode{} );
  }

}


// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::ReactionFieldParms;
  using exaStamp::reaction_field_compute_energy;
  using onika::physics::Quantity;

  // cutoff_mode: truncated | shifted_force | { mode: switching, rs: <distance> }
  template<> struct convert<exaStamp::ReactionFieldCutoffMode>
  {
    static bool decode(const Node& node, exaStamp::ReactionFieldCutoffMode& v)
    {
      v = exaStamp::ReactionFieldCutoffMode{};
      if( node.IsScalar() )
      {
        v.mode = exaStamp::reaction_field_cutoff_mode_from_name( node.as<std::string>() );
      }
      else if( node.IsMap() )
      {
        exaStamp::check_yaml_keys( node , { "mode", "rs" } , "reaction field cutoff_mode" );
        if( !node["mode"] ) { onika::fatal_error() << "reaction field cutoff_mode: missing key 'mode'" << std::endl; }
        v.mode = exaStamp::reaction_field_cutoff_mode_from_name( node["mode"].as<std::string>() );
        if( node["rs"] ) v.rs = node["rs"].as<Quantity>().convert();
      }
      else { return false; }
      if( v.mode == exaStamp::RF_CUTOFF_SWITCHING && v.rs <= 0.0 )
      {
        onika::fatal_error() << "reaction field cutoff_mode: switching needs the switching start rs, e.g. { mode: switching, rs: 11.5 ang }" << std::endl;
      }
      if( v.mode != exaStamp::RF_CUTOFF_SWITCHING && v.rs != 0.0 )
      {
        onika::fatal_error() << "reaction field cutoff_mode: rs is only used with mode: switching" << std::endl;
      }
      return true;
    }
  };

  template<> struct convert<ReactionFieldParms>
  {
    static bool decode(const Node& node, ReactionFieldParms& v)
    {
      if( !node.IsMap() ) { return false; }
      exaStamp::check_yaml_keys( node , { "rc", "epsilon", "cutoff_mode" } , "reaction field parameters" );
      init_rf( v, node["rc"].as<Quantity>().convert() , node["epsilon"].as<double>() );
      if( node["cutoff_mode"] ) reaction_field_set_cutoff_mode( v , node["cutoff_mode"].as<exaStamp::ReactionFieldCutoffMode>() );
      return true;
    }
  };
}

