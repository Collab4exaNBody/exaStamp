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
#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>
#include <exaStamp/potential/reaction_field/reaction_field.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // assembled Lennard-Jones & Reaction Field parameters
  struct LJRFParms
  {
    double lj_rcut = 0.0;
    double lj_ecut = 0.0;
    LennardJonesParms lj;
    ReactionFieldParms rf;
    bool shiftlj = true;
  };

//# pragma omp declare simd uniform(p,ppp) notinbranch
  ONIKA_HOST_DEVICE_FUNC inline void ljrf_energy(const LJRFParms& p, const PairPotentialMinimalParameters& pp, double r, double& e, double& de)
  {
    double lj_e=0.0, lj_de=0.0;
    if( r <= p.lj_rcut )
    {
      lj_compute_energy( p.lj , pp , r , lj_e , lj_de );
      lj_e -= p.lj_ecut;
      //printf("LJ: r=%g, rcut=%g, e=%g, de=%g\n",r,p.lj_rcut,lj_e,lj_de);
    }

    double rf_e=0.0, rf_de=0.0;
    if( r <= p.rf.rc )
    {
      reaction_field_compute_energy( p.rf, pp.m_atom_a.m_charge * pp.m_atom_b.m_charge, r, rf_e, rf_de );
      //printf("RF: r=%g, rcut=%g, e=%g, de=%g\n",r,p.rf.rc,rf_e,rf_de);
    }
    
    e = lj_e + rf_e;
    de = lj_de + rf_de;

    //printf("r=% .6e : e=% .6e , de=% .6e : c1=% .6e , c2=% .6e\n",r,e,de, pp.m_atom_a.m_charge, pp.m_atom_b.m_charge);
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::LJRFParms;
  using exaStamp::LennardJonesParms;
  using exaStamp::ReactionFieldParms;
  using onika::physics::Quantity;
  using exaStamp::lj_compute_energy;
  using exaStamp::PairPotentialMinimalParameters;

  template<> struct convert<LJRFParms>
  {
    static bool decode(const Node& node, LJRFParms& v)
    {
      if( !node.IsMap() ) { return false; }
      
      if( node["lj_rcut"] ) v.lj_rcut = node["lj_rcut"].as<Quantity>().convert();
      else v.lj_rcut = 0.0;
      
      if( node["lj"] ) v.lj = node["lj"].as<LennardJonesParms>();
      else v.lj = LennardJonesParms{ 0. , 0. };
      
      if( node["shiftlj"] ) v.shiftlj = node["shiftlj"].as<bool>();
      
      if (v.shiftlj) {
	double lj_e=0.0, lj_de=0.0;	
	if( v.lj_rcut > 0 ) lj_compute_energy( v.lj , PairPotentialMinimalParameters{} , v.lj_rcut , lj_e , lj_de );
	v.lj_ecut = lj_e;
      } 

      if( node["rf"] ) v.rf = node["rf"].as<ReactionFieldParms>();
      else v.rf = ReactionFieldParms{};

      return true;
    }
  };
}

