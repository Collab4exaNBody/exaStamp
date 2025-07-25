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

#include <exaStamp/potential/pair_potentials/exp6/exp6.h>
#include <exaStamp/potential/reaction_field/reaction_field.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // assembled Lennard-Jones & Reaction Field parameters
  struct Exp6RFParms
  {
    double exp6_rcut = 0.0;
    double exp6_ecut = 0.0;
    Exp6Parms exp6;
    // double rf_rcut;
    ReactionFieldParms rf;
  };

//# pragma omp declare simd uniform(p,ppp) notinbranch
  ONIKA_HOST_DEVICE_FUNC inline void exp6rf_energy(const Exp6RFParms& p, const PairPotentialMinimalParameters& pp, double r, double& e, double& de)
  {
    double exp6_e=0.0, exp6_de=0.0;
    if( r <= p.exp6_rcut )
    {
      exp6_compute_energy( p.exp6 , pp , r , exp6_e , exp6_de );
      exp6_e -= p.exp6_ecut;
    }

    double rf_e=0.0, rf_de=0.0;
    if( r <= p.rf.rc /*_rcut*/ )
    {
      reaction_field_compute_energy( p.rf, pp.m_atom_a.m_charge * pp.m_atom_b.m_charge, r, rf_e, rf_de );
    }

    e = exp6_e + rf_e;
    de = exp6_de + rf_de;
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::Exp6RFParms;
  using exaStamp::Exp6Parms;
  using exaStamp::ReactionFieldParms;
  using onika::physics::Quantity;
  using exaStamp::exp6_compute_energy;
  using exaStamp::PairPotentialMinimalParameters;

  template<> struct convert<Exp6RFParms>
  {
    static bool decode(const Node& node, Exp6RFParms& v)
    {
      if( !node.IsMap() ) { return false; }
      
      if( node["exp6_rcut"] ) v.exp6_rcut = node["exp6_rcut"].as<Quantity>().convert();
      else v.exp6_rcut = 0.0;
      
      if( node["exp6"] ) v.exp6 = node["exp6"].as<Exp6Parms>();
      else v.exp6 = Exp6Parms{ 0. , 1. , 0. , 0. };
      v.exp6_ecut = 0.0;
      
      double exp6_e=0.0, exp6_de=0.0;
      exp6_compute_energy( v.exp6 , PairPotentialMinimalParameters{} , v.exp6_rcut , exp6_e , exp6_de );
      v.exp6_ecut = exp6_e;

      if( node["rf"] ) v.rf = node["rf"].as<ReactionFieldParms>();
      else v.rf = ReactionFieldParms{ 0., 0. , 0. , 0. , 0. };

      return true;
    }
  };
}

