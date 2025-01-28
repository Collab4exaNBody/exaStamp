#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>

#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>
#include <exaStamp/potential/pair_potentials/coul_wolf_pair/coul_wolf_pair.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // assembled Lennard-Jones & Reaction Field parameters
  struct LJWOLFParms
  {
    double lj_rcut = 0.0;
    double lj_ecut = 0.0;
    LennardJonesParms lj;
    CoulWolfParms cw;
    bool shiftlj;
  };

//# pragma omp declare simd uniform(p,ppp) notinbranch
  ONIKA_HOST_DEVICE_FUNC inline void ljwolf_energy(const LJWOLFParms& p, const PairPotentialMinimalParameters& pp, double r, double& e, double& de)
  {
    double lj_e=0.0, lj_de=0.0;
    if( r <= p.lj_rcut )
    {
      lj_compute_energy( p.lj , pp , r , lj_e , lj_de );
      lj_e -= p.lj_ecut;
      //printf("LJ: r=%g, rcut=%g, e=%g, de=%g\n",r,p.lj_rcut,lj_e,lj_de);
    }

    double cw_e=0.0, cw_de=0.0;
    if( r <= p.cw.rc )
    {
      coul_wolf_pair_energy( p.cw, pp, r, cw_e, cw_de );
      //printf("RF: r=%g, rcut=%g, e=%g, de=%g\n",r,p.rf.rc,rf_e,rf_de);
    }
    
    e = lj_e + cw_e;
    de = lj_de + cw_de;

    //printf("r=% .6e : e=% .6e , de=% .6e : c1=% .6e , c2=% .6e\n",r,e,de, pp.m_atom_a.m_charge, pp.m_atom_b.m_charge);
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::LJWOLFParms;
  using exaStamp::LennardJonesParms;
  using exaStamp::CoulWolfParms;
  using exanb::Quantity;
  using exaStamp::lj_compute_energy;
  using exaStamp::PairPotentialMinimalParameters;

  template<> struct convert<LJWOLFParms>
  {
    static bool decode(const Node& node, LJWOLFParms& v)
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

      if( node["cw"] ) v.cw = node["cw"].as<CoulWolfParms>();
      else v.cw = CoulWolfParms{};

      return true;
    }
  };
}

