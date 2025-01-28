#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>

#include <onika/physics/units.h>
#include <onika/physics/units.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <onika/physics/constants.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // Coulwolfpair Parameters
  struct CoulWolfParms
  {
    double alpha = 0.0;
    double rc = 0.0;
    double qqrd2e = 14.399645;
    double e_shift = 0.0;
    double f_shift = 0.0;    
  };
  
  // core computation kernel for coul wolf potential
  ONIKA_HOST_DEVICE_FUNC inline void coul_wolf_pair_energy(const CoulWolfParms& p, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    assert( r > 0. );

    double c1 = p_pair.m_atom_a.m_charge;
    double c2 = p_pair.m_atom_b.m_charge;
    
    // LAMMPS
    double prefactor = p.qqrd2e * c1 * c2 / r;
    double erfcc = erfc(p.alpha * r);
    double erfcd = exp(-p.alpha * p.alpha * r * r);
    double v_sh = (erfcc - p.e_shift * r) * prefactor;
    double dvdrr = (erfcc / ( r * r ) + 2.0 * p.alpha / sqrt(M_PI) * erfcd / r) + p.f_shift;
    double forcecoul = dvdrr * r * r * prefactor;
    double fpair = -forcecoul / r;

#   ifdef EXANB_UNITS_V2
    e = EXANB_QUANTITY( v_sh * eV ).convert();
    de = EXANB_QUANTITY( fpair * eV / ang ).convert();
#   else
    e = UnityConverterHelper::convert(v_sh, "eV"); // WARNING : this isn't Cuda compatible
    de = UnityConverterHelper::convert(fpair, "eV/ang");
#   endif

  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::CoulWolfParms>
  {
    static bool decode(const Node& node, exaStamp::CoulWolfParms& v)
    {
      using onika::physics::Quantity;
      if( !node.IsMap() ) { return false; }
      v.alpha = node["alpha"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      v.e_shift = erfc(v.alpha * v.rc) / v.rc;
      v.f_shift = -(v.e_shift + 2.0 * v.alpha / sqrt(M_PI) * exp(-v.alpha * v.alpha * v.rc * v.rc)) / v.rc;      
      return true;
    }
  };
}

