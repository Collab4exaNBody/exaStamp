#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>

#include <onika/physics/units.h>
#include <exanb/core/unityConverterHelper.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exanb/core/physics_constants.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // Coul_cut Parameters
  struct CoulCutParms
  {
    double dielectric = 1.0;
    bool shift = true;
  };

  ONIKA_HOST_DEVICE_FUNC inline void coul_cut_energy(const CoulCutParms& p, const PairPotentialMinimalParameters& p_pair, double r, double& e, double& de)
  {
    assert( r > 0. );
    double qqr2e = 14.399645; // units = eV.ang/e-^2 from LAMMPS

#   ifdef EXANB_UNITS_V2
    double qqrd2e = EXANB_QUANTITY( qqr2e * eV * ang / (ec^2) ).convert() / p.dielectric; 
#   else
    double qqrd2e = UnityConverterHelper::convert(qqr2e, "eV.ang/e-^2") / p.dielectric;
#   endif
    double C = qqrd2e;
    e = C * p_pair.m_atom_a.m_charge * p_pair.m_atom_b.m_charge / r;
    de = -e/r;
  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::CoulCutParms>
  {
    static bool decode(const Node& node, exaStamp::CoulCutParms& v)
    {
      using exanb::Quantity;

      if( !node.IsMap() ) { return false; }
      v.dielectric = node["dielectric"].as<Quantity>().convert();
      v.shift = node["shift"].as<bool>();
      std::cout << "shift val = " << v.shift << std::endl;
      return true;
    }
  };
}

