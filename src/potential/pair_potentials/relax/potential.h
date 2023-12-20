#pragma once

#include <cmath>
#include <limits>
#include <utility>
#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exanb/core/physics_constants.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // Relax Parameters
  struct RelaxParms
  {
    double r1 = 0.1;
    double rc = 2.0;
  };
  
  ONIKA_HOST_DEVICE_FUNC inline void relax_compute_energy(const RelaxParms& p, const PairPotentialMinimalParameters& pair, double rij, double& e, double& de)
  {
    double r = rij;
    if( r < p.r1 ) r = p.r1;
    if( r > p.rc ) r = p.rc;
    e = ( p.rc / r ) - 1.0 ;
    de = - e;
  }

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<exaStamp::RelaxParms>
  {
    static bool decode(const Node& node, exaStamp::RelaxParms& v)
    {
      using exanb::Quantity;
      v = exaStamp::RelaxParms{};
      if( !node.IsMap() ) { return false; }
      if( ! node["r1"] ) { return false; }
      if( ! node["rc"] ) { return false; }
      v.r1 = node["r1"].as<Quantity>().convert();
      v.rc = node["rc"].as<Quantity>().convert();
      return true;
    }
  };
}

#define USTAMP_POTENTIAL_NAME     relax
#define USTAMP_POTENTIAL_PARAMS   RelaxParms
#define USTAMP_POTENTIAL_COMPUTE  relax_compute_energy
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) char(0)
#define USTAMP_POTENTIAL_ENABLE_CUDA 1

