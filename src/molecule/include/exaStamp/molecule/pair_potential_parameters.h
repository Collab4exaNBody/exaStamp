#pragma once

#include <yaml-cpp/yaml.h>
#include <vector>
#include <string>

#include <exanb/core/quantity_yaml.h>
#include <exanb/core/log.h>

#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>
#include <exaStamp/potential/pair_potentials/exp6/exp6.h>
#include <exaStamp/potential/reaction_field/reaction_field.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  struct IntramolecularRFParam
  {
    ReactionFieldParms m_param = {};
    double m_c1 = 0.0;
    double m_c2 = 0.0;
    inline bool operator < ( const IntramolecularRFParam& r ) const
    {
      // ecut and ecut_rf omitted on prupose
      const std::array<double,7> A = { m_param.rc, m_param.RF0, m_param.RF1, m_param.RF2, m_param.ecut, m_c1, m_c2 };
      const std::array<double,7> B = { r.m_param.rc, r.m_param.RF0, r.m_param.RF1, r.m_param.RF2, r.m_param.ecut, r.m_c1, r.m_c2 };
      return A < B;
    }

  };

  struct IntramolecularLJExp6Param
  {
    double m_A = 0.0;         // only meaningful for Exp6
    double m_B_ISEXP6 = 1.0;  // if == 0.0, then it's a LJ potential
    double m_C_EPSILON = 0.0; // either Exp6's parameter C or LJ's parameter Epsilon
    double m_D_SIGMA = 0.0;   // either Exp6's parameter D or LJ's parameter Sigma
    double m_rcut = 0.0;
    double m_ecut = 0.0;
    
    ONIKA_HOST_DEVICE_FUNC inline bool is_lj() const
    {
      return m_B_ISEXP6==0.0;
    }
    ONIKA_HOST_DEVICE_FUNC inline bool is_exp6() const
    {
      return m_B_ISEXP6!=0.0;
    }
    
    inline bool operator < ( const IntramolecularLJExp6Param& r ) const
    {
      // ecut and ecut_rf omitted on prupose
      const std::array<double,6> A = { m_A, m_B_ISEXP6, m_C_EPSILON, m_D_SIGMA, m_rcut, m_ecut };
      const std::array<double,6> B = { r.m_A, r.m_B_ISEXP6, r.m_C_EPSILON, r.m_D_SIGMA, r.m_rcut, r.m_ecut };
      return A < B;
    }
    
    ONIKA_HOST_DEVICE_FUNC inline void compute_energy(double r,double &e, double& de) const
    {
      if( is_lj() )
      {
        lj_compute_energy( LennardJonesParms{m_C_EPSILON,m_D_SIGMA}, PairPotentialMinimalParameters{}, r, e, de );
      }
      else
      {
        assert( is_exp6() );
        exp6_compute_energy( Exp6Parms{m_A,m_B_ISEXP6,m_C_EPSILON,m_D_SIGMA}, PairPotentialMinimalParameters{}, r, e, de );
      }
      e -= m_ecut;
    }
    ONIKA_HOST_DEVICE_FUNC inline void update_ecut()
    {
      double e=0.0, de=0.0;
      m_ecut = 0.0;
      compute_energy(m_rcut,e,de);
      m_ecut = e;
    }
  };
  
  struct IntramolecularPairUserParam
  {
    std::string m_type_a;
    std::string m_type_b;
    IntramolecularRFParam m_rf = {};
    IntramolecularLJExp6Param m_ljexp6 = {};
  };

  using IntramolecularPairUserParamVector = std::vector<IntramolecularPairUserParam>;

}

// Yaml conversion operators, allows to read bonds potentials parameters from config file
namespace YAML
{
  using exanb::fatal_error;
  using exaStamp::IntramolecularPairUserParam;
  using exaStamp::ReactionFieldParms;
  using exaStamp::Exp6Parms;
  using exaStamp::LennardJonesParms;
  using exaStamp::IntramolecularRFParam;
  using exaStamp::IntramolecularLJExp6Param;

  template<> struct convert<IntramolecularPairUserParam>
  {
    static bool decode(const Node& node, IntramolecularPairUserParam& p)
    {
      if(!node["type_a"])
      {
        fatal_error() << "atom type a missing" << std::endl;
      }

      if(!node["type_b"])
      {
        fatal_error() << "atom type b missing" << std::endl;
      }

      p.m_type_a = node["type_a"].as<std::string>();
      p.m_type_b = node["type_b"].as<std::string>();
      p.m_rf = IntramolecularRFParam{};
      p.m_ljexp6 = IntramolecularLJExp6Param{};

      if( node["parameters"] )
      {
        if( node["parameters"]["lj"] )
        {
          if( ! node["parameters"]["lj_rcut"] ) { fatal_error()<<"missing rcut for LJ"<<std::endl; }
          auto lj = node["parameters"]["lj"].as<LennardJonesParms>();
          p.m_ljexp6.m_B_ISEXP6 = 0.0;
          p.m_ljexp6.m_C_EPSILON = lj.epsilon;
          p.m_ljexp6.m_D_SIGMA = lj.sigma;
          assert( p.m_ljexp6.is_lj() );
          p.m_ljexp6.m_rcut = node["parameters"]["lj_rcut"].as<Quantity>().convert();
        }
        else if( node["parameters"]["exp6"] )
        {
          if( ! node["parameters"]["exp6_rcut"] ) { fatal_error()<<"missing rcut for Exp6"<<std::endl; }
          auto exp6 = node["parameters"]["exp6"].as<Exp6Parms>();
          p.m_ljexp6.m_A = exp6.A;
          p.m_ljexp6.m_B_ISEXP6 = exp6.B;
          p.m_ljexp6.m_C_EPSILON = exp6.C;
          p.m_ljexp6.m_D_SIGMA = exp6.D;          
          assert( p.m_ljexp6.is_exp6() );
          p.m_ljexp6.m_rcut = node["parameters"]["exp6_rcut"].as<Quantity>().convert();
        }
        if( node["parameters"]["rf"] )
        {
          p.m_rf.m_param = node["parameters"]["rf"].as<ReactionFieldParms>();
        }
      }
      return true;
    }
  };
}
