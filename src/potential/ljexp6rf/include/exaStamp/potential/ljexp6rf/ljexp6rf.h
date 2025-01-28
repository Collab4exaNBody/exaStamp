#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>

#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential/pair_potentials/exp6/exp6.h>
#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>
#include <exaStamp/potential/reaction_field/reaction_field.h>
#include <exaStamp/compute/force_energy.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  // assembled Lennard-Jones & Reaction Field parameters
  struct LJExp6RFParms
  {
    double m_A = 0.0;         // only meaningful for Exp6
    double m_B_ISEXP6 = 1.0;  // if == 0.0, then it's a LJ potential
    double m_C_EPSILON = 0.0; // either Exp6's parameter C or LJ's parameter Epsilon
    double m_D_SIGMA = 0.0;   // either Exp6's parameter D or LJ's parameter Sigma
    double m_rcut = 0.0;
    double m_ecut = 0.0;

    ReactionFieldParms m_rf = {};
    
    ONIKA_HOST_DEVICE_FUNC inline bool pair_is_null() const
    {
      return m_A==0.0 && m_C_EPSILON==0.0 && m_D_SIGMA==0.0 && m_ecut==0.0;
    }

    ONIKA_HOST_DEVICE_FUNC inline bool is_null() const
    {
      return pair_is_null() && m_rf.is_null();
    }

    ONIKA_HOST_DEVICE_FUNC inline bool is_lj() const
    {
      return m_B_ISEXP6==0.0;
    }

    ONIKA_HOST_DEVICE_FUNC inline bool is_exp6() const
    {
      return m_B_ISEXP6!=0.0;
    }

    inline void set_lj_parameters(double sigma, double epsilon, double rcut)
    {
      m_B_ISEXP6 = 0.0;
      m_C_EPSILON = epsilon;
      m_D_SIGMA = sigma;
      m_rcut = rcut;
      m_ecut = 0.0;
      assert( is_lj() );
    }
    
    inline void set_exp6_parameters(double A, double B, double C, double D, double rcut)
    {
      m_A = A;
      m_B_ISEXP6 = B;
      m_C_EPSILON = C;
      m_D_SIGMA = D;
      m_rcut = rcut;
      m_ecut = 0.0;
      assert( is_exp6() );
    }

    inline void set_reaction_field_parameters(double rc, double epsilon)
    {
      init_rf( m_rf , rc , epsilon );
    }

    inline LJExp6RFParms scale(double pair_weight, double rf_weight) const
    {
      // LJ linearly scales with epsilon
      // Exp6 linearly scales with A,C and D terms (not B)
      // Reaction Field linearly scales with RF0, RF1, RF2 and ecut (ecut optional, can be recomputed from RF0, RF1 and RF2 as well)
      LJExp6RFParms p = *this;
      if( p.is_lj() )
      {
        p.m_C_EPSILON *= pair_weight;
      }
      else
      {
        p.m_A         *= pair_weight;
        p.m_C_EPSILON *= pair_weight;
        p.m_D_SIGMA   *= pair_weight;
      }
      p.m_ecut        *= pair_weight;
      
      p.m_rf.RF0      *= rf_weight;
      p.m_rf.RF1      *= rf_weight;
      p.m_rf.RF2      *= rf_weight;
      p.m_rf.ecut     *= rf_weight;
      
      return p;
    }

    inline bool operator < ( const LJExp6RFParms& r ) const
    {
      // ecut and ecut_rf omitted on prupose
      const std::array<double,11> A = {   m_A,   m_B_ISEXP6,   m_C_EPSILON,   m_D_SIGMA,   m_rcut,   m_ecut,   m_rf.rc,   m_rf.RF0,   m_rf.RF1,   m_rf.RF2,   m_rf.ecut };
      const std::array<double,11> B = { r.m_A, r.m_B_ISEXP6, r.m_C_EPSILON, r.m_D_SIGMA, r.m_rcut, r.m_ecut, r.m_rf.rc, r.m_rf.RF0, r.m_rf.RF1, r.m_rf.RF2, r.m_rf.ecut };
      return A < B;
    }
    inline bool operator == ( const LJExp6RFParms& r ) const
    {
      return ! ( *this < r )
          && ! ( r < *this );
    }

    ONIKA_HOST_DEVICE_FUNC ForceEnergy compute_force_energy(double r, double charge_product , double pair_weight=1.0, double rf_weight=1.0 ) const
    {
      // pair potential part
      double pair_e = 0.0;
      double pair_de = 0.0;
      if( is_lj() )
      {
        lj_compute_energy( LennardJonesParms{m_C_EPSILON,m_D_SIGMA}, PairPotentialMinimalParameters{}, r, pair_e, pair_de );
      }
      else
      {
        assert( is_exp6() );
        exp6_compute_energy( Exp6Parms{m_A,m_B_ISEXP6,m_C_EPSILON,m_D_SIGMA}, PairPotentialMinimalParameters{}, r, pair_e, pair_de );
      }
      // energy shift at cutoff for pair potential
      pair_e -= m_ecut;
      pair_e *= pair_weight;
      pair_de *= pair_weight;

      // reaction field part
      double rf_e = 0.0;
      double rf_de = 0.0;
      reaction_field_compute_energy( m_rf,  charge_product, r, rf_e, rf_de );
      rf_e *= rf_weight;
      rf_de *= rf_weight;
      //printf("r=%.5e pair_de=%.5e rf_de=%.5e\n",r,pair_de,rf_de);

      return { pair_de + rf_de , pair_e + rf_e };
    }
    
    ONIKA_HOST_DEVICE_FUNC inline void update_ecut()
    {
      m_ecut = 0.0;
      const auto [ de , e ] = compute_force_energy(m_rcut,0.0);
      m_ecut = e;
    }

    template<class StreamT>
    inline StreamT& to_stream(StreamT& out) const
    {
      const char* sep = "";
      out << std::setprecision(4);
      if( m_rcut>0.0 ) 
      {
        if(is_lj()) out<<"LJ(epsilon="<<m_C_EPSILON<<",sigma="<<m_D_SIGMA;
        else out<<"Exp6(A="<<m_A<<",B="<<m_B_ISEXP6<<",C="<<m_C_EPSILON<<",D="<<m_D_SIGMA;
        out<<",rcut="<<m_rcut<<",ecut="<<m_ecut<<")";
        sep = "+";
      }
      if(m_rf.rc>0.0) out << sep << "RF(rc="<<m_rf.rc<<",RF0="<<m_rf.RF0<<",RF1="<<m_rf.RF1<<",RF2="<<m_rf.RF2<<",ecut="<<m_rf.ecut<<")";
      return out;
    }

  };

  inline std::ostream& operator << (std::ostream& out , const LJExp6RFParms& p) { return p.to_stream(out); }

  struct LJExp6RFMultiParmsPair
  {
    std::string m_type_a;
    std::string m_type_b;
    LJExp6RFParms m_params;
  };
  
  struct LJExp6RFMultiParms
  {
    std::vector<LJExp6RFMultiParmsPair> m_potentials;
    inline LJExp6RFParms* params_for_pair(const std::string& type_a, const std::string& type_b)
    {
      for(auto& elem : m_potentials)
      {
        if( ( elem.m_type_a == type_a && elem.m_type_b == type_b ) || ( elem.m_type_b == type_a && elem.m_type_a == type_b ) ) return & elem.m_params;
      }
      return nullptr;
    }
  };

}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  
  template<> struct convert<exaStamp::LJExp6RFParms>
  {
    static bool decode(const Node& node, exaStamp::LJExp6RFParms& v)
    {
      if( !node.IsMap() ) { return false; }
      if( node["exp6"] && node["lj"] )
      {
        exanb::fatal_error()<<"Specifying both Exp6 and LJ parameters is forbidden. Choose one or another."<<std::endl;
      }
      
      v = exaStamp::LJExp6RFParms{};
      
      double rcut = 0.0;
      if( node["rcut"] )
      {
        rcut = node["rcut"].as<exanb::Quantity>().convert();
      }

      if( node["exp6"] )
      {
        auto p = node["exp6"].as<exaStamp::Exp6Parms>();
        v.set_exp6_parameters( p.A, p.B, p.C, p.D , rcut );
      }
      else if( node["lj"] )
      {
        auto p = node["lj"].as<exaStamp::LennardJonesParms>();
        v.set_lj_parameters( p.sigma, p.epsilon , rcut );
      }

      if( node["rf"] )
      {
        v.m_rf = node["rf"].as<exaStamp::ReactionFieldParms>();
      }
      
      v.update_ecut();

      return true;
    }
  };
  
  template<> struct convert<exaStamp::LJExp6RFMultiParms>
  {
    static bool decode(const Node& node, exaStamp::LJExp6RFMultiParms& v)
    {
      v.m_potentials.clear();
      if( !node.IsSequence() )
      {
        exanb::fatal_error() << "LJExp6RFMultiParms is not a list as expected" << std::endl;
        return false;
      }
      for(auto pairmat : node)
      {
        if( ! pairmat.IsMap() )
        {
          exanb::fatal_error() << "LJExp6RFMultiParms pair entry is not a map as expected" << std::endl;
          return false;
        }
        if( !pairmat["type_a"] ) { exanb::fatal_error()<<"missing keyword type_a in material pair description"<<std::endl; }
        if( !pairmat["type_b"] ) { exanb::fatal_error()<<"missing keyword type_b in material pair description"<<std::endl; }
        if( !pairmat["potential"] ) { exanb::fatal_error()<<"missing keyword potential in material pair description"<<std::endl; }
        v.m_potentials.push_back( exaStamp::LJExp6RFMultiParmsPair{ pairmat["type_a"].as<std::string>() , pairmat["type_b"].as<std::string>() , pairmat["potential"].as<exaStamp::LJExp6RFParms>() } );
      }
      return true;
    }
  };

  
}

