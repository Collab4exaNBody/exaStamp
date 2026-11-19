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
#include <limits>
#include <yaml-cpp/yaml.h>
#include <exaStamp/yaml_check_keys.h>
#include <onika/physics/units.h>

#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/potential/pair_potentials/exp6/exp6.h>
#include <exaStamp/potential/pair_potentials/lennard_jones/lennard_jones.h>
#include <exaStamp/potential/pair_potentials/coul_rf/coul_rf.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  struct ForceEnergy
  {
    double m_force = 0.0;
    double m_energy = 0.0;
  };
  
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

    // true when the pair has a LJ or Exp6 term (an rf-only pair has none)
    ONIKA_HOST_DEVICE_FUNC inline bool has_pair_term() const
    {
      return m_A!=0.0 || m_C_EPSILON!=0.0 || m_D_SIGMA!=0.0;
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
      //      p.m_ecut        *= pair_weight;
      
      p.m_rf.RF0      *= rf_weight;
      p.m_rf.RF1      *= rf_weight;
      p.m_rf.RF2      *= rf_weight;
      //      p.m_rf.ecut     *= rf_weight;
      
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

    // unshifted LJ or Exp6 term, regardless of cutoff
    ONIKA_HOST_DEVICE_FUNC inline void compute_pair_term(double r, double& pair_e, double& pair_de) const
    {
      if( is_lj() )
      {
        lj_compute_energy( LennardJonesParms{m_C_EPSILON,m_D_SIGMA}, PairPotentialMinimalParameters{}, r, pair_e, pair_de );
      }
      else
      {
        assert( is_exp6() );
        exp6_compute_energy( Exp6Parms{m_A,m_B_ISEXP6,m_C_EPSILON,m_D_SIGMA}, PairPotentialMinimalParameters{}, r, pair_e, pair_de );
      }
    }

    // Each term only contributes below its own cutoff (m_rcut for LJ/Exp6, m_rf.rc for reaction field), since the
    // neighbor distance cutoff is the largest of all pair and reaction field cutoffs. A zero cutoff means that the
    // caller applies the cutoff itself.
    ONIKA_HOST_DEVICE_FUNC ForceEnergy compute_force_energy(double r, double charge_product , double pair_weight=1.0, double rf_weight=1.0 ) const
    {
      // pair potential part
      double pair_e = 0.0;
      double pair_de = 0.0;
      if( has_pair_term() && ( m_rcut <= 0.0 || r < m_rcut ) )
      {
        compute_pair_term( r , pair_e , pair_de );
        // energy shift at cutoff for pair potential
        pair_e -= m_ecut;
      }
      pair_e *= pair_weight;
      pair_de *= pair_weight;

      // reaction field part
      double rf_e = 0.0;
      double rf_de = 0.0;
      if( m_rf.rc <= 0.0 || r < m_rf.rc )
      {
        reaction_field_compute_energy_weighted( m_rf,  charge_product, r, rf_e, rf_de , rf_weight);
      }

      return { pair_de + rf_de , pair_e + rf_e };
    }

    ONIKA_HOST_DEVICE_FUNC double compute_self_energy(double charge_product) const
    {
      double rf_self_e;
      reaction_field_compute_self_energy( m_rf,  charge_product, rf_self_e);
      return rf_self_e;
    }    
    
    // energy shift of the LJ/Exp6 term at its own cutoff. Without an rcut inside the parameters (multi-species
    // and multimat inputs give it on the pair entry) there is no shift here : the pair templates shift the energy
    // at the pair cutoff themselves. Never evaluated at r = 0 (would give NaN).
    ONIKA_HOST_DEVICE_FUNC inline void update_ecut()
    {
      m_ecut = 0.0;
      if( m_rcut > 0.0 && has_pair_term() )
      {
        double de = 0.0;
        compute_pair_term( m_rcut , m_ecut , de );
      }
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

  // Single precision copy of LJExp6RFParms, for the mixed precision kernel (ljexp6rf_pc_fp32): the pair terms are
  // evaluated in FP32, the caller accumulates forces, energies and virials in FP64.
  // Cutoffs are stored squared (FLT_MAX when the term has no cutoff of its own) and 12/B is precomputed, so that the
  // only division is one reciprocal per pair.
  struct LJExp6RFParmsFP32
  {
    float m_A = 0.0f;
    float m_B = 0.0f;           // Exp6 B, 0 for LJ
    float m_12_over_B = 0.0f;   // Exp6 only
    float m_C_EPSILON = 0.0f;
    float m_D_SIGMA = 0.0f;     // Exp6 D, or LJ sigma
    float m_sigma2 = 0.0f;      // LJ only
    float m_rcut2 = 0.0f;
    float m_ecut = 0.0f;
    float m_rf_rc2 = 0.0f;
    float m_RF0 = 0.0f;
    float m_RF1 = 0.0f;
    float m_RF2 = 0.0f;
    double m_RF2_fp64 = 0.0;    // self energy, constant per particle : kept exact
    bool m_is_lj = false;
    bool m_has_pair_term = false;

    LJExp6RFParmsFP32() = default;

    inline LJExp6RFParmsFP32(const LJExp6RFParms& p)
    {
      constexpr float no_cutoff = std::numeric_limits<float>::max();
      m_is_lj = p.is_lj();
      m_has_pair_term = p.has_pair_term();
      m_A = p.m_A;
      m_B = m_is_lj ? 0.0f : static_cast<float>( p.m_B_ISEXP6 );
      m_12_over_B = m_is_lj ? 0.0f : static_cast<float>( 12.0 / p.m_B_ISEXP6 );
      m_C_EPSILON = p.m_C_EPSILON;
      m_D_SIGMA = p.m_D_SIGMA;
      m_sigma2 = p.m_D_SIGMA * p.m_D_SIGMA;
      m_rcut2 = p.m_rcut > 0.0 ? static_cast<float>( p.m_rcut * p.m_rcut ) : no_cutoff;
      m_ecut = p.m_ecut;
      m_rf_rc2 = p.m_rf.rc > 0.0 ? static_cast<float>( p.m_rf.rc * p.m_rf.rc ) : no_cutoff;
      m_RF0 = p.m_rf.RF0;
      m_RF1 = p.m_rf.RF1;
      m_RF2 = p.m_rf.RF2;
      m_RF2_fp64 = p.m_rf.RF2;
    }

    ONIKA_HOST_DEVICE_FUNC inline double compute_self_energy(double charge_product) const
    {
      return - charge_product * m_RF2_fp64;
    }

    // returns { (dE/dr)/r , E } at squared distance r2, same terms and cutoffs as LJExp6RFParms::compute_force_energy
    ONIKA_HOST_DEVICE_FUNC inline void compute_force_energy(float r2, float charge_product, float& de_over_r, float& e) const
    {
      const float inv_r2 = 1.0f / r2;
      const float r = sqrtf( r2 );
      const float inv_r = r * inv_r2;
      de_over_r = 0.0f;
      e = 0.0f;
      if( m_has_pair_term && r2 < m_rcut2 )
      {
        if( m_is_lj )
        {
          const float s2 = m_sigma2 * inv_r2;
          const float s6 = s2 * s2 * s2;
          const float s12 = s6 * s6;
          e = 4.0f * m_C_EPSILON * ( s12 - s6 );
          de_over_r = -24.0f * m_C_EPSILON * ( 2.0f * s12 - s6 ) * inv_r2;
        }
        else
        {
          const float inv_r6 = inv_r2 * inv_r2 * inv_r2;
          const float Cr6 = m_C_EPSILON * inv_r6;
          const float t = m_12_over_B * inv_r;
          const float t2 = t * t;
          const float t4 = t2 * t2;
          const float Dt12 = m_D_SIGMA * ( t4 * t4 * t4 );
          const float Aexp = m_A * expf( - m_B * r );
          e = Aexp - Cr6 + Dt12;
          de_over_r = ( - m_B * Aexp ) * inv_r + ( 6.0f * Cr6 - 12.0f * Dt12 ) * inv_r2;
        }
        e -= m_ecut;
      }
      if( r2 < m_rf_rc2 )
      {
        e += charge_product * ( m_RF0 * inv_r + m_RF1 * r2 - m_RF2 );
        de_over_r += charge_product * ( 2.0f * m_RF1 - m_RF0 * inv_r2 * inv_r );
      }
    }
  };

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
      exaStamp::check_yaml_keys( node , { "exp6", "lj", "rcut", "rf" } , "ljexp6rf parameters" );
      if( node["exp6"] && node["lj"] )
      {
        exanb::fatal_error()<<"Specifying both Exp6 and LJ parameters is forbidden. Choose one or another."<<std::endl;
      }
      
      v = exaStamp::LJExp6RFParms{};
      
      double rcut = 0.0;
      if( node["rcut"] )
      {
        rcut = node["rcut"].as<onika::physics::Quantity>().convert();
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

