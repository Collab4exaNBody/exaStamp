#include "pair_potential_template.h"
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>

#include <yaml-cpp/yaml.h>
#include <exanb/core/quantity_yaml.h>
#include <exanb/core/basic_types_yaml.h>

#include <exanb/core/quaternion_to_matrix.h>

#ifndef PRIV_NAMESPACE_NAME
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)
#endif

// default definitions doing nothing
#ifndef DEBUG_ADDITIONAL_FIELDS
#define DEBUG_ADDITIONAL_FIELDS /**/
#endif

#ifndef DEBUG_ADDITIONAL_PARAMETERS
#define DEBUG_ADDITIONAL_PARAMETERS /**/
#endif

#ifndef DEBUG_ADDITIONAL_PARAMETER_NAMES
#define DEBUG_ADDITIONAL_PARAMETER_NAMES /**/
#endif

namespace exaStamp
{
  using namespace exanb;

  namespace PRIV_NAMESPACE_NAME
  {
  
    struct PotentialPairParam
    {
      onika::cuda::ro_shallow_copy_t< USTAMP_POTENTIAL_PARAMS > p;
      double rcut = 0.0;
      PairPotentialMinimalParameters pair_params;
      double ecut = 0.0;
    };
  
    struct PotentialMultiParameters
    {
      using UserPotParams = std::map< std::pair<std::string,std::string> , PotentialPairParam >;
      static constexpr size_t MAX_TYPE_PAIR_IDS = 256;
      PotentialPairParam m_pair_params[MAX_TYPE_PAIR_IDS]; // indexed by type pair id
      size_t m_nb_pair_params = 0;
      UserPotParams * m_user_pot_parameters = nullptr;
    };
  
    // helper functor, with templated call operator (with/without weights)
    struct PairMultiForceOp
    {
      const PotentialMultiParameters* p = nullptr;

      // ComputeBuffer less computation with virial
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        Vec3d dr,double d2,
        double& ep,
        double& fx,double& fy,double& fz,
        Mat3d& vir,
        unsigned int type_a, 
        DEBUG_ADDITIONAL_PARAMETERS 
        CellParticlesT* cells,size_t cell_b,size_t p_b,
        double weight ) const
      {
        assert( cells != nullptr );
        const double r = sqrt(d2);
        const unsigned int type_b = cells[cell_b][field::type][p_b];
        const unsigned int pair_id = unique_pair_id(type_a,type_b);
        assert( pair_id < p->m_nb_pair_params );
        const auto& pp = p->m_pair_params[pair_id];
        if( r <= pp.rcut )
        {
          double e=0.0, de=0.0;
          auto tp = pp.pair_params; if( type_a > type_b ) { auto x=tp.m_atom_a; tp.m_atom_a=tp.m_atom_b; tp.m_atom_b=x; }

#         if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
          USTAMP_POTENTIAL_COMPUTE( pp.p, tp, r, e, de , weight );
#         else
          USTAMP_POTENTIAL_COMPUTE( pp.p, tp, r, e, de );
          e *= weight;
          de *= weight;
#         endif
          e -= pp.ecut * weight;
          de /= r;
          
          const Vec3d fe = { de * dr.x , de * dr.y , de * dr.z };
          fx += fe.x;
          fy += fe.y;
          fz += fe.z;
          ep += .5 * e;
          vir += tensor(fe,dr) * -0.5;
        }
      }
      
      // ComputeBuffer less computation with virial
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        Vec3d dr,double d2,
        double& ep,
        double& fx,double& fy,double& fz,
        unsigned int type_a,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT* cells,size_t cell_b,size_t p_b,
        double weight ) const
      {
        assert( cells != nullptr );
        const double r = sqrt(d2);
        const auto type_b = cells[cell_b][field::type][p_b];
        const unsigned int pair_id = unique_pair_id(type_a,type_b);
        const auto& pp = p->m_pair_params[pair_id];      
        if( r <= pp.rcut )
        {
          double e=0.0, de=0.0;
          auto tp = pp.pair_params; if( type_a > type_b ) { auto x=tp.m_atom_a; tp.m_atom_a=tp.m_atom_b; tp.m_atom_b=x; }
          USTAMP_POTENTIAL_COMPUTE( pp.p, tp , r, e, de );
	  e -= pp.ecut;
          de /= r;
          e *= weight;
          de *= weight;        
          const Vec3d fe = { de * dr.x , de * dr.y , de * dr.z };
          fx += fe.x;
          fy += fe.y;
          fz += fe.z;
          ep += .5 * e;
        }
      }

      // ComputeBuffer less computation with virial
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        Vec3d dr,double d2,
        double& fx,double& fy,double& fz,
        Mat3d& vir,
        unsigned int type_a, 
        DEBUG_ADDITIONAL_PARAMETERS 
        CellParticlesT* cells,size_t cell_b,size_t p_b,
        double weight ) const
      {
        assert( cells != nullptr );
        double ep = 0.0;
        const double r = sqrt(d2);
        const unsigned int type_b = cells[cell_b][field::type][p_b];
        const unsigned int pair_id = unique_pair_id(type_a,type_b);
        const auto& pp = p->m_pair_params[pair_id];
        if( r <= pp.rcut )
        {
          double e=0.0, de=0.0;
          auto tp = pp.pair_params; if( type_a > type_b ) { auto x=tp.m_atom_a; tp.m_atom_a=tp.m_atom_b; tp.m_atom_b=x; }
          USTAMP_POTENTIAL_COMPUTE( pp.p, tp , r, e, de );
	  e -= pp.ecut;
          de /= r;
          e *= weight;
          de *= weight;        
          const Vec3d fe = { de * dr.x , de * dr.y , de * dr.z };
          fx += fe.x;
          fy += fe.y;
          fz += fe.z;
          ep += .5 * e;
          vir += tensor(fe,dr) * -0.5;
        }
      }
      
      // ComputeBuffer less computation with virial
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        Vec3d dr,double d2,
        double& fx,double& fy,double& fz,
        unsigned int type_a,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT* cells,size_t cell_b,size_t p_b,
        double weight ) const
      {
        assert( cells != nullptr );
        double ep = 0.0;
        const double r = sqrt(d2);
        const auto type_b = cells[cell_b][field::type][p_b];
        const unsigned int pair_id = unique_pair_id(type_a,type_b);
        const auto& pp = p->m_pair_params[pair_id];      
        if( r <= pp.rcut )
        {
          double e=0.0, de=0.0;
          auto tp = pp.pair_params; if( type_a > type_b ) { auto x=tp.m_atom_a; tp.m_atom_a=tp.m_atom_b; tp.m_atom_b=x; }
          USTAMP_POTENTIAL_COMPUTE( pp.p, tp , r, e, de );
	  e -= pp.ecut;
          de /= r;
          e *= weight;
          de *= weight;        
          const Vec3d fe = { de * dr.x , de * dr.y , de * dr.z };
          fx += fe.x;
          fy += fe.y;
          fz += fe.z;
          ep += .5 * e;
        }
      }

     
    };
  }

}

namespace exanb
{
  template<> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::PairMultiForceOp>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible = true;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif
  };

}

namespace YAML
{
  using exaStamp::PRIV_NAMESPACE_NAME::PotentialMultiParameters;

  template<> struct convert< PotentialMultiParameters >
  {
    static inline bool decode(const Node& node, PotentialMultiParameters& rmpp )
    {
      using namespace exaStamp;
      using namespace exanb;
      using UserPotParams = PotentialMultiParameters::UserPotParams;

      if( ! node.IsSequence() )
      {
        ::exanb::lerr_stream() << "PotentialMultiParameters is not a sequence as expected"<<std::endl;
        return false;
      }
      for(auto pp: node)
      {
        if( ! pp.IsMap() )
        {
          ::exanb::lerr_stream() << "PotentialMultiParameters item is not a map as expected"<<std::endl;
          return false;
        }
        if( ! pp["type_a"] || ! pp["type_b"] || ! pp["rcut"] || ! pp["parameters"] )
        {
          ::exanb::lerr_stream() << "PotentialMultiParameters : missing informations. required information keys are : type_a, type_b, rcut and parameters."<<std::endl;
          return false;
        }
        auto type_a =  pp["type_a"].as<std::string>();
        auto type_b =  pp["type_b"].as<std::string>();
        if( rmpp.m_user_pot_parameters == nullptr ) rmpp.m_user_pot_parameters = new /*( rmpp.m_user_pot_data )*/ UserPotParams {};
        auto & pair_pot = (* rmpp.m_user_pot_parameters) [ std::make_pair(type_a,type_b) ];
        pair_pot.p = pp["parameters"].as<USTAMP_POTENTIAL_PARAMS>();
        pair_pot.rcut = pp["rcut"].as<Quantity>().convert();
        pair_pot.ecut = 0.0;
      }
      return true;
    }
  };
}


