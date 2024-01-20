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

  struct RigidMoleculePairContext
  {
    Vec3d mol_atom_ra[MAX_RIGID_MOLECULE_ATOMS];
  };

  namespace PRIV_NAMESPACE_NAME
  {
  
    struct RigidMolPotentialPairParam
    {
      onika::cuda::ro_shallow_copy_t< USTAMP_POTENTIAL_PARAMS > p;
      double rcut = 0.0;
      PairPotentialMinimalParameters pair_params;
      double ecut = 0.0;
    };
  
    struct RigidMolPotentialParameters
    {
      using UserPotParams = std::map< std::pair<std::string,std::string> , RigidMolPotentialPairParam >;
      static constexpr size_t MAX_TYPE_PAIR_IDS = 16;
      RigidMolPotentialPairParam m_pair_params[MAX_TYPE_PAIR_IDS]; // indexed by type pair id
      RigidMoleculeAtom m_atoms[MAX_RIGID_MOLECULE_ATOMS];
      unsigned int m_n_atoms = 0;
      unsigned int m_nb_pair_params = 0;
      UserPotParams * m_user_pot_parameters = nullptr;
      //uint8_t m_user_pot_data[ sizeof(UserPotParams) ];
    };
  
    // helper functor, with templated call operator (with/without weights)
    struct RigidMolForceOp
    {
      const RigidMolPotentialParameters* __restrict__ p;

     // without virial computation
      template<bool UseWeights, class CellParticlesT, class Mat3dT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        const Mat3dT& /* _virial */,
        const Quaternion& orient_a,
        Vec3d& _couple,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT cells
        ) const
      {
        const auto & p = *(this->p);
        
        Vec3d mol_atom_ra[MAX_RIGID_MOLECULE_ATOMS];

        // pre-compute atom a force site positions
        if( n > 0 )
        {
          const Mat3d mat_a = quaternion_to_matrix( orient_a );        
          for(unsigned int sub_atom_a=0; sub_atom_a<p.m_n_atoms; ++sub_atom_a)
          {
            mol_atom_ra[sub_atom_a] = mat_a * p.m_atoms[sub_atom_a].m_pos;
          }
        }
        
        //int nPairs = 0;
        for(size_t i=0;i<n;i++)
        {
          size_t cell_b=0, p_b=0;
          tab.nbh.get(i,cell_b,p_b);
          const auto orient_b = cells[cell_b][field::orient][p_b];
          const Mat3d mat_b = quaternion_to_matrix( orient_b );
          const Vec3d mol_b_r = { tab.drx[i] , tab.dry[i] , tab.drz[i] };
          const auto weight = tab.nbh_data.get(i);
/*
          Vec3d mol_atom_rb[MAX_RIGID_MOLECULE_ATOMS];
          for(unsigned int sub_atom_b=0; sub_atom_b<p.m_n_atoms; ++sub_atom_b)
          {
            mol_atom_rb[sub_atom_b] = mat_b * p.m_atoms[sub_atom_b].m_pos;
          }
*/        
          for(unsigned int sub_atom_a=0; sub_atom_a<p.m_n_atoms; ++sub_atom_a)
          {
            const unsigned int type_a = p.m_atoms[sub_atom_a].m_atom_type;
            const Vec3d ra = mol_atom_ra[sub_atom_a]; // mat_a * p.m_atoms[sub_atom_a].m_pos;
            Vec3d site_a_F = { 0., 0., 0. };
            double site_a_ep = 0.;
          
            // calcul sub-a
            for(unsigned int sub_atom_b=0; sub_atom_b<p.m_n_atoms; ++sub_atom_b)
            {
              const unsigned int type_b = p.m_atoms[sub_atom_b].m_atom_type;
              const unsigned int pair_id = unique_pair_id(type_a,type_b);
              const auto& pp = p.m_pair_params[pair_id]; 
              const Vec3d rb =  mol_b_r + /*mol_atom_rb[sub_atom_b]*/ mat_b * p.m_atoms[sub_atom_b].m_pos;
              const Vec3d dr = rb - ra;
              const double rcut2 = pp.rcut * pp.rcut;
              const double d2 = norm2(dr);
              if( d2 <= rcut2 && weight > 0.0 )
              {
                //++ nPairs;
                const double r = sqrt(d2);
                double e=0.0, de=0.0;

#               if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
                USTAMP_POTENTIAL_COMPUTE( pp.p, pp.pair_params, r, e, de , weight );
#               else
                USTAMP_POTENTIAL_COMPUTE( pp.p, pp.pair_params, r, e, de );
                e *= weight;
                de *= weight;
#               endif
                e -= pp.ecut * weight;
                de /= r;
                
                site_a_F += de * dr;
                site_a_ep += e * 0.5;
              }
            }
            // _virial += ???
            _ep += site_a_ep;
            _fx += site_a_F.x;
            _fy += site_a_F.y;
            _fz += site_a_F.z;
            _couple += cross( ra , site_a_F );
          }
        }
        //printf("%d nbh, %d interactions\n",int(n),int(nPairs));
      }

      // without virial computation
      template<bool UseWeights, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& ep,
        double& fx,
        double& fy,
        double& fz,
        const Quaternion& orient_a,
        Vec3d& couple_a,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT cells
        ) const
      {
        //using Mat3d = ::exanb::FakeMat3d;
        this->operator () ( n, tab, ep, fx, fy, fz, FakeMat3d{}, orient_a, couple_a, DEBUG_ADDITIONAL_PARAMETER_NAMES cells );
      }


      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        RigidMoleculePairContext& p_a_ctx,
        CellParticlesT* cells,
        size_t cell_a,
        size_t p_a,
        exanb::ComputePairParticleContextStart
        ) const
      {
        const auto & p = *(this->p);
        const Mat3d mat_a = quaternion_to_matrix( cells[cell_a][field::orient][p_a] );
        for(unsigned int sub_atom_a=0; sub_atom_a<p.m_n_atoms; ++sub_atom_a)
        {
          p_a_ctx.mol_atom_ra[sub_atom_a] = mat_a * p.m_atoms[sub_atom_a].m_pos;
        }
      }

      template<class CellParticlesT, class Mat3dT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        const RigidMoleculePairContext& p_a_ctx,
        const Vec3d& dr, double d2,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        const Mat3dT& /* _virial */,
        const Quaternion& orient_a,
        Vec3d& _couple,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT cells,size_t cell_b,size_t p_b,
        double weight
        ) const
      {
        const auto & p = *(this->p);
        //int nPairs = 0;
        //for(size_t i=0;i<n;i++)
        //{
        //  size_t cell_b=0, p_b=0;
        //  tab.nbh.get(i,cell_b,p_b);
          const auto orient_b = cells[cell_b][field::orient][p_b];
          const Mat3d mat_b = quaternion_to_matrix( orient_b );
          const Vec3d mol_b_r = dr; // { tab.drx[i] , tab.dry[i] , tab.drz[i] };
        //  const auto weight = tab.nbh_data.get(i);
        
          for(unsigned int sub_atom_a=0; sub_atom_a<p.m_n_atoms; ++sub_atom_a)
          {
            const unsigned int type_a = p.m_atoms[sub_atom_a].m_atom_type;
            const Vec3d ra = p_a_ctx.mol_atom_ra[sub_atom_a]; // mat_a * p.m_atoms[sub_atom_a].m_pos;
            Vec3d site_a_F = { 0., 0., 0. };
            double site_a_ep = 0.;
          
            // calcul sub-a
            for(unsigned int sub_atom_b=0; sub_atom_b<p.m_n_atoms; ++sub_atom_b)
            {
              const unsigned int type_b = p.m_atoms[sub_atom_b].m_atom_type;
              const unsigned int pair_id = unique_pair_id(type_a,type_b);
              const auto& pp = p.m_pair_params[pair_id]; 
              const Vec3d rb =  mol_b_r + mat_b * p.m_atoms[sub_atom_b].m_pos;
              const Vec3d dr = rb - ra;
              const double rcut2 = pp.rcut * pp.rcut;
              const double d2 = norm2(dr);
              if( d2 <= rcut2 && weight > 0.0 )
              {
                //++ nPairs;
                const double r = sqrt(d2);
                double e=0.0, de=0.0;
                
#               if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
                USTAMP_POTENTIAL_COMPUTE( pp.p, pp.pair_params, r, e, de , weight );
#               else
                USTAMP_POTENTIAL_COMPUTE( pp.p, pp.pair_params, r, e, de );
                e *= weight;
                de *= weight;
#               endif
                e -= pp.ecut * weight;
                de /= r;
                
                site_a_F += de * dr;
                site_a_ep += e * 0.5;
              }
            }
            // _virial += ???
            _ep += site_a_ep;
            _fx += site_a_F.x;
            _fy += site_a_F.y;
            _fz += site_a_F.z;
            _couple += cross( ra , site_a_F );
          }
        //}
        //printf("%d nbh, %d interactions\n",int(n),int(nPairs));
      }

      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        const RigidMoleculePairContext& p_a_ctx,
        const Vec3d& dr, double d2,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        //Mat3dT& /* _virial */,
        const Quaternion& orient_a,
        Vec3d& _couple,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT cells,size_t cell_b,size_t p_b,
        double weight
        ) const
      {
        this->operator() ( p_a_ctx, dr, d2, _ep, _fx, _fy, _fz, FakeMat3d{}, orient_a, _couple, DEBUG_ADDITIONAL_PARAMETER_NAMES cells, cell_b, p_b, weight );
      }

    };
  }

}

namespace exanb
{
  template<> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::RigidMolForceOp>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = false;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif
  };

/*
  template<> struct ComputePairParticleContextTraits<exaStamp::PRIV_NAMESPACE_NAME::RigidMolForceOp>
  {
    using ParticleContext = exaStamp::RigidMoleculePairContext;    
    static inline constexpr bool HasParticleContextStart = true;
    static inline constexpr bool HasParticleContextStop = false;    
  };
*/

}

namespace YAML
{
  using exaStamp::PRIV_NAMESPACE_NAME::RigidMolPotentialParameters;

  template<> struct convert< RigidMolPotentialParameters >
  {
    static inline bool decode(const Node& node, RigidMolPotentialParameters& rmpp )
    {
      using namespace exaStamp;
      using namespace exanb;
      using UserPotParams = RigidMolPotentialParameters::UserPotParams;
      if( ! node.IsSequence() )
      {
        ::exanb::lerr_stream() << "RigidMolPotentialParameters is not a sequence as expected"<<std::endl;
        return false;
      }
      for(auto pp: node)
      {
        if( ! pp.IsMap() )
        {
          ::exanb::lerr_stream() << "RigidMolPotentialParameters item is not a map as expected"<<std::endl;
          return false;
        }
        if( ! pp["type_a"] || ! pp["type_b"] || ! pp["rcut"] || ! pp["parameters"] )
        {
          ::exanb::lerr_stream() << "RigidMolPotentialParameters : missing informations. required information keys are : type_a, type_b, rcut and parameters."<<std::endl;
          return false;
        }
        auto type_a =  pp["type_a"].as<std::string>();
        auto type_b =  pp["type_b"].as<std::string>();
        if( rmpp.m_user_pot_parameters == nullptr ) rmpp.m_user_pot_parameters = new /*( rmpp.m_user_pot_data )*/ UserPotParams {};
        auto & pair_pot = (* rmpp.m_user_pot_parameters) [ std::make_pair(type_a,type_b) ];
        pair_pot.p = pp["parameters"].as<USTAMP_POTENTIAL_PARAMS>();
        pair_pot.rcut = pp["rcut"].as<Quantity>().convert();
      }
      return true;
    }
  };
}


