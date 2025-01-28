#include "pair_potential_template.h"
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_operators.h>
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

#if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
# define USTAMP_POTENTIAL_COMPUTE_WEIGHT(PotParam,PairParam,r,e,de,weight) USTAMP_POTENTIAL_COMPUTE(PotParam,PairParam,r,e,de,weight)
#else
# define USTAMP_POTENTIAL_COMPUTE_WEIGHT(PotParam,PairParam,r,e,de,weight) USTAMP_POTENTIAL_COMPUTE(PotParam,PairParam,r,e,de); e*=weight; de*=weight
#endif

namespace exaStamp
{
  using namespace exanb;

  struct RigidMoleculePairContext
  {
    Vec3d mol_atom_ra[MAX_RIGID_MOLECULE_ATOMS];
    Vec3d force;
    Vec3d couple;
    double ep;
    uint16_t type_a;
    uint16_t n_atoms_a;
    uint16_t atoms_start_a;
  };

  namespace PRIV_NAMESPACE_NAME
  {
  
    struct RigidMolPotentialPairParam
    {
      onika::cuda::ro_shallow_copy_t< USTAMP_POTENTIAL_PARAMS > p;
      PairPotentialMinimalParameters pair_params;
      double rcut = 0.0;
      double ecut = 0.0;
    };
  
    struct RigidMolPotentialParameters
    {
      using UserPotParams = std::map< std::pair<std::string,std::string> , RigidMolPotentialPairParam >;

      static constexpr size_t MAX_TYPE_PAIR_IDS = 15; // this lets room for 5 single atom species
      static constexpr size_t MAX_RIGIDMOL_ATOM_TYPES = 6; // single atoms and rigid molecule species together
      static constexpr size_t MAX_ALL_MOLECULE_ATOMS = 4;  // all rigid molecule atoms (position and type for each atom in each rigid molecule with more than 1 atom)
      RigidMolPotentialPairParam m_pair_params[MAX_TYPE_PAIR_IDS]; // indexed by type pair id
      RigidMoleculeAtom m_atoms[MAX_ALL_MOLECULE_ATOMS];
      uint16_t m_n_atoms[MAX_RIGIDMOL_ATOM_TYPES];
      uint16_t m_atoms_start[MAX_RIGIDMOL_ATOM_TYPES];
      
      unsigned int m_nb_pair_params = 0;
      UserPotParams * m_user_pot_parameters = nullptr;
    };
  
    // helper functor, with templated call operator (with/without weights)
    struct RigidMolForceOp
    {
      const RigidMolPotentialParameters* __restrict__ p;

      template<class CPBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        CPBufferT& cpbuf,
        CellParticlesT* cells,
        size_t cell_a,
        size_t p_a,
        exanb::ComputePairParticleContextStart
        ) const
      {
        const auto & p = *(this->p);
        cpbuf.ext.type_a = cells[cell_a][field::type][p_a];
        cpbuf.ext.n_atoms_a = p.m_n_atoms[cpbuf.ext.type_a];
        cpbuf.ext.atoms_start_a = p.m_atoms_start[cpbuf.ext.type_a];
        if( cpbuf.ext.n_atoms_a > 1 )
        {
          const Mat3d mat_a = quaternion_to_matrix( cells[cell_a][field::orient][p_a] );
          for(unsigned int sub_atom_a=0; sub_atom_a < cpbuf.ext.n_atoms_a; ++sub_atom_a)
          {
            cpbuf.ext.mol_atom_ra[sub_atom_a] = mat_a * p.m_atoms[ cpbuf.ext.atoms_start_a + sub_atom_a ].m_pos;
          }
          cpbuf.ext.couple = Vec3d{ 0. , 0. , 0. };
        }
        cpbuf.ext.force = Vec3d{ 0. , 0. , 0. };
        cpbuf.ext.ep = 0.0;
      }

      template<class CPBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        CPBufferT& cpbuf,
        const Vec3d& mol_b_r, double d2,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT cells,size_t cell_b,size_t p_b,
        double weight
        ) const
      {
        const auto & p = *(this->p);

        const unsigned int rm_type_b = cells[cell_b][field::type][p_b];
        const unsigned int n_atoms_b = p.m_n_atoms[rm_type_b];
        const unsigned int atoms_start_b = p.m_atoms_start[rm_type_b];

        if( cpbuf.ext.n_atoms_a > 1 )
        {
          if( n_atoms_b > 1 )
          {
            const auto orient_b = cells[cell_b][field::orient][p_b];
            const Mat3d mat_b = quaternion_to_matrix( orient_b );
            for(unsigned int sub_atom_a=0; sub_atom_a < cpbuf.ext.n_atoms_a ; ++sub_atom_a)
            {
              const unsigned int type_a = p.m_atoms[ cpbuf.ext.atoms_start_a + sub_atom_a ].m_atom_type;
              const Vec3d ra = cpbuf.ext.mol_atom_ra[sub_atom_a];
              Vec3d site_a_F = { 0., 0., 0. };
              double site_a_ep = 0.;
              for(unsigned int sub_atom_b=0; sub_atom_b<n_atoms_b; ++sub_atom_b)
              {
                const unsigned int type_b = p.m_atoms[ atoms_start_b + sub_atom_b ].m_atom_type;
                const unsigned int pair_id = unique_pair_id(type_a,type_b);
                const auto& pp = p.m_pair_params[pair_id]; 
                const Vec3d rb =  mol_b_r + mat_b * p.m_atoms[ atoms_start_b + sub_atom_b ].m_pos;
                const Vec3d dr = rb - ra;
                const double rcut2 = pp.rcut * pp.rcut;
                const double d2 = norm2(dr);
                if( d2 <= rcut2 && weight > 0.0 )
                {
                  const double r = sqrt(d2);
                  double e=0.0, de=0.0;                  
                  USTAMP_POTENTIAL_COMPUTE_WEIGHT( pp.p, pp.pair_params, r, e, de , weight );
                  e -= pp.ecut * weight;
                  de /= r;
                  site_a_F += de * dr;
                  site_a_ep += e * 0.5;
                }
              }
              cpbuf.ext.force  += site_a_F;
              cpbuf.ext.ep     += site_a_ep;
              cpbuf.ext.couple += cross( ra , site_a_F );
            }
          }
          else // particle B is a single atom
          {
            for(unsigned int sub_atom_a=0; sub_atom_a < cpbuf.ext.n_atoms_a ; ++sub_atom_a)
            {
              const unsigned int type_a = p.m_atoms[ cpbuf.ext.atoms_start_a + sub_atom_a ].m_atom_type;
              const Vec3d ra = cpbuf.ext.mol_atom_ra[ sub_atom_a ]; // mat_a * p.m_atoms[sub_atom_a].m_pos;
              const unsigned int pair_id = unique_pair_id(type_a,rm_type_b);
              const auto& pp = p.m_pair_params[pair_id]; 
              const Vec3d dr = mol_b_r - ra;
              const double rcut2 = pp.rcut * pp.rcut;
              const double d2 = norm2(dr);
              if( d2 <= rcut2 && weight > 0.0 )
              {
                const double r = sqrt(d2);
                double e=0.0, de=0.0;                  
                USTAMP_POTENTIAL_COMPUTE_WEIGHT( pp.p, pp.pair_params, r, e, de , weight );
                e -= pp.ecut * weight;
                de /= r;
                const Vec3d site_a_F = de * dr;
                cpbuf.ext.force  += site_a_F;
                cpbuf.ext.couple += cross( ra , site_a_F );
                cpbuf.ext.ep     += e * 0.5;
              }
            }
          }
          
        }
        else // central particle a is a single atom
        {
          if( n_atoms_b > 1 )
          {
            const auto orient_b = cells[cell_b][field::orient][p_b];
            const Mat3d mat_b = quaternion_to_matrix( orient_b );
            for(unsigned int sub_atom_b=0; sub_atom_b<n_atoms_b; ++sub_atom_b)
            {
              const unsigned int type_b = p.m_atoms[ atoms_start_b + sub_atom_b ].m_atom_type;
              const unsigned int pair_id = unique_pair_id(cpbuf.ext.type_a,type_b);
              const auto& pp = p.m_pair_params[pair_id]; 
              const Vec3d dr = mol_b_r + mat_b * p.m_atoms[ atoms_start_b + sub_atom_b ].m_pos;
              const double rcut2 = pp.rcut * pp.rcut;
              const double d2 = norm2(dr);
              if( d2 <= rcut2 && weight > 0.0 )
              {
                const double r = sqrt(d2);
                double e=0.0, de=0.0;                  
                USTAMP_POTENTIAL_COMPUTE_WEIGHT( pp.p, pp.pair_params, r, e, de , weight );
                e -= pp.ecut * weight;
                de /= r;
                cpbuf.ext.force += de * dr;
                cpbuf.ext.ep    += e * 0.5;
              }
            }
          }
          else // single atom pair, just a basic pair potential
          {
            const unsigned int pair_id = unique_pair_id(cpbuf.ext.type_a,rm_type_b);
            const auto& pp = p.m_pair_params[pair_id]; 
            const Vec3d dr = mol_b_r;
            const double rcut2 = pp.rcut * pp.rcut;
            const double d2 = norm2(dr);
            if( d2 <= rcut2 && weight > 0.0 )
            {
              const double r = sqrt(d2);
              double e=0.0, de=0.0;                  
              USTAMP_POTENTIAL_COMPUTE_WEIGHT( pp.p, pp.pair_params, r, e, de , weight );
              e -= pp.ecut * weight;
              de /= r;
              cpbuf.ext.force += de * dr;
              cpbuf.ext.ep    += e * 0.5;
            }
          }
        }
      }

      template<class CPBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        CPBufferT& cpbuf,
        CellParticlesT* cells,
        size_t cell_a,
        size_t p_a,
        exanb::ComputePairParticleContextStop
        ) const
      {
        cells[cell_a][field::fx][p_a]     += cpbuf.ext.force.x;
        cells[cell_a][field::fy][p_a]     += cpbuf.ext.force.y;
        cells[cell_a][field::fz][p_a]     += cpbuf.ext.force.z;
        cells[cell_a][field::ep][p_a]     += cpbuf.ext.ep;
        cells[cell_a][field::couple][p_a] += cpbuf.ext.couple;
      }

    };
  }

}

namespace exanb
{
  template<> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::RigidMolForceOp>
  {
    //static inline constexpr bool HasParticleContextStart = true;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif

    static inline constexpr bool ComputeBufferCompatible = false;
    static inline constexpr bool BufferLessCompatible    = true;

    static inline constexpr bool HasParticleContextStart      = true;    
    static inline constexpr bool HasParticleContext           = true;
    static inline constexpr bool HasParticleContextStop       = true;
  };

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


