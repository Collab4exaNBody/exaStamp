#pragma once

#include <cmath>
#include <functional>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exaStamp/molecule/potential_functional.h>
#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exanb/core/concurent_add_contributions.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <onika/integral_constant.h>
#include <onika/parallel/parallel_for.h>

namespace exaStamp
{

  template<class CellsT, class VirialFieldT, bool ComputeVirial = false>
  struct BondForceOp
  {
    CellsT m_cells;
    exanb::spin_mutex_array * m_particle_locks = nullptr;
    const MoleculeGenericFuncParam * const __restrict__ m_func_params = nullptr;;
    const int * const __restrict__ m_bond_param_idx = nullptr;
    const ChemicalBond * const __restrict__ m_bond_list = nullptr;
    const Vec3d m_size_box = {0.,0.,0.};
    const double m_half_min_size_box = 0.;
    const Mat3d m_xform = {};
    const bool m_xform_is_identity = false;
    VirialFieldT m_virial_field;
    
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t i ) const
    {
      static constexpr bool CPAA = gpu_device_execution();
      static constexpr bool LOCK = ! gpu_device_execution();
      exanb::ComputePairOptionalLocks<LOCK> cp_locks = {};
      if constexpr ( LOCK ) cp_locks.m_locks = m_particle_locks;
      using ParticleLockT = decltype( cp_locks[0][0] );
      
      const int param_idx = m_bond_param_idx[i];
      
      if( param_idx >= 0 )
      {    
        const uint64_t atom_to_decode_a = m_bond_list[i][0]; //atom_from_idmap( bonds_list[i][0] ,id_map, id_map_ghosts ); 
        const uint64_t atom_to_decode_b = m_bond_list[i][1]; //atom_from_idmap( bonds_list[i][1] ,id_map, id_map_ghosts ); 
        
        size_t cell_a=0, p_a=0, cell_b=0, p_b=0;
        unsigned int type_a=0, type_b=0;

        decode_cell_particle(atom_to_decode_a, cell_a, p_a, type_a );
        decode_cell_particle(atom_to_decode_b, cell_b, p_b, type_b );
        const Vec3d ra = Vec3d { m_cells[cell_a][field::rx][p_a] , m_cells[cell_a][field::ry][p_a] , m_cells[cell_a][field::rz][p_a] };
        const Vec3d rb = Vec3d { m_cells[cell_b][field::rx][p_b] , m_cells[cell_b][field::ry][p_b] , m_cells[cell_b][field::rz][p_b] };
        Vec3d r = periodic_r_delta( ra , rb , m_size_box , m_half_min_size_box );

        if( ! m_xform_is_identity ) { r = m_xform * r; }

        double norm_r = norm(r);
        assert(norm_r>0);

        // compute force and energy
        const auto fe = intramolecular_func( m_func_params[ param_idx ] , norm_r );

        double e = fe.second;

        // Compute forces
        // F = dEp/dr nij
        double dEp_dr = fe.first;
        Vec3d F = ( r * dEp_dr ) / norm_r;

        if constexpr ( ComputeVirial )
        {
          auto virial = tensor(F,r) * 0.5;
          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_a][p_a]
            , m_cells[cell_a][field::fx][p_a], m_cells[cell_a][field::fy][p_a], m_cells[cell_a][field::fz][p_a], m_cells[cell_a][field::ep][p_a], m_cells[cell_a][m_virial_field][p_a]
            , F.x, F.y, F.z, e*0.5, virial );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_b][p_b]
            , m_cells[cell_b][field::fx][p_b], m_cells[cell_b][field::fy][p_b], m_cells[cell_b][field::fz][p_b], m_cells[cell_b][field::ep][p_b], m_cells[cell_b][m_virial_field][p_b]
            , -F.x, -F.y, -F.z, e*0.5, virial );
        }
        else
        {
          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_a][p_a]
            , m_cells[cell_a][field::fx][p_a], m_cells[cell_a][field::fy][p_a], m_cells[cell_a][field::fz][p_a], m_cells[cell_a][field::ep][p_a]
            , F.x, F.y, F.z, e*0.5 );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_b][p_b]
            , m_cells[cell_b][field::fx][p_b], m_cells[cell_b][field::fy][p_b], m_cells[cell_b][field::fz][p_b], m_cells[cell_b][field::ep][p_b]
            , -F.x, -F.y, -F.z, e*0.5 );
        }
      }
    }
  };

}


namespace onika
{
  namespace parallel
  {
    template<class CellsT, class VirialFieldT, bool vir> struct ParallelForFunctorTraits< exaStamp::BondForceOp<CellsT,VirialFieldT,vir> >
    {
      static inline constexpr bool CudaCompatible = true;
    };
  }
}

