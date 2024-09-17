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
  struct AngleForceOp
  {
    CellsT m_cells;
    exanb::spin_mutex_array * m_particle_locks = nullptr;
    const MoleculeGenericFuncParam * const __restrict__ m_func_params = nullptr;;
    const int * const __restrict__ m_angle_param_idx = nullptr;
    const ChemicalAngle * const __restrict__ m_angle_list = nullptr;
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
      
      const int param_idx = m_angle_param_idx[i];
      
      if( param_idx >= 0 )
      {
        const uint64_t atom_to_decode_a = m_angle_list[i][0]; 
        const uint64_t atom_to_decode_b = m_angle_list[i][1];
        const uint64_t atom_to_decode_c = m_angle_list[i][2];
        
        size_t cell_a=0, p_a=0; unsigned int type_a=0;
        size_t cell_b=0, p_b=0; unsigned int type_b=0;
        size_t cell_c=0, p_c=0; unsigned int type_c=0;

        decode_cell_particle(atom_to_decode_a, cell_a, p_a, type_a );
        decode_cell_particle(atom_to_decode_b, cell_b, p_b, type_b );
        decode_cell_particle(atom_to_decode_c, cell_c, p_c, type_c );
        
        const Vec3d ra = Vec3d { m_cells[cell_a][field::rx][p_a] , m_cells[cell_a][field::ry][p_a] , m_cells[cell_a][field::rz][p_a] };
        const Vec3d rb = Vec3d { m_cells[cell_b][field::rx][p_b] , m_cells[cell_b][field::ry][p_b] , m_cells[cell_b][field::rz][p_b] };
        const Vec3d rc = Vec3d { m_cells[cell_c][field::rx][p_c] , m_cells[cell_c][field::ry][p_c] , m_cells[cell_c][field::rz][p_c] };

        // first harm
        Vec3d r1 = periodic_r_delta( rb , ra , m_size_box , m_half_min_size_box );

        // second harm
        Vec3d r2 = periodic_r_delta( rb , rc , m_size_box , m_half_min_size_box ); // rc - rb;
        
        //ldbg << "r1 avant : " << r1 << std::endl;
        if( ! m_xform_is_identity )
        {
          r1 = m_xform * r1;
          r2 = m_xform * r2;
        }

        //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
        const double norm_r1 = norm(r1);
        const double norm_r2 = norm(r2);
        assert(norm_r1>0);
        assert(norm_r2>0);

        //---------------------------------------------------------------------------------------------
        // Compute directions of Fa and Fb
        Vec3d tmp = cross(r1, r2);
        assert( tmp.x!=0. || tmp.y!=0. || tmp.z!=0. );
        const Vec3d dir_a = cross(r1 , tmp);
        const Vec3d dir_b = cross(tmp, r2 );
        const double norm_pa = norm(dir_a);
        const Vec3d pa_r1 = dir_a / norm_pa / norm_r1;
        const double norm_pb = norm(dir_b);
        const Vec3d pb_r2 = dir_b / norm_pb / norm_r2;
        //---------------------------------------------------------------------------------------------

        //---------------------------------------------------------------------------------------------
        // Compute energy
        const double theta = angle(r1,r2);
        const auto fe = intramolecular_func( m_func_params[ param_idx ] , theta ); //it->second->force_energy( theta );
        const double e = fe.second; //b.f_energy(theta);
        //---------------------------------------------------------------------------------------------

        //---------------------------------------------------------------------------------------------
        // Compute forces
        // F = - grad(Ep)
        const double dep_on_dtheta = fe.first; //b.f_forces(theta);

        const Vec3d F1 = pa_r1 * (-dep_on_dtheta);
        const Vec3d F2 = pb_r2 * (-dep_on_dtheta);
        const Vec3d Fo = - ( F1 + F2 );

        if constexpr ( ComputeVirial )
        {
          const Mat3d virial1 = tensor (F1,r1) * 0.5;
          const Mat3d virial2 = tensor (F2,r2) * 0.5;
          const Mat3d virialo = virial1 + virial2;

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_a][p_a]
            , m_cells[cell_a][field::fx][p_a], m_cells[cell_a][field::fy][p_a], m_cells[cell_a][field::fz][p_a], m_cells[cell_a][field::ep][p_a], m_cells[cell_a][m_virial_field][p_a]
            , F1.x, F1.y, F1.z, e/3, virial1 );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_b][p_b]
            , m_cells[cell_b][field::fx][p_b], m_cells[cell_b][field::fy][p_b], m_cells[cell_b][field::fz][p_b], m_cells[cell_b][field::ep][p_b], m_cells[cell_b][m_virial_field][p_b]
            , Fo.x, Fo.y, Fo.z, e/3, virialo );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_c][p_c]
            , m_cells[cell_c][field::fx][p_c], m_cells[cell_c][field::fy][p_c], m_cells[cell_c][field::fz][p_c], m_cells[cell_c][field::ep][p_c], m_cells[cell_c][m_virial_field][p_c]
            , F2.x, F2.y, F2.z, e/3, virial2 );
        }
        else
        {
          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_a][p_a]
            , m_cells[cell_a][field::fx][p_a], m_cells[cell_a][field::fy][p_a], m_cells[cell_a][field::fz][p_a], m_cells[cell_a][field::ep][p_a]
            , F1.x, F1.y, F1.z, e/3 );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_b][p_b]
            , m_cells[cell_b][field::fx][p_b], m_cells[cell_b][field::fy][p_b], m_cells[cell_b][field::fz][p_b], m_cells[cell_b][field::ep][p_b]
            , Fo.x, Fo.y, Fo.z, e/3 );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_c][p_c]
            , m_cells[cell_c][field::fx][p_c], m_cells[cell_c][field::fy][p_c], m_cells[cell_c][field::fz][p_c], m_cells[cell_c][field::ep][p_c]
            , F2.x, F2.y, F2.z, e/3 );
        }      
      }
    }
  };

}


namespace onika
{
  namespace parallel
  {
    template<class CellsT, class VirialFieldT, bool vir> struct ParallelForFunctorTraits< exaStamp::AngleForceOp<CellsT,VirialFieldT,vir> >
    {
      static inline constexpr bool CudaCompatible = true;
    };
  }
}


