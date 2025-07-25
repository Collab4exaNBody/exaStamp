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
  struct ImproperForceOp
  {
    CellsT m_cells;
    exanb::spin_mutex_array * m_particle_locks = nullptr;
    const MoleculeGenericFuncParam * const __restrict__ m_func_params = nullptr;;
    const int * const __restrict__ m_improper_param_idx = nullptr;
    const ChemicalImproper * const __restrict__ m_improper_list = nullptr;
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
      
      const int param_idx = m_improper_param_idx[i];
      if( param_idx >= 0 )
      {    
        const uint64_t atom_to_decode_a = m_improper_list[i][0]; 
        const uint64_t atom_to_decode_b = m_improper_list[i][1];
        const uint64_t atom_to_decode_c = m_improper_list[i][2];
        const uint64_t atom_to_decode_d = m_improper_list[i][3];
        
        size_t cell_a=0, p_a=0; unsigned int type_a=0;
        size_t cell_b=0, p_b=0; unsigned int type_b=0;
        size_t cell_c=0, p_c=0; unsigned int type_c=0;
        size_t cell_d=0, p_d=0; unsigned int type_d=0;

        decode_cell_particle(atom_to_decode_a, cell_a, p_a, type_a );
        decode_cell_particle(atom_to_decode_b, cell_b, p_b, type_b );
        decode_cell_particle(atom_to_decode_c, cell_c, p_c, type_c );
        decode_cell_particle(atom_to_decode_d, cell_d, p_d, type_d );
        
        const Vec3d r0 = Vec3d { m_cells[cell_a][field::rx][p_a] , m_cells[cell_a][field::ry][p_a] , m_cells[cell_a][field::rz][p_a] };
        const Vec3d r1 = Vec3d { m_cells[cell_b][field::rx][p_b] , m_cells[cell_b][field::ry][p_b] , m_cells[cell_b][field::rz][p_b] };
        const Vec3d r2 = Vec3d { m_cells[cell_c][field::rx][p_c] , m_cells[cell_c][field::ry][p_c] , m_cells[cell_c][field::rz][p_c] };
        const Vec3d r3 = Vec3d { m_cells[cell_d][field::rx][p_d] , m_cells[cell_d][field::ry][p_d] , m_cells[cell_d][field::rz][p_d] };

        Vec3d r10 = periodic_r_delta( r1 , r0 , m_size_box , m_half_min_size_box );
        Vec3d r12 = periodic_r_delta( r1 , r2 , m_size_box , m_half_min_size_box );
        Vec3d r23 = periodic_r_delta( r2 , r3 , m_size_box , m_half_min_size_box );

        if( ! m_xform_is_identity )
        {
          r10 = m_xform * r10;
          r12 = m_xform * r12;
          r23 = m_xform * r23;
        }

        //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
        // Compute energy
        const Vec3d V = cross(r10,r12);
        const Vec3d W = cross(r23,r12);

        assert(V!=(Vec3d{0,0,0}));//case r10 and r12 aligned
        assert(W!=(Vec3d{0,0,0}));//case r12 and r23 aligned

        const double phi = signum(dot(cross(V,W), r12)) * angle(V, W);
        const auto force_energy = intramolecular_func( m_func_params[ param_idx ] , phi );
        const double e = force_energy.second / 4.; //t.f_energy(phi);
        //---------------------------------------------------------------------------------------------

        const double norm_r2 = norm(r12);
        assert(norm_r2>0);

        //pa = V/norm(V)
        const double sqr_norm_V = dot(V,V);
        const Vec3d pa_r1sintheta1 = V * norm_r2 / sqr_norm_V;

        //pd = m_W/norm(m_W)
        const double sqr_norm_W = dot(W,W);
        const Vec3d pd_r3sintheta2 = - W * norm_r2 / sqr_norm_W ;

        //---------------------------------------------------------------------------------------------
        // Compute forces
        // F = - grad(Ep)
        const double dep_on_dphi = force_energy.first; //t.f_forces(phi);

        //Force on A
        const Vec3d Fa = -dep_on_dphi * pa_r1sintheta1;
        
        //Force on D
        const Vec3d Fd = -dep_on_dphi * pd_r3sintheta2;

        //Force on C
        //-1/norm(r12)^2 ((r23 + r12) ^ FD + r1 ^ FA) ^ r12    
        const double sqr_norm_r2 = dot(r12,r12);
        const Vec3d  Fc = - cross( cross(r23+r12,Fd) + cross(r10,Fa) , r12 ) / sqr_norm_r2;

        //e /= 4.;
        const Vec3d Fo = - (Fa+Fc+Fd);

        if constexpr ( ComputeVirial )
        {
          const Mat3d vira = tensor (Fa,r10) * 0.5;
          const Mat3d virc = tensor (Fc,r12) * 0.5;
          const Mat3d vird = tensor (Fd,r23) * 0.5;      
          const Mat3d virac = vira + virc;
          const Mat3d vircd = virc + vird;

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_a][p_a]
            , m_cells[cell_a][field::fx][p_a], m_cells[cell_a][field::fy][p_a], m_cells[cell_a][field::fz][p_a], m_cells[cell_a][field::ep][p_a], m_cells[cell_a][m_virial_field][p_a]
            , Fa.x, Fa.y, Fa.z, e, vira );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_b][p_b]
            , m_cells[cell_b][field::fx][p_b], m_cells[cell_b][field::fy][p_b], m_cells[cell_b][field::fz][p_b], m_cells[cell_b][field::ep][p_b], m_cells[cell_b][m_virial_field][p_b]
            , Fo.x, Fo.y, Fo.z, e, virac );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_c][p_c]
            , m_cells[cell_c][field::fx][p_c], m_cells[cell_c][field::fy][p_c], m_cells[cell_c][field::fz][p_c], m_cells[cell_c][field::ep][p_c], m_cells[cell_c][m_virial_field][p_c]
            , Fc.x, Fc.y, Fc.z, e, vircd );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double,Mat3d> (
              cp_locks[cell_d][p_d]
            , m_cells[cell_d][field::fx][p_d], m_cells[cell_d][field::fy][p_d], m_cells[cell_d][field::fz][p_d], m_cells[cell_d][field::ep][p_d], m_cells[cell_d][m_virial_field][p_d]
            , Fd.x, Fd.y, Fd.z, e, vird );
        }
        else
        {
          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_a][p_a]
            , m_cells[cell_a][field::fx][p_a], m_cells[cell_a][field::fy][p_a], m_cells[cell_a][field::fz][p_a], m_cells[cell_a][field::ep][p_a]
            , Fa.x, Fa.y, Fa.z, e );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_b][p_b]
            , m_cells[cell_b][field::fx][p_b], m_cells[cell_b][field::fy][p_b], m_cells[cell_b][field::fz][p_b], m_cells[cell_b][field::ep][p_b]
            , Fo.x, Fo.y, Fo.z, e );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_c][p_c]
            , m_cells[cell_c][field::fx][p_c], m_cells[cell_c][field::fy][p_c], m_cells[cell_c][field::fz][p_c], m_cells[cell_c][field::ep][p_c]
            , Fc.x, Fc.y, Fc.z, e );

          concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,double> (
              cp_locks[cell_d][p_d]
            , m_cells[cell_d][field::fx][p_d], m_cells[cell_d][field::fy][p_d], m_cells[cell_d][field::fz][p_d], m_cells[cell_d][field::ep][p_d]
            , Fd.x, Fd.y, Fd.z, e );
        }
      }
    }
  };

}


namespace onika
{
  namespace parallel
  {
    template<class CellsT, class VirialFieldT, bool vir> struct ParallelForFunctorTraits< exaStamp::ImproperForceOp<CellsT,VirialFieldT,vir> >
    {
      static inline constexpr bool CudaCompatible = true;
    };
  }
}


