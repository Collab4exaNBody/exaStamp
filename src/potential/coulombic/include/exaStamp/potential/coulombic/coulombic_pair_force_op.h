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

// Shared pair force functor and execution helper for coulombic pair operators (wolf, dsf, rf, ewald short range).
// A kernel is a trivially copyable functor : kernel( qi*qj , r , e , de ) with de = de/dr.

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/config.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/core/concurent_add_contributions.h>
#include <onika/log.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <type_traits>

namespace exaStamp
{
  using namespace exanb;

  template<bool _ComputeEnergy, bool _ComputeVirial>
  struct CoulombicPairComputeContext
  {
    double charge_a = 0.0;
    Vec3d f = {0.,0.,0.};
  };

  template<> struct CoulombicPairComputeContext<true,false>
  {
    double charge_a = 0.0;
    Vec3d f = {0.,0.,0.};
    double ep = 0.0;
  };

  template<> struct CoulombicPairComputeContext<true,true>
  {
    double charge_a = 0.0;
    Vec3d f = {0.,0.,0.};
    double ep = 0.0;
    Mat3d virial = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  };

  template<class KernelT, class CPLocksT, class ChargeFieldT, class VirialFieldT, bool _UseSymetry, bool _ComputeEnergy, bool _ComputeVirial>
  struct CoulombicPairForceOp
  {
    static inline constexpr bool ComputeEnergy = _ComputeEnergy;
    static inline constexpr bool ComputeVirial = _ComputeVirial;
    static inline constexpr bool UseSymetry = _UseSymetry;

    static_assert( !ComputeVirial || ComputeEnergy );
    
    const KernelT m_kernel;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    CPLocksT m_locks; // held by value, so that device code never dereferences a host object
    ChargeFieldT m_charge_field;
    VirialFieldT m_virial_field;
    bool m_per_atom_charge = true;

    using ParticleLockT = std::remove_reference_t< decltype( m_locks[0][0] ) >;

    template<class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline double particle_charge( CellParticlesT cells, size_t cell, size_t p ) const
    {
      if( m_per_atom_charge ) return cells[cell][m_charge_field][p];
      else return m_species[ cells[cell][field::type][p] ].m_charge;
    }

    template<class LockT, bool CPAA, bool LOCK, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void add_to_particle_impl( LockT& lock, CellParticlesT cells, size_t cell, size_t p, const Vec3d& f, double e, const Mat3d& vir ) const
    {
      if constexpr ( ComputeEnergy && ComputeVirial )
      {
        concurent_add_contributions<LockT,CPAA,LOCK,double,double,double,double,Mat3d> ( lock
          , cells[cell][field::fx][p], cells[cell][field::fy][p], cells[cell][field::fz][p], cells[cell][field::ep][p], cells[cell][m_virial_field][p]
          , f.x, f.y, f.z, e, vir );
      }
      if constexpr ( ComputeEnergy && !ComputeVirial )
      {
        concurent_add_contributions<LockT,CPAA,LOCK,double,double,double,double> ( lock
          , cells[cell][field::fx][p], cells[cell][field::fy][p], cells[cell][field::fz][p], cells[cell][field::ep][p]
          , f.x, f.y, f.z, e );
      }
      if constexpr ( !ComputeEnergy && !ComputeVirial )
      {
        concurent_add_contributions<LockT,CPAA,LOCK,double,double,double> ( lock
          , cells[cell][field::fx][p], cells[cell][field::fy][p], cells[cell][field::fz][p]
          , f.x, f.y, f.z );
      }
    }

    // adds contributions to a particle. In symmetric mode, a particle may be concurently updated by another thread :
    // atomic adds on GPU, locks on CPU. Device code never touches the (host) lock array.
    template<class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC inline void add_to_particle( CellParticlesT cells, size_t cell, size_t p, const Vec3d& f, double e, const Mat3d& vir ) const
    {
      static constexpr bool CPAA = UseSymetry &&   gpu_device_execution();
      static constexpr bool LOCK = UseSymetry && ! gpu_device_execution() && CPLocksT::use_locks;
      if constexpr ( LOCK )
      {
        add_to_particle_impl<ParticleLockT,CPAA,true>( m_locks[cell][p], cells, cell, p, f, e, vir );
      }
      else
      {
        FakeParticleLock fake_lock;
        add_to_particle_impl<FakeParticleLock,CPAA,false>( fake_lock, cells, cell, p, f, e, vir );
      }
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a , size_t p_a, exanb::ComputePairParticleContextStart ) const
    {
      ctx.ext.f = Vec3d{0.,0.,0.};
      if constexpr ( ComputeEnergy )
      {
        ctx.ext.ep = 0.0;
        if constexpr ( ComputeVirial ) ctx.ext.virial = Mat3d{0.,0.,0.,0.,0.,0.,0.,0.,0.};
      }
      ctx.ext.charge_a = particle_charge( cells, cell_a, p_a );
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () (ComputeBufferT& ctx, CellParticlesT cells, size_t cell_a, size_t p_a, exanb::ComputePairParticleContextStop ) const
    {
      double e = 0.0;
      Mat3d vir = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      if constexpr ( ComputeEnergy ) e = ctx.ext.ep;
      if constexpr ( ComputeVirial ) vir = ctx.ext.virial;
      add_to_particle( cells, cell_a, p_a, ctx.ext.f, e, vir );
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( ComputeBufferT& ctx, Vec3d dr, double d2, CellParticlesT cells, size_t cell_b, size_t p_b, double weight ) const
    {
      const double charge_b = particle_charge( cells, cell_b, p_b );
      const double r = sqrt(d2);
      double e=0.0, de=0.0;
      m_kernel( ctx.ext.charge_a * charge_b, r, e, de );
      e *= weight; de *= weight;
      de /= r;
      const Vec3d dr_fe = de * dr;
      ctx.ext.f += dr_fe;
      Mat3d virial = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      if constexpr ( ComputeEnergy )
      {
        ctx.ext.ep += .5 * e;
        if constexpr ( ComputeVirial )
        {
          virial = tensor( dr_fe, dr ) * -0.5;
          ctx.ext.virial += virial;
        }
      }
      if constexpr ( UseSymetry )
      {
        add_to_particle( cells, cell_b, p_b, -dr_fe, .5*e, virial );
      }
    }
  };

}

namespace exanb
{
  template<class KernelT, class CPLocksT, class ChargeFieldT, class VirialFieldT, bool _UseSymetry, bool _ComputeEnergy, bool _ComputeVirial>
  struct ComputePairTraits< exaStamp::CoulombicPairForceOp<KernelT,CPLocksT,ChargeFieldT,VirialFieldT,_UseSymetry,_ComputeEnergy,_ComputeVirial> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible      = false;
    static inline constexpr bool BufferLessCompatible         = true;
    static inline constexpr bool HasParticleContextStart      = true;    
    static inline constexpr bool HasParticleContext           = true;
    static inline constexpr bool HasParticleContextStop       = true;
    static inline constexpr bool CudaCompatible               = true;
  };
}

namespace exaStamp
{
  struct CoulombicPairOptions
  {
    double rcut = 0.0;
    bool per_atom_charge = true;  // read charges from field::charge, otherwise from species charge
    bool use_symmetry = false;    // neighbor lists are symmetric (half lists), contributions added to both particles
    bool log_energy = false;      // compute per particle energy and virial
    bool ghost_fold_back = false; // symmetric mode only : pairs from owned cells only, ghost particles' contributions are added back by update_(virial_)force_energy_from_ghost
    const CompactGridPairWeights * weights = nullptr;
    GridParticleLocks * particle_locks = nullptr;
  };

  // runs the pair computation for a given kernel. exec_ctx_func() returns a parallel execution context
  template<class GridT, class KernelT, class ExecCtxFuncT>
  inline void coulombic_pair_compute( GridT& grid, const exanb::GridChunkNeighbors& chunk_neighbors, const Domain& domain, const ParticleSpecies& species,
                                      const KernelT& kernel, const CoulombicPairOptions& opt, ExecCtxFuncT exec_ctx_func )
  {
    if( opt.use_symmetry && opt.particle_locks == nullptr )
    {
      fatal_error() << "use_symmetry requires particle_locks" << std::endl;
    }
    if( opt.ghost_fold_back && ! opt.use_symmetry )
    {
      fatal_error() << "ghost_fold_back requires use_symmetry" << std::endl;
    }

    using ChargeFieldT = decltype( grid.field_accessor( field::charge ) );
    using VirialFieldT = decltype( grid.field_accessor( field::virial ) );
    const ChargeFieldT charge_field = grid.field_accessor( field::charge );
    const VirialFieldT virial_field = grid.field_accessor( field::virial );

    auto run = [&]( auto cp_locks , auto cp_weight , auto use_sym_tag , auto energy_tag )
    {
      static constexpr bool UseSym = decltype(use_sym_tag)::value;
      static constexpr bool Energy = decltype(energy_tag)::value;
      using ForceOp = CoulombicPairForceOp<KernelT,decltype(cp_locks),ChargeFieldT,VirialFieldT,UseSym,Energy,Energy>;
      using CPBufT = ComputePairBuffer2<false,false, CoulombicPairComputeContext<Energy,Energy> >;
      exanb::GridChunkNeighborsLightWeightIt<UseSym> nbh_it{ chunk_neighbors };
      LinearXForm cp_xform { domain.xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
      static constexpr std::true_type use_cells_accessor = {};
      ForceOp force_op { kernel, species.data(), cp_locks, charge_field, virial_field, opt.per_atom_charge };
      compute_cell_particle_pairs2( grid, opt.rcut, UseSym && ! opt.ghost_fold_back, optional, make_compute_pair_buffer<CPBufT>(), force_op, onika::FlatTuple<>{}, DefaultPositionFields{}, exec_ctx_func(), use_cells_accessor );
    };

    auto run_opt_energy = [&]( auto cp_locks , auto cp_weight , auto use_sym_tag )
    {
      if( opt.log_energy ) run( cp_locks, cp_weight, use_sym_tag, std::true_type{} );
      else                 run( cp_locks, cp_weight, use_sym_tag, std::false_type{} );
    };

    auto run_opt_weights = [&]( auto cp_locks , auto use_sym_tag )
    {
      if( opt.weights != nullptr ) run_opt_energy( cp_locks, CompactPairWeightIterator{ opt.weights->m_cell_weights.data() }, use_sym_tag );
      else                         run_opt_energy( cp_locks, ComputePairNullWeightIterator{}, use_sym_tag );
    };

    if( opt.use_symmetry ) run_opt_weights( ComputePairOptionalLocks<true>{ opt.particle_locks->data() } , std::true_type{} );
    else                   run_opt_weights( ComputePairOptionalLocks<false>{} , std::false_type{} );
  }

}
