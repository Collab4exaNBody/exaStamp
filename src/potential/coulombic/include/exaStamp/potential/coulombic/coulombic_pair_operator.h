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

// Generic operators for coulombic pair potentials whose parameters have a cutoff member 'rc' (wolf, dsf, rf),
// and for the one-body self energy term (wolf, dsf). Kernels provide operator()(c,r,e,de), a static documentation()
// and, when the method has a self energy, self_energy(q).

#include <exaStamp/potential/coulombic/coulombic_pair_force_op.h>
#include <exanb/compute/compute_cell_particles.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/flat_tuple.h>
#include <exaStamp/unit_system.h>

namespace exaStamp
{
  using namespace exanb;

  // self energy functor : ep += kernel.self_energy(q)
  template<class KernelT, bool PerAtomCharge>
  struct CoulombicSelfEnergyFunc
  {
    const KernelT m_kernel;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double& ep, ChargeOrTypeT ct) const
    {
      double q = 0.0;
      if constexpr ( PerAtomCharge ) q = ct;
      else q = m_species[ct].m_charge;
      ep += m_kernel.self_energy( q );
    }
  };
}

namespace exanb
{
  template<class KernelT, bool PerAtomCharge> struct ComputeCellParticlesTraits< exaStamp::CoulombicSelfEnergyFunc<KernelT,PerAtomCharge> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  template<class KernelT> concept CoulombicKernelWithSelfEnergy = requires( const KernelT& k ) { k.self_energy( 0.0 ); };

  // adds the one-body self energy of owned particles to ep (energy only, no force)
  template<class GridT, class KernelT, class ExecCtxFuncT>
  inline void coulombic_self_energy_compute( GridT& grid, const ParticleSpecies& species, const KernelT& kernel, bool per_atom_charge, ExecCtxFuncT exec_ctx )
  {
    auto ep = grid.field_accessor( field::ep );
    if( per_atom_charge )
    {
      CoulombicSelfEnergyFunc<KernelT,true> func = { kernel , species.data() };
      compute_cell_particles( grid , false , func , onika::make_flat_tuple( ep, grid.field_accessor( field::charge ) ) , exec_ctx() );
    }
    else
    {
      CoulombicSelfEnergyFunc<KernelT,false> func = { kernel , species.data() };
      compute_cell_particles( grid , false , func , onika::make_flat_tuple( ep, grid.field_accessor( field::type ) ) , exec_ctx() );
    }
  }
}

namespace exaStamp
{
  template< class GridT, class ParamsT, class KernelT, class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz > >
  class CoulombicPairPC : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( ParamsT                   , parameters          , INPUT , REQUIRED );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors     , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( CompactGridPairWeights    , compact_nbh_weight  , INPUT , OPTIONAL );
    ADD_SLOT( bool                      , enable_pair_weights , INPUT , true );
    ADD_SLOT( bool                      , per_atom_charge     , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( bool                      , use_symmetry        , INPUT , false , DocString{"must match the symmetric setting of neighbor lists"} );
    ADD_SLOT( bool                      , ghost_fold_back     , INPUT , false , DocString{"with use_symmetry : compute pairs from owned cells only (faster) ; ghost contributions must be added back by update_virial_force_energy_from_ghost in compute_force_epilog, after zero_force_energy: { ghost: true } in compute_force_prolog"} );
    ADD_SLOT( bool                      , self_energy         , INPUT , true , DocString{"add the one-body self energy term (wolf, dsf ; ignored for rf, which has none) to per particle energies"} );
    ADD_SLOT( bool                      , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( Domain                    , domain              , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies           , species             , INPUT , REQUIRED );    
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( bool                      , coulombic_self_energy_included , OUTPUT , DocString{"true when this operator computes the self energy term (read by coulombic_*_self to prevent double counting)"} );
    ADD_SLOT( double                    , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( GridParticleLocks         , particle_locks      , INPUT_OUTPUT , OPTIONAL , DocString{"particle spin locks"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      CoulombicPairOptions opt;
      opt.rcut = parameters->rc;
      *rcut_max = std::max( *rcut_max , opt.rcut );
      const bool with_self_energy = CoulombicKernelWithSelfEnergy<KernelT> && *self_energy;
      *coulombic_self_energy_included = with_self_energy;

      // nothing to compute : usefull when compute_force is called at the very first to initialize rcut_max
      if( grid->number_of_cells() == 0 ) return ;

      opt.log_energy = trigger_thermo_state.has_value() ? *trigger_thermo_state : false;
      opt.per_atom_charge = *per_atom_charge;
      opt.use_symmetry = *use_symmetry;
      opt.ghost_fold_back = *ghost_fold_back;
      if( compact_nbh_weight.has_value() && *enable_pair_weights ) opt.weights = compact_nbh_weight.get_pointer();
      if( particle_locks.has_value() ) opt.particle_locks = particle_locks.get_pointer();

      ldbg << std::boolalpha << name() << ": rc="<< opt.rcut <<" , pair_weights="<< (opt.weights!=nullptr) <<" , log_energy="<< opt.log_energy
           <<" , use_symmetry="<< opt.use_symmetry <<" , ghost_fold_back="<< opt.ghost_fold_back <<" , per_atom_charge="<< opt.per_atom_charge << std::endl;

      coulombic_pair_compute( *grid, *chunk_neighbors, *domain, *species, KernelT{ *parameters }, opt, [this](){ return parallel_execution_context(); } );

      if constexpr ( CoulombicKernelWithSelfEnergy<KernelT> )
      {
        if( with_self_energy && opt.log_energy )
        {
          coulombic_self_energy_compute( *grid, *species, KernelT{ *parameters }, opt.per_atom_charge, [this](){ return parallel_execution_context("self_energy"); } );
        }
      }
    }

    inline std::string documentation() const override final
    {
      return std::string( KernelT::documentation() ) + R"EOF(
Charges are read from the per particle charge field (per_atom_charge: true, default) or from the species.
Forces are always computed ; per particle energies and virial when trigger_thermo_state is true.
With use_symmetry: true, ghost_fold_back: true computes pairs from owned cells only (see the slot documentation).
)EOF";
    }
  };
}

namespace exaStamp
{
  // one-body self energy term (energy only, no force), computed when thermodynamic state is requested
  template< class GridT, class ParamsT, class KernelT, class = AssertGridHasFields< GridT, field::_ep > >
  class CoulombicSelfPC : public OperatorNode
  {
    ADD_SLOT( GridT             , grid                , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies   , species             , INPUT , REQUIRED );
    ADD_SLOT( ParamsT           , parameters          , INPUT , REQUIRED );
    ADD_SLOT( bool              , per_atom_charge     , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( bool              , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( bool              , coulombic_self_energy_included , INPUT , OPTIONAL , DocString{"set by coulombic_wolf / coulombic_dsf when they already compute the self energy"} );

  public:
    inline void execute () override final
    {
      if( coulombic_self_energy_included.has_value() && *coulombic_self_energy_included )
      {
        fatal_error() << name() << " : the self energy is already computed by the coulombic pair operator. "
                      << "Remove "<< name() <<", or set self_energy: false on the pair operator." << std::endl;
      }
      const bool log_energy = trigger_thermo_state.has_value() ? *trigger_thermo_state : false;
      if( ! log_energy || grid->number_of_cells() == 0 ) return;
      coulombic_self_energy_compute( *grid, *species, KernelT{ *parameters }, *per_atom_charge, [this](){ return parallel_execution_context(); } );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Self energy term -(e_shift/2 + alpha/sqrt(pi)).q^2/(4.pi.epsilon0) of damped shifted coulombic potentials (wolf, dsf).
Energy only (no force), computed when trigger_thermo_state is true. For the pair styles (coul_wolf, coul_dsf, ljwolf),
which do not include it, with per_atom_charge: false. coulombic_wolf and coulombic_dsf already include it (self_energy slot).
)EOF";
    }
  };

}
