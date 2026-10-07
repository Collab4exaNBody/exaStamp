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
// and for the one-body self energy term (wolf, dsf).

#include <exaStamp/potential/coulombic/coulombic_pair_force_op.h>
#include <exanb/compute/compute_cell_particles.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/flat_tuple.h>
#include <exaStamp/unit_system.h>

namespace exaStamp
{
  using namespace exanb;

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
    ADD_SLOT( bool                      , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( Domain                    , domain              , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies           , species             , INPUT , REQUIRED );    
    ADD_SLOT( GridT                     , grid                , INPUT_OUTPUT );
    ADD_SLOT( double                    , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( GridParticleLocks         , particle_locks      , INPUT_OUTPUT , OPTIONAL , DocString{"particle spin locks"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      CoulombicPairOptions opt;
      opt.rcut = parameters->rc;
      *rcut_max = std::max( *rcut_max , opt.rcut );
      
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
    }
  };

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
  // one-body self energy term (energy only, no force), computed when thermodynamic state is requested
  template< class GridT, class ParamsT, class KernelT, class = AssertGridHasFields< GridT, field::_ep > >
  class CoulombicSelfPC : public OperatorNode
  {
    ADD_SLOT( GridT             , grid                , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies   , species             , INPUT , REQUIRED );
    ADD_SLOT( ParamsT           , parameters          , INPUT , REQUIRED );
    ADD_SLOT( bool              , per_atom_charge     , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( bool              , trigger_thermo_state, INPUT , OPTIONAL );

  public:
    inline void execute () override final
    {
      const bool log_energy = trigger_thermo_state.has_value() ? *trigger_thermo_state : false;
      if( ! log_energy || grid->number_of_cells() == 0 ) return;
      auto ep = grid->field_accessor( field::ep );
      if( *per_atom_charge )
      {
        CoulombicSelfEnergyFunc<KernelT,true> func = { KernelT{ *parameters } , species->data() };
        compute_cell_particles( *grid , false , func , onika::make_flat_tuple( ep, grid->field_accessor( field::charge ) ) , parallel_execution_context() );
      }
      else
      {
        CoulombicSelfEnergyFunc<KernelT,false> func = { KernelT{ *parameters } , species->data() };
        compute_cell_particles( *grid , false , func , onika::make_flat_tuple( ep, grid->field_accessor( field::type ) ) , parallel_execution_context() );
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Self energy term -(e_shift/2 + alpha/sqrt(pi)).q^2/(4.pi.epsilon0) of damped shifted coulombic potentials (wolf, dsf).
Must be used together with the corresponding pair operator (coulombic_wolf, coulombic_dsf, coul_wolf_pair...), which does not include it.
)EOF";
    }
  };

}
