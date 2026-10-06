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

#include <exanb/core/grid.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exaStamp/potential/coulombic/coulombic_pair_force_op.h>
#include <exaStamp/potential/coulombic/ewald.h>

namespace exaStamp
{
inline namespace coulombic_ewald // distinct symbols from the legacy ewald plugin (plugins are loaded RTLD_GLOBAL)
{
  using namespace exanb;

  struct EwaldShortRangeKernel
  {
    ReadOnlyEwaldParameters m_params;
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double c, double r, double& e, double& de) const { ewald_compute_energy( m_params, c, r, e, de ); }
  };

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class EwaldShortRangePC : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( EwaldParameters           , ewald_config        , INPUT , REQUIRED );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors     , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( CompactGridPairWeights    , compact_nbh_weight  , INPUT , OPTIONAL );
    ADD_SLOT( bool                      , enable_pair_weights , INPUT , true );
    ADD_SLOT( bool                      , per_atom_charge     , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( bool                      , use_symmetry        , INPUT , false , DocString{"must match the symmetric setting of neighbor lists"} );
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
      opt.rcut = ewald_config->radius;
      *rcut_max = std::max( *rcut_max , opt.rcut );

      // nothing to compute : usefull when compute_force is called at the very first to initialize rcut_max
      if( grid->number_of_cells() == 0 ) return ;

      opt.log_energy = trigger_thermo_state.has_value() ? *trigger_thermo_state : false;
      opt.per_atom_charge = *per_atom_charge;
      opt.use_symmetry = *use_symmetry;
      if( compact_nbh_weight.has_value() && *enable_pair_weights ) opt.weights = compact_nbh_weight.get_pointer();
      if( particle_locks.has_value() ) opt.particle_locks = particle_locks.get_pointer();

      ldbg << std::boolalpha << "coulombic_ewald_short_range: rc="<< opt.rcut <<" , g_ewald="<< ewald_config->g_ewald <<" , pair_weights="<< (opt.weights!=nullptr)
           <<" , log_energy="<< opt.log_energy <<" , use_symmetry="<< opt.use_symmetry <<" , per_atom_charge="<< opt.per_atom_charge << std::endl;

      const EwaldShortRangeKernel kernel = { ReadOnlyEwaldParameters( *ewald_config ) };
      coulombic_pair_compute( *grid, *chunk_neighbors, *domain, *species, kernel, opt, [this](){ return parallel_execution_context(); } );
    }
  };

  template<class GridT> using EwaldShortRangePCTmpl = EwaldShortRangePC<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(coulombic_ewald_short_range)
  {
    OperatorNodeFactory::instance()->register_factory( "coulombic_ewald_short_range" , make_grid_variant_operator<EwaldShortRangePCTmpl> );
  }

}
}
